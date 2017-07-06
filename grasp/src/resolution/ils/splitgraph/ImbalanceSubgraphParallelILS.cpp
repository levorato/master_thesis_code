/*
 * ImbalanceSubgraphParallelILS.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#include "./include/ImbalanceSubgraphParallelILS.h"
#include "util/include/ProcessUtil.h"
#include "util/parallel/include/MPIUtil.h"
#include "../include/CUDAILS.h"
#include "problem/include/RCCProblem.h"

#include "../../construction/include/GainFunctionFactory.h"
#include "../../construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "util/include/RandomUtil.h"
#include "problem/include/ClusteringProblemFactory.h"

#include "../../validation/include/CCImbalanceCalculator.h"

#include <iomanip>
#include <cstring>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <cfloat>
#include <algorithm>
#include <boost/math/special_functions/round.hpp>
// Enable PBGL interfaces to BGL algorithms
#include <boost/graph/use_mpi.hpp>
// Communicate via MPI
#include <boost/graph/distributed/mpi_process_group.hpp>

// Distributed adjacency list
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/property_map/property_map_iterator.hpp>
#include <boost/graph/distributed/vertex_list_adaptor.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/distributed/distributed_graph_utility.hpp>
#include <boost/graph/distributed/vertex_list_adaptor.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/property_map_iterator.hpp>
#include <boost/graph/distributed/vertex_list_adaptor.hpp>
#include <boost/graph/iteration_macros.hpp>
// #include <boost/graph/distributed/adjlist/redistribute.hpp>


namespace ublas = boost::numeric::ublas::detail;

namespace resolution {
namespace ils {

using namespace resolution::construction;
using namespace boost::algorithm;
using namespace std;
using namespace util;
using namespace util::parallel;

//#define is_overloaded_process(list) ((list.size() >= 2*(long)ceil(g->getGlobalN() / (double)numberOfProcesses)) and (g->getN() <= 110000))
#define is_overloaded_process(list) (true)

ImbalanceSubgraphParallelILS::ImbalanceSubgraphParallelILS(const unsigned long& seed, const int& allocationStrategy, const int& slaves, const int& searchSlaves,
		const bool& split, const bool& cuda, const bool& pgraph) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split), cudaEnabled(cuda), parallelgraph(pgraph), vertexImbalance(), timeSpentAtIteration(), util(),
		numberOfFrustratedSolutions(0), randomSeed(seed) {

}

ImbalanceSubgraphParallelILS::~ImbalanceSubgraphParallelILS() {
	// TODO Auto-generated destructor stub
}

clusteringgraph::Clustering ImbalanceSubgraphParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		clusteringgraph::SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info, std::vector<long>& verticesInLeaderProcess) {

	stringstream constructivePhaseResults;
	stringstream iterationResults;
	stringstream iterationTimeSpent;
	unsigned int numberOfProcesses = numberOfSlaves + 1;
	this->verticesInLeaderProcess = verticesInLeaderProcess;

	BOOST_LOG_TRIVIAL(info) << "Starting split graph ILS (local ILS time limit = " << LOCAL_ILS_TIME_LIMIT << " s)...";

	// *** STEP A => PREPROCESSING PHASE OF DISTRIBUTED METAHEURISTIC: New split graph partitioning here
	// Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// TODO IMPLEMENTAR PARTICIONAMENTO INICIAL PARA A PARALLEL BGL
	ClusterArray splitgraphClusterArray(g->getGlobalN());
	// preenchimento inicial baseado na distribuicao da parallel bgl

	ProcessClustering tmpsplitgraphClustering(g, problem, splitgraphClusterArray); // = preProcessSplitgraphPartitioning(g, problem, true);

	// *** STEP B => Calls the individual ILS processing for each subgraph, invoking Parallel ILS with MPI
	// WARNING: THE SPLITGRAPHCLUSTERARRAY IS A SEPARATE DATA STRUCTURE, DIFFERENT THAN THE CURRENT CLUSTERING
	// IT CONTROLS THE PARTITIONING BETWEEN PROCESSORS (SPLIT GRAPH)
	clusteringgraph::Clustering Cc(splitgraphClusterArray, *g, problem, 0.0, 0.0);
	Cc = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info,
			tmpsplitgraphClustering, Cc);

	// FULL CALCULATION OF PROCESS IMBALANCE MATRIZ, BROUGHT FROM CONTRUCTIVE PHASE
	ClusterArray splitClusterArray = tmpsplitgraphClustering.getClusterArray();
	ImbalanceMatrix processClusterImbMatrix = util.calculateProcessToProcessImbalanceMatrix(*g, splitClusterArray,
					this->vertexImbalance, numberOfSlaves + 1);
	ProcessClustering splitgraphClustering(g, problem, splitClusterArray, processClusterImbMatrix);
	ClusterArray initialSplitgraphClusterArray = splitgraphClustering.getClusterArray();

	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timeSpentInILS = timeSpentInConstruction;
	timer.resume();
	start_time = timer.elapsed();

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Initial solution found: I(P) = " <<
			Cc.getImbalance().getValue();
	constructivePhaseResults << 1 << "," << Cc.getImbalance().getValue() << ","
					<< Cc.getImbalance().getPositiveValue()
					<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << Cc.getImbalance().getValue()
						<< "\n";
	Cc.printClustering(g->getGlobalN());
	clusteringgraph::Clustering bestClustering(Cc);
	ProcessClustering bestSplitgraphClustering(splitgraphClustering);

	// STEP 2: TRY TO IMPROVE THE GLOBAL SOLUTION THROUGH A DISTRIBUTED METAHEURISTIC WHICH
	// EXCHANGES VERTICES BETWEEN PROCESSES
	iterationResults << "VNDCount,BestSol,BestSol+,BestSol-,k,Time,NumImprovements,NumFrustrations\n";
	int improvementCount = 0, improvements = 0, executionCount = 0;
	long totalFrustratedSolutions = 0;
	do {
		// 1. O processo mestre será o responsável por manter duas estruturas de dados de controle do imbalance:
		//		Um vetor com os vértices que mais contribuem para I(P) (soma de imbalance de cada vertice)
		//		Uma matriz com a soma do imbalance entre processos

		// Distributed VND
		improvements = variableNeighborhoodDescent(g, bestSplitgraphClustering,
				bestClustering, numberOfProcesses,
				construct, vnd,
				iter, iterMaxILS, perturbationLevelMax,
				problem, info, timeSpentInILS, executionCount);
		executionCount++;
		improvementCount += improvements;
		totalFrustratedSolutions += numberOfFrustratedSolutions;

		// Last step. Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();

		// Write the results into ostream os, using csv format
		// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		// iteration results are only being shown when the solution improved
		iterationResults << (executionCount) << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
				<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
				<< "," << fixed << setprecision(4) << timeSpentInILS
				<< "," << improvements << "," << numberOfFrustratedSolutions << "\n";
		timer.resume();
		start_time = timer.elapsed();
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInILS >= vnd->getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}
	} while(improvements > 0 and executionCount < 1);

	// prints all the info about time spent in each distributed MH invocation
	for(int it = 0; it < timeSpentAtIteration.size(); it++) {
		for(int px = 0; px < numberOfProcesses; px++) {
			iterationTimeSpent << timeSpentAtIteration[it][px] << ", ";
		}
		iterationTimeSpent << "\n";
	}

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of improvements to initial solution: " << improvementCount;
	// *** ONLY IN THE END OF THE DISTRIBUTED METAHEURISTIC: runs a ILS iteration (or a local search?) over the merged global solution

	// 6. Runs a local search over the merged global solution, trying to improve it
	// TODO Run a full ILS iteration over the split graph initial solution
	/*  DISABLED FOR THE DISTRIBUTED GRAPH, FOR OBVIOUS REASONS
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Applying local search over the initial solution.";
	NeighborhoodSearch* neigborhoodSearch;
	GainFunctionFactory functionFactory(g);
	NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, this->numberOfSlaves, this->numberOfSearchSlaves);
	neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	VariableNeighborhoodDescent vnd2(*neigborhoodSearch, vnd->getRandomSeed(), vnd->getNeighborhoodSize(),
			false, vnd->getTimeLimit());
	bestClustering = vnd2.localSearch(g, bestClustering, 0, problem, timeSpentInILS, info.processRank);
	CBest = bestClustering;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after applying local search: I(P) = " << bestClustering.getImbalance().getValue();

	ConstructClustering NoConstruct(functionFactory.build(construct->getGainFunctionType()), vnd->getRandomSeed(),
			construct->getAlpha(), &bestClustering);
	CUDAILS cudails;
	//bestClustering = cudails.executeILS(&NoConstruct, vnd, g, iter, 1, perturbationLevelMax, problem, info);
	*/

	// 7. Stops the timer and stores the elapsed time
	timer.stop();
	end_time = timer.elapsed();

	// 8. Write the results into ostream os, using csv format
	// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
	timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
	iterationResults << executionCount << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
			<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
			<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction)
			<< "," << improvementCount << "," << totalFrustratedSolutions << "\n";
	Imbalance bestValue = bestClustering.getImbalance();
	int iterationValue = 1;
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
				<< "," << bestValue.getPositiveValue()
				<< "," << bestValue.getNegativeValue()
				<< setprecision(0)
				<< "," << bestClustering.getNumberOfClusters()
				<< "," << (iterationValue+1)
				<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) // timeSpentOnBestSolution
				<< "," << improvementCount
				<< "," << totalFrustratedSolutions << "\n";

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after 1 ILS parallel graph iteration: I(P) = " << bestClustering.getImbalance().getValue();

	// Saves the iteration results to csv file
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Save the time spent at each iteration by each process to csv file
	generateOutputFile(problem, iterationTimeSpent, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timespent-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);

	return bestClustering;
}

// TODO: REIMPLEMENT THIS WHOLE METHOD IN ACCORDANCE TO BOOST PARALLEL BGL
ProcessClustering ImbalanceSubgraphParallelILS::preProcessSplitgraphPartitioning(clusteringgraph::SignedGraph *g,
		problem::ClusteringProblem& problem, bool partitionByVertex) {
	// 1. Divide the graph into (numberOfSlaves + 1) non-overlapping parts

	unsigned int numberOfProcesses = numberOfSlaves + 1;
	string strNP = boost::lexical_cast<std::string>(numberOfProcesses);
	BOOST_LOG_TRIVIAL(info) << "Invoking partitioning for spit graph (numberOfProcesses = " << strNP << ")...";

	/**  ETAPA DE PRE-PROCESSAMENTO - GERACAO DA PARTICAO INICIAL DO GRAFO:
	 n: numero de vertices
	 p: numero de processos
	P0: Gr = G (para o processo 0);
	Loop sequencial, cada processador executa o loop, um processador pi por vez:
	{
		Si = vazio;
		Escolha um nó com menor cadinalidade de arestas negativas entre do grafo residual Gr (o que vai aumentar menos o imbalance)
		(opção 2: maior cardinalidade de arestas positivas)

		Enquanto a cardinalidade do conjunto Si for menor que o n / p:
		{
			Escolha um nó ni pert (Gr \ Si) com maior cardinalidade de arestas positivas entre ni e Si;
			Si = Si + ni;
		}
		Gr = Gr – Si
		Envia Gr para o processo pi+1
	}
	Rodar o ILS dentro de cada subgrafo Gi;
	Merge das soluçoes;
	Busca 1 iteração do ILS sobre a solucao global;  */

	// *** INITIALIZATION: Gr = G (in the begginning, the residual graph equals the whole graph);
	long n = g->getN();
	long m = g->getM();
	// initializes a vector containing vertex degree info (vertex number and its pos and neg degrees)
	std::vector<VertexDegree> Gr;
	for(long i = 0; i < n; i++) {
		Gr.push_back(VertexDegree(i, 0.0, 0.0));
	}
	long desiredCardinality = 0;
	// a list containing the partition number each vertex belongs to
	ClusterArray splitgraphClusterArray(n);
	if(partitionByVertex) {
		desiredCardinality = long(floor(n / (double)numberOfProcesses));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " vertices.";
	} else {
		desiredCardinality = long(floor(m / (double)numberOfProcesses));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " edges.";
	}

	// *** STEP A => PREPROCESSING PHASE OF DISTRIBUTED METAHEURISTIC: New split graph partitioning here
	// Triggers local processing time calculation
	double timeSpent = 0.0;
	double totalTimeSpent = 0.0;
	boost::timer::cpu_timer ftimer;
	ftimer.start();
	boost::timer::cpu_times fstart_time = ftimer.elapsed();

	// TODO estudar como usar o grafo distribuido por inteiro: http://www.boost.org/doc/libs/1_60_0/libs/graph_parallel/doc/html/vertex_list_adaptor.html
	// For each process pi
	for(int pi = 0; pi < numberOfProcesses; pi++) {
		BOOST_LOG_TRIVIAL(debug) << "Processing partition for processor " << pi;
		// *** STEP 1: Chooses a node with the smallest negative cardinality INSIDE the residual graph Gr
		//  (will probably result in the smallest increase in the imbalance)
		clusteringgraph::Clustering GrCluster(*g);  // uses Clustering data structure to control the elements in the set
		assert(Gr.size() > 0);
		GrCluster.addCluster(*g, problem, Gr[0].id, false);
		for(int i = 1; i < Gr.size(); i++) {
			GrCluster.addNodeToCluster(*g, problem, Gr[i].id, 0, false);
		}
		BOOST_LOG_TRIVIAL(debug) << "Created initial residual graph list of size " << GrCluster.getClusterSize(0);

		// updates the negative edge sum of each vertex ni in Gr (neg edge sum between ni and Si)
		for(int i = 0; i < Gr.size(); i++) {
			long ni = Gr[i].id;
			Gr[i].negativeDegree = fabs(g->getNegativeEdgeSumBetweenVertexAndClustering(ni, GrCluster.getClusterArray()));
		}

		// Obtains the vertex in Gr which has the SMALLEST absolute negative edge sum
		std::vector<VertexDegree>::iterator it_min = std::min_element(Gr.begin(), Gr.end(), neg_degree_ordering_asc());  // O(n)
		long vertex_y = it_min->id;
		BOOST_LOG_TRIVIAL(debug) << "The vertex in Gr with the smallest neg edge sum is " << vertex_y << " with sum = " << it_min->negativeDegree;
		// the Si clustering will have only 1 cluster, cluster zero
		clusteringgraph::Clustering Si(*g);  // Si = empty, uses Clustering data structure to control the elements in the set
		Si.addCluster(*g, problem, vertex_y, false);
		Gr.erase(it_min);

		ClusterArray cTemp = Si.getClusterArray();
		long numberOfEdges = 0;
		for(int i = 0; i < Gr.size(); i++) {
			long ni = Gr[i].id;
			Gr[i].positiveDegree = g->getPositiveEdgeSumBetweenVertexAndClustering(ni, Si.getClusterArray());
		}

		while(Gr.size() > 0) {
			std::vector<VertexDegree>::iterator it_max = std::max_element(Gr.begin(), Gr.end(), pos_degree_ordering_asc());  // O(n)
			BOOST_LOG_TRIVIAL(debug) << "ClusterSize = " << Si.getClusterSize(0) << ": The vertex in (Gr - Si) with the biggest pos edge sum is " << it_max->id << " with sum = " << it_max->positiveDegree;

			// Chooses a node max_ni belonging to (Gr \ Si) with the biggest cardinality of positive edges between ni and Si;
			long max_ni = it_max->id;
			// Si = Si + ni
			Si.addNodeToCluster(*g, problem, max_ni, 0, false);
			// Gr = Gr - ni
			Gr.erase(it_max);
			numberOfEdges += g->getOutDegree(max_ni);

			// incremental update of Gr
			// updates the positive edge sum of each vertex ni in Gr (pos edge sum between ni and Si)
			boost::timer::cpu_timer timer;
			timer.start();
			boost::timer::cpu_times start_time = timer.elapsed();

			boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g->graph));
			ParallelGraph::edge_descriptor e;
			ParallelGraph::out_edge_iterator f2, l2;
			std::vector<double> deltaPosDegree(n, (double)0);
			// keeps an array containing every positive edge between ni_max and every other vertex -> O(n)
			for (boost::tie(f2, l2) = out_edges(vertex(max_ni, *(g->graph)), *(g->graph)); f2 != l2; ++f2) {
				e = *f2;
				long j = target(*f2, *(g->graph)).local;
				if(ew[e].weight > 0) {
					deltaPosDegree[j] = ew[e].weight;
				}
			}
			// updates the positiveDegree of between the vertices in Gr and the cluster Si
			for(int i = 0; i < Gr.size(); i++) {
				long ni = Gr[i].id;
				Gr[i].positiveDegree += deltaPosDegree[ni];
			}
			// Last step. Stops the timer and stores the elapsed time
			timer.stop();
			boost::timer::cpu_times end_time = timer.elapsed();
			timeSpent += (end_time.wall - start_time.wall) / double(1000000000);

			if(partitionByVertex) {
				if(Si.getClusterSize(0) >= desiredCardinality)  break;
			} else {
				if(numberOfEdges >= desiredCardinality)  break;
			}
		}
		BOOST_LOG_TRIVIAL(info) << "Partition ready.";
		BOOST_LOG_TRIVIAL(info) << "Time spent on recalculation: " << timeSpent;
		// Gr = Gr – Si  TODO: subtrair a lista partialGraph da lista residualGraph -> ja feito acima
		// Envia Si para o processo pi+1
		ClusterArray SiArray = Si.getClusterArray();
		for(int i = 0; i < n; i++) {
			if(SiArray[i] >= 0) {
				splitgraphClusterArray[i] = pi;
			}
		}
	}

	// Last step. Stops the timer and stores the elapsed time
	ftimer.stop();
	boost::timer::cpu_times fend_time = ftimer.elapsed();
	totalTimeSpent = (fend_time.wall - fstart_time.wall) / double(1000000000);
	BOOST_LOG_TRIVIAL(info) << "Total time spent on construction: " << totalTimeSpent << " s.";
	BOOST_LOG_TRIVIAL(info) << "Time spent on recalculation of degrees: " << timeSpent << " s ( " << (100 * timeSpent / totalTimeSpent) << "\% ).";

	// processes the initial solution generated by the split graph partitioning
	BOOST_LOG_TRIVIAL(info) << "Generated a cluster array of " << splitgraphClusterArray.size() << " vertices.";

	ImbalanceMatrix processClusterImbMatrix = util.calculateProcessToProcessImbalanceMatrix(*g, splitgraphClusterArray,
				this->vertexImbalance, numberOfSlaves + 1);
	ProcessClustering Cc(g, problem, splitgraphClusterArray, processClusterImbMatrix);

	return Cc;
}

clusteringgraph::Clustering ImbalanceSubgraphParallelILS::runDistributedILS(ConstructClustering *construct,
		VariableNeighborhoodDescent *vnd, clusteringgraph::SignedGraph *g, const int& iter, const int& iterMaxILS,
		const int& perturbationLevelMax, problem::ClusteringProblem& problem, ExecutionInfo& info,
		ProcessClustering& splitgraphClustering, clusteringgraph::Clustering& currentSolution) {

	clusteringgraph::Clustering c = distributeSubgraphsBetweenProcessesAndRunILS(construct, vnd, g,
			iter, iterMaxILS, perturbationLevelMax, problem, info, splitgraphClustering);
	int retries = 1;

	// TODO return the last global solution found by distributed ILS or the less worse one?
	while((c.getImbalance().getValue() >= currentSolution.getImbalance().getValue()) and (retries <= MAX_ILS_RETRIES)) {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Found worse global solution with I(P) = "
				<< c.getImbalance().getValue() << " (best known I(P) = " << currentSolution.getImbalance().getValue()
				<< "). Retry " << retries << "...";
		retries++;
		c = distributeSubgraphsBetweenProcessesAndRunILS(construct, vnd, g,
				iter, iterMaxILS, perturbationLevelMax, problem, info, splitgraphClustering);
	}
	return c;
}

clusteringgraph::Clustering ImbalanceSubgraphParallelILS::distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
		VariableNeighborhoodDescent *vnd, clusteringgraph::SignedGraph *g, const int& iter, const int& iterMaxILS,
		const int& perturbationLevelMax, problem::ClusteringProblem& problem, ExecutionInfo& info,
		ProcessClustering& splitgraphClustering) {

	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	// long n = g->getN();
	unsigned int numberOfProcesses = numberOfSlaves + 1;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == problem::ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}
	// obtain the current distribution of vertices between processes
	// OLD CODE: std::vector< std::vector< long > > verticesInProcess = splitgraphClustering.getListOfVerticesInEachProcess(g);
	// typedef graph_traits<ParallelGraph>::vertex_descriptor Key;
	// VertexProcessorMap to_processor_map = get(vertex_rank, *g->graph);
	// TODO obter a propriedade certa da PBGL que retorna o processador a que o vertice pertence
	BOOST_LOG_TRIVIAL(info) << "distributeSubgraphsBetweenProcessesAndRunILS started...";
	BOOST_LOG_TRIVIAL(info) << "Distributing the graph between " << (numberOfSlaves + 1) << " processes.";

	// There are numberOfProcesses subgraphs
	/*
	std::vector<SubGraph> subgraphList;
	// each subgraph will have a subset of the main graph's nodes and edges, based on the previous clustering
	for(int p = 0; p < numberOfProcesses; p++) {
		SubGraph sg = (*(g->graph)).create_subgraph(); //verticesInProcess[k].begin(), verticesInProcess[k].end());
		subgraphList.push_back(sg);  // --> SUBGRAPH COPY CTOR NOT WORKING!!!
		for(std::vector<long>::iterator it = verticesInProcess[p].begin(); it != verticesInProcess[p].end(); it++) {
			add_vertex(*it, subgraphList.back());
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << p;
		}

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_edges = " <<
				num_edges(subgraphList.back()) << " , num_vertices = " << num_vertices(subgraphList.back());
	} */
	std::vector<InputMessageParallelILS> messagesSent;
	for(int i = 0; i < numberOfSlaves + 1; i++) {
		// 2. Distribute numberOfSlaves graph parts between the ILS Slave processes
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
							problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
							numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k, true, NULL, true,
							/*runILS=*/true, /*redistributeVertices=*/false);
		if(k < 0) {
			imsg.setClustering(&CCclustering);
		}
		//imsg.setVertexList(verticesInProcess[i+1]); NAO EH MAIS NECESSARIO
		// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
		// TODO obter verticesInProcess !
		imsg.cudaEnabled = is_overloaded_process(verticesInProcess[i]);
		messagesSent.push_back(imsg);
		if(i > 0) {  // only sends messages to worker processes, the leader process will process his own message after this loop
			world.send(slaveList[i - 1], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
			BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Message sent to process " << slaveList[i - 1];
			BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Size of ILS Input Message: " << (sizeof(imsg)/1024.0) << "kB.";
		}
	}
	// 2.1. the leader does its part of the work: runs ILS using the first part of the divided graph
	VariableNeighborhoodDescent localVND = *vnd;
	localVND.setTimeLimit(LOCAL_ILS_TIME_LIMIT);
	OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(messagesSent[0], g);

	// 3. the leader receives the processing results
	OutputMessage omsg;
	std::map<int, OutputMessage> messageMap;
	leaderProcessingMessage.globalVertexId = this->verticesInLeaderProcess;
	messageMap[0] = leaderProcessingMessage;
	// N: the number of vertices in the global graph
	long N = leaderProcessingMessage.num_vertices;

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
	for(int i = 0; i < numberOfSlaves; i++) {
		OutputMessage omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ". Obj = " <<
				omsg.clustering.getImbalance().getValue();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// stores the time spent by each process
		// timeSpent[i+1] = omsg.timeSpent;
		messageMap[procNum] = omsg;
		N += omsg.num_vertices;
	}
	assert(N == g->getGlobalN());

	// Global cluster array
	ClusterArray globalClusterArray(N, 0);
	ClusterArray splitClusterArray = splitgraphClustering.getClusterArray();
	// g->setGlobalN(N);
	BOOST_LOG_TRIVIAL(debug) << "Global cluster size is " << N;
	long clusterOffset = 0;
	double internalImbalancePosSum = 0.0;
	double internalImbalanceNegSum = 0.0;
	std::vector<unsigned int> clusterProcessOrigin;  // which process a cluster belongs to?
	std::vector<Imbalance> internalProcessImbalance(numberOfProcesses, Imbalance());  // the internal imbalance calculated by each process

	 // Create a global index map, which takes vertex descriptors and
	// turns them into numbers in the range [0, N), where N is the
	// total number of vertices in the distributed graph.
	/* typedef boost::parallel::global_index_map<vtkGraphDistributedVertexIndexMap, vtkVertexGlobalMap> GlobalIndexMap; */

	// https://github.com/Kitware/VTK/blob/c3ec2495b183e3327820e927af7f8f90d34c3474/Infovis/Parallel/vtkPBGLGraphAdapter.h
	// http://www.osl.iu.edu/research/pbgl/documentation/vertex_list_adaptor.html#id2
	/* GlobalIndexMap globalIndexMap(process_group(*(g->graph)),
								  g->getN(),
								  MakeDistributedVertexIndexMap(*(g->graph)),
								  get(boost::vertex_global, *(g->graph))); */

	/* vertex_list_adaptor<ParallelGraph> vla(*(g->graph), globalIndexMap); */

	/*
	BGL_FORALL_VERTICES(v, *(g->graph), ParallelGraph)
	{
		int p_owner = owner(v); // process_id_type
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Owner of vertex " << v.local << ": " << p_owner;
	} */

	for(int procNum = 0; procNum < numberOfSlaves + 1; procNum++) {
		OutputMessage omsg = messageMap[procNum];
		// stores the time spent by each process
		//timeSpent[i+1] = omsg.timeSpent;
		// 4. Merge the partial solutions into a global solution for the whole graph
		ClusterArray localClusterArray = omsg.clustering.getClusterArray();
		Imbalance internalImbalance = omsg.clustering.getImbalance();
		internalProcessImbalance[procNum] = internalImbalance;
		internalImbalancePosSum += internalImbalance.getPositiveValue();
		internalImbalanceNegSum += internalImbalance.getNegativeValue();

		for(long v = 0; v < localClusterArray.size(); v++) {
			// Obtains vertex v's number in the global graph
			long vglobal = omsg.globalVertexId[v];
			splitClusterArray[vglobal] = procNum;
		}

		long msg_nc = omsg.clustering.getNumberOfClusters();
		if(msg_nc > 0) {
			BOOST_LOG_TRIVIAL(info) << "localClusterArray.size() = " << localClusterArray.size();
			BOOST_LOG_TRIVIAL(info) << "omsg.globalVertexId.size() = " << omsg.globalVertexId.size();
			assert(omsg.num_vertices == omsg.globalVertexId.size());
			assert(localClusterArray.size() == omsg.globalVertexId.size());
			for(long v = 0; v < localClusterArray.size(); v++) {
				// Obtains vertex v's number in the global graph
				long vglobal = omsg.globalVertexId[v];
				globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
				/*
				BOOST_LOG_TRIVIAL(info) << "Vertex " << v << " in local solution of process P" << procNum <<
						" becomes vertex " << vglobal <<
						" in global solution and belongs to cluster " << localClusterArray[v] <<
						", that is, global cluster " << globalClusterArray[vglobal];
				*/
			}
			clusterOffset += msg_nc;
			// all the clusters in this interval belong to process 'stat.source()'
			for(int clusterCount = 0; clusterCount < msg_nc; clusterCount++) {
				clusterProcessOrigin.push_back(procNum);
			}
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << procNum << ": num_edges = " <<
						omsg.num_edges << " , num_vertices = " << omsg.num_vertices << ", I(P) = " << omsg.clustering.getImbalance().getValue()
						 << ", k = " << msg_nc;
	}

	// Calculates the external imbalance sum (between processes)
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Calculating external process imbalance matrix...";
	ImbalanceMatrix processClusterImbMatrix = util.calculateProcessToProcessImbalanceMatrix(*g, splitClusterArray,
					this->vertexImbalance, numberOfSlaves + 1);
	ProcessClustering Cc(g, problem, splitClusterArray, processClusterImbMatrix);
	splitgraphClustering.setInterProcessImbalanceMatrix(processClusterImbMatrix);
	Imbalance externalImbalance = util.calculateExternalImbalanceSumBetweenProcesses(processClusterImbMatrix);

	// 5. Builds the clustering with the merge of each process local ILS result
	// EVITA O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Building global clustering object...";
	clusteringgraph::Clustering globalClustering(globalClusterArray, *g, problem, internalImbalancePosSum + externalImbalance.getPositiveValue(),
			internalImbalanceNegSum + externalImbalance.getNegativeValue(), clusterProcessOrigin, internalProcessImbalance);

	clusteringgraph::validation::CCImbalanceCalculator calc = CCImbalanceCalculator::instance(g->getGraphFileLocation());
	Imbalance validation = calc.objectiveFunction(globalClusterArray);
	if(fabs(validation.getValue() - globalClustering.getImbalance().getValue()) > EPS) {
		BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] distributeSubgraphsBetweenProcessesAndRunILS: I(P) DOES NOT MATCH! Correct I(P) = "
						<< std::setprecision(5) << std::fixed << validation.getValue()
						<< " vs I(P) = " << globalClustering.getImbalance().getValue();
	} else {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] distributeSubgraphsBetweenProcessesAndRunILS: I(P) OK!";
	}

	// FIXME remover codigo de validacao
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** VALIDATION *** Calculating internal process imbalance matrix...";
	std::vector<Imbalance> internalProcessImbalanceValidation = util.calculateProcessInternalImbalance(*g,
					splitClusterArray, globalClusterArray, numberOfSlaves + 1);
	Imbalance internalImbalanceValidation = util.calculateInternalImbalanceSumOfAllProcesses(internalProcessImbalanceValidation);
	Imbalance tmpinternalImbalance = util.calculateInternalImbalanceSumOfAllProcesses(internalProcessImbalance);
	bool imbIgual = tmpinternalImbalance == internalImbalanceValidation;
	BOOST_LOG_TRIVIAL(info) << "OpX *** Os imbalances internos sao iguais: " << imbIgual;
	if(not imbIgual) {
		stringstream ss1, ss2;
		ss1 << "Vetor imbalance algoritmo: ";
		ss2 << "Vetor imbalance VALIDACAO: ";
		for(int p = 0; p < numberOfSlaves + 1; p++) {
			ss1 << internalProcessImbalance[p].getValue() << " ";
			ss2 << internalProcessImbalanceValidation[p].getValue() << " ";
		}
		BOOST_LOG_TRIVIAL(info) << ss1.str();
		BOOST_LOG_TRIVIAL(info) << ss2.str();
	}
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** VALIDATION *** Calculating external process imbalance matrix...";
	ImbalanceMatrix processClusterImbMatrixValidation = util.calculateProcessToProcessImbalanceMatrix(*g, splitClusterArray,
						this->vertexImbalance, numberOfSlaves + 1);
	bool igual = ublas::equals(processClusterImbMatrix.pos, processClusterImbMatrixValidation.pos, 0.001, 0.1);
	bool igual2 = ublas::equals(processClusterImbMatrix.neg, processClusterImbMatrixValidation.neg, 0.001, 0.1);
	BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2;
	processClusterImbMatrix = processClusterImbMatrixValidation;


	// 6. Stores the time spent in this iteration of distributed MH
	// timeSpentAtIteration.push_back(timeSpent);
	util.validaSplitgraphArray(*g, splitgraphClustering, globalClustering);
	clusteringgraph::Clustering cl(splitClusterArray, *g, problem, 0.0, 0.0);
	splitgraphClustering.setSplitgraphClustering(cl);

	return globalClustering;
}

bool ImbalanceSubgraphParallelILS::moveCluster1opt(clusteringgraph::SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
		clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses,
		ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   Move  cluster  1 - opt  ***********************************
	// b) The master process chooses a cluster to be moved to another process
	// select a cluster c to move - c is the cluster that causes the biggest imbalance
	// and, at the same time, giving priority to moving a cluster from an overloaded process to a less-loaded one
	// (priotitize vertex load balancing between processes)
	// TODO the biggest imbalance in the whole graph or only between the pair of processes x and y?
	// for now it is the vertex c that belongs to the processPair and causes the biggest imbalance in the whole graph
	// Tries to move the top 25% of most imbalanced clusters
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	const ImbalanceMatrix initialProcessClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best global solution so far: I(P) = " << bestClustering.getImbalance().getValue();
	assert(splitgraphClusterArray.size() == g->getGlobalN());
	assert(bestClustering.getClusterArray().size() == g->getGlobalN());
	long nc = bestClustering.getNumberOfClusters();
	std::vector<unsigned int> clusterProcessOrigin = bestClustering.getClusterProcessOrigin();
	// ----------------------------------------------------------------------
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	std::vector<Coordinate> clusterImbalanceList = util.obtainListOfImbalancedClusters(*g, bestClustering);
	// sorts the list of imbalanced clusters according to the percentage of imbalance (value)
	std::sort(clusterImbalanceList.begin(), clusterImbalanceList.end(), coordinate_ordering_desc());
	long numberOfImbalancedClustersToBeMoved = (int)ceil(clusterImbalanceList.size() * PERCENTAGE_OF_MOST_IMBALANCED_CLUSTERS_TO_BE_MOVED);
	if(numberOfImbalancedClustersToBeMoved == 0) {
		numberOfImbalancedClustersToBeMoved = 1;
	}
	bool foundBetterSolution = false;
	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	long n = g->getGlobalN();
	// unsigned int numberOfProcesses = numberOfSlaves + 1;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == problem::ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}
	// obtains the number of clusters from each process 'procNum'
	int np = bestSplitgraphClustering.getNumberOfClusters();
	assert(np == world.size());
	std::vector<long> numberOfClustersInProcess;
	for(int px = 0; px < np; px++) {
		numberOfClustersInProcess.push_back(util.obtainListOfClustersFromProcess(*g, bestClustering, px).size());
	}
	const clusteringgraph::Clustering initialGlobalClustering = bestClustering;
	ProcessClustering initialSplitgraphClustering = bestSplitgraphClustering;

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** moveCluster1opt";
	// stores the list of clusters where x is the cluster that will be moved to another process y

	// *** STEP 2: POPULATE LIST OF MOVEMENTS
	std::vector<Coordinate> movementList;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of potential clusters to be moved: " << numberOfImbalancedClustersToBeMoved;
	long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0);
	for(long nic = 0; nic < numberOfImbalancedClustersToBeMoved; nic++) {
	// for(std::vector<long>::iterator clIter = movementList.begin(); clIter != movementList.end(); clIter++) {
	// 	    int clusterNum = *clIter;
		long clusterToMove = clusterImbalanceList[nic].x;
		// finds out to which process the cluster belongs to
		int sourceProcess = clusterProcessOrigin[clusterToMove];
		// avoid moving the cluster in case it is the only one inside the process
		BOOST_LOG_TRIVIAL(trace) << "[Parallel ILS SplitGraph] Testing potential cluster " << clusterToMove << " from process " << sourceProcess
				<< " which has " << numberOfClustersInProcess[sourceProcess] << " clusters in total.";
		if(numberOfClustersInProcess[sourceProcess] <= 1)  continue;
		std::vector<long> listOfMovedVertices = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMove);
		// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Global cluster " << clusterToMove << " of size " <<
		//		bestClustering.getClusterSize(clusterToMove) << " and imbalance = " << clusterImbalanceList[nic].value << ".";
		// TODO VERIFICAR SE COM O NOVO PROCESS GROUP DA BOOST PARALLEL BGL O NUMERO DO VERTICE TEM QUE COMECAR COM 1!
		for(int destinationProcess = 0; destinationProcess < np; destinationProcess++) {
			// TODO evaluate the effect of restricting the number of vertices in each process (avoid process overload)
			if( (destinationProcess != sourceProcess) and ( bestSplitgraphClustering.getClusterSize(destinationProcess) +
								listOfMovedVertices.size() < maxVerticesAllowedInProcess ) ) {
				Coordinate movement(clusterToMove, destinationProcess);
				movementList.push_back(movement);
			}
		}
	}
	// -------------------------------------------

	// *** STEP 3: EVALUATE SEVERAL MOVEMENTS IN PARALLEL (DISTRIBUTED ALGORITHM)
	ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	while(not movementList.empty()) {
		Coordinate movement = movementList.back();
		movementList.pop_back();
		long clusterToMove = movement.x;
		int destinationProcess = movement.y;

		unsigned int sourceProcess = clusterProcessOrigin[clusterToMove];
		std::vector<long> listOfModifiedVertices = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMove);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] numClusters of sourceProcess is " << numberOfClustersInProcess[sourceProcess];

		// STEP 1: sends the movements to participating processes through MPI
		std::vector< long > verticesInSourceProcess;
		// Realiza a movimentacao dos vertices de um cluster especifico (cluster move)
		ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
		for(long elem = 0; elem < listOfModifiedVertices.size(); elem++) {
			tempSplitgraphClusterArray[listOfModifiedVertices[elem]] = destinationProcess;
		}
		std::vector< long > verticesInDestinationProcess;
		verticesInSourceProcess.clear();
		std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
		for(long i = 0; i < n; i++) {
			long k = tempSplitgraphClusterArray[i];
			assert(k < np);
			verticesInProcess[k].push_back(i);
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] vertex " << i << " goes to process " << k;
			if(k == destinationProcess) {
				verticesInDestinationProcess.push_back(i);
			} else if(k == sourceProcess) {
				verticesInSourceProcess.push_back(i);
			}
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to move global cluster " << clusterToMove
				<< " from process " << sourceProcess << " to process " << destinationProcess;
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] " << listOfModifiedVertices.size() << " vertices will be moved.";

		// Each process will execute 1 ILS procedure for each subgraph
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
							problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
							numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k, true, NULL, true,
							/*runILS=*/true, /*redistributeVertices=*/true);
		if(k < 0) {  imsg.setClustering(&CCclustering);  }
		imsg.splitgraphClusterArray = tempSplitgraphClusterArray;
		InputMessageParallelILS imsg2 = imsg;
		// sends the modified subgraphs that will be solved by ILS
		imsg.setVertexList(verticesInDestinationProcess);
		imsg2.setVertexList(verticesInSourceProcess);
		// FIXME remover o codigo de debug abaixo
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Vertices to be in source process (P " << sourceProcess << "): ";
		Print(verticesInSourceProcess);
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Vertices to be in destination process (P " << destinationProcess << "): ";
		Print(verticesInDestinationProcess);

		// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
		imsg.cudaEnabled = is_overloaded_process(verticesInDestinationProcess);
		imsg2.cudaEnabled = true;
		std::vector<int> participatingProcessList;
		std::vector<InputMessageParallelILS> messagesSent;
		int leaderParticipates = 0, leaderProcess = 0;
		if(destinationProcess != leaderProcess) {
			world.send(destinationProcess, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		} else { leaderParticipates = 1; }
		if(sourceProcess != leaderProcess) {
			world.send(sourceProcess, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg2);
		} else { leaderParticipates = 2; }
		messagesSent.push_back(imsg);
		messagesSent.push_back(imsg2);

		for(int pr = 1; pr < np; pr++) {
			if((pr != destinationProcess) and (pr != sourceProcess)) {
				InputMessageParallelILS imsg4 = imsg;
				imsg4.runILS = true;  // FIXME trocar de volta para false
				imsg4.setVertexList(verticesInProcess[pr]);
				world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg4);
				messagesSent.push_back(imsg4);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << pr << " (runILS = false)";
			}
		}
		participatingProcessList.push_back(destinationProcess);
		participatingProcessList.push_back(sourceProcess);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to processes " << sourceProcess << " (source) and " <<
				destinationProcess << "(destination).";
		// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Size of ILS Input Message: " << (sizeof(imsg)/1024.0) << "kB.";

		// The leader may participate in the source-process movement
		std::map<int, OutputMessage> messageMap;
		if(leaderParticipates) {
			// 2.1. the leader does its part of the work: runs ILS to solve the source subgraph
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking local ILS on leader process (P0).";
			VariableNeighborhoodDescent localVND = *vnd;
			localVND.setTimeLimit(LOCAL_ILS_TIME_LIMIT);
			InputMessageParallelILS inputMsg = messagesSent[leaderParticipates - 1];
			// moves the affected vertices between processes
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(inputMsg, g);
			  // construct, &localVND, g, iter, iterMaxILS, perturbationLevelMax,
			  //		problem, info, verticesInSourceProcess);
			numberOfTestedCombinations += leaderProcessingMessage.numberOfTestedCombinations;
			messageMap[leaderProcess] = leaderProcessingMessage;
		} else {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking redistribution on leader process (P0).";
			InputMessageParallelILS imsg4 = imsg;
			imsg4.runILS = false;
			imsg4.setVertexList(verticesInProcess[0]);
			imsg4.splitgraphClusterArray = tempSplitgraphClusterArray;
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg4, g);
			messageMap[leaderProcess] = leaderProcessingMessage;
		}

		// STEP 3: the leader receives the processing results
		// Maps the messages according to the movement executed
		int totalMessagesSent = messagesSent.size();
		if(leaderParticipates) {
			totalMessagesSent--;
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for " << totalMessagesSent << " return messages...";
		for(int i = 0; i < totalMessagesSent; i++) {
			OutputMessage omsg;
			mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
			int tag = stat.tag();
			int procNum = stat.source();
			if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
				BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
				i++;
				while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
					stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
					i++;
				}
				throw std::invalid_argument( "Error message received from slave. See error logs." );
			}
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			// process the result of the execution of process i
			// sums the total number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// stores the time spent by each process
			// timeSpent[i+1] = omsg.timeSpent;
			messageMap[procNum] = omsg;
		}

		// STEP 4: For each movement, merge the partial solutions into a global solution for the whole graph
		// Merges the partial solutions into a global solution for the whole graph
		// LEADER: includes the leader's processing result as well
		// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
		// structure containing the vertex id in the global graph / clustering
		// process the move to each destination process
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Analyzing movement to process " << destinationProcess;

		// the list containing the internal imbalance of each process local clustering
		std::vector<Imbalance> newInternalProcessImbalance = initialGlobalClustering.getInternalProcessImbalance();
		tempSplitgraphClusterArray = currentSplitgraphClusterArray;
		// removes all the clusters which belong to the processes participating in this cluster movement
		// and calculates the new number of clusters, to be used in the offset below
		clusteringgraph::Clustering tempClustering = initialGlobalClustering;
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Removing all clusters from participating processes...";
		tempClustering.removeAllClustersFromProcess(g, sourceProcess);
		tempClustering.removeAllClustersFromProcess(g, destinationProcess);
		long clusterOffset = tempClustering.getNumberOfClusters();
		std::vector<unsigned int> clusterProcessOrigin = tempClustering.getClusterProcessOrigin();
		ClusterArray globalClusterArray = tempClustering.getClusterArray();

		for(int workerProcess = 0; workerProcess < 2; workerProcess++) {
			BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Worker process " << workerProcess;
			int process = participatingProcessList[workerProcess];
			BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Processing clusters from participating process " << process << "...";
			OutputMessage msg = messageMap[process];
			ClusterArray localClusterArray = msg.clustering.getClusterArray();

			long msg_nc = msg.clustering.getNumberOfClusters();
			if(msg_nc > 0) {
				assert(localClusterArray.size() == msg.globalVertexId.size());
				for(long v = 0; v < localClusterArray.size(); v++) {
					// Obtains vertex v's number in the global graph
					long vglobal = msg.globalVertexId[v];
					// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Processing global vertex " << vglobal;
					// BOOST_LOG_TRIVIAL(debug) << "vglobal = " << vglobal;
					assert(vglobal < globalClusterArray.size());
					globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
				}
				clusterOffset += msg_nc;
				// all clusters in this interval belong to the destination process in the 1-move-cluster
				for(int clusterCount = 0; clusterCount < msg_nc; clusterCount++) {
					clusterProcessOrigin.push_back(process);
				}
			}
			// Atualizar newInternalProcessImbalance com os novos imbalances de cada processo
			Imbalance internalImbalance = msg.clustering.getImbalance();
			newInternalProcessImbalance[process] = internalImbalance;
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph P" << (process) << ": num_edges = " << msg.num_edges
					<< " , num_vertices = " << msg.num_vertices	<< ", I(P) = " << msg.clustering.getImbalance().getValue() << ", k = " << msg_nc;
		}

		// Realiza a movimentacao dos vertices de um cluster especifico (cluster move) para o processo de destino
		for(long elem = 0; elem < listOfModifiedVertices.size(); elem++) {
			// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Global vertex " << listOfModifiedVertices[elem] << " goes to process " << destinationProcess;
			tempSplitgraphClusterArray[listOfModifiedVertices[elem]] = destinationProcess;
		}
		// O RECALCULO DA MATRIZ ABAIXO EH FEITO DE FORMA INCREMENTAL, reutilizando a matrix processClusterImbMatrix
		ImbalanceMatrix tempProcessClusterImbMatrix = initialProcessClusterImbMatrix;
		// Realiza o calculo do delta da matriz de imbalance entre processos
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Updating ProcessToProcessImbalanceMatrix...";
		util.updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
				listOfModifiedVertices, tempProcessClusterImbMatrix, numberOfSlaves + 1);
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";



		// FIXME remover codigo de validacao
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** VALIDATION *** Calculating internal process imbalance matrix...";
		std::vector<Imbalance> internalProcessImbalanceValidation = util.calculateProcessInternalImbalance(*g,
									tempSplitgraphClusterArray, globalClusterArray, numberOfSlaves + 1);
		Imbalance internalImbalanceValidation = util.calculateInternalImbalanceSumOfAllProcesses(internalProcessImbalanceValidation);
		Imbalance tmpinternalImbalance = util.calculateInternalImbalanceSumOfAllProcesses(newInternalProcessImbalance);
		bool imbIgual = tmpinternalImbalance == internalImbalanceValidation;
		BOOST_LOG_TRIVIAL(info) << "OpX *** Os imbalances internos sao iguais: " << imbIgual;
		if(not imbIgual) {
			stringstream ss1, ss2;
			ss1 << "Vetor imbalance algoritmo: ";
			ss2 << "Vetor imbalance VALIDACAO: ";
			for(int p = 0; p < np; p++) {
				ss1 << newInternalProcessImbalance[p].getValue() << " ";
				ss2 << internalProcessImbalanceValidation[p].getValue() << " ";
			}
			BOOST_LOG_TRIVIAL(info) << ss1.str();
			BOOST_LOG_TRIVIAL(info) << ss2.str();
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** SPLITGRAPH ARRAY VALIDATION ***";
		for(int p = 0; p < np; p++) {
			OutputMessage omsg = messageMap[p];
			ClusterArray splitgraphArray = omsg.splitgraphClusterArray;
			BOOST_LOG_TRIVIAL(info) << "omsg.splitgraphClusterArray size is :" << omsg.splitgraphClusterArray.size();
			int count = 0;
			for(int i = 0; i < g->getGlobalN(); i++) {
				if(splitgraphArray[i] >= 0) {
					count++;
					if(splitgraphArray[i] != tempSplitgraphClusterArray[i]) {
						BOOST_LOG_TRIVIAL(error) << "Splitgraph array invalid for i = " << i << " (tempSplitgraphClusterArray[i] = "
								<< tempSplitgraphClusterArray[i] << ") / (omsg.splitgraphClusterArray[i]" << splitgraphArray[i] << ").";
					}
				}
			}
			BOOST_LOG_TRIVIAL(info) << "Total vertex count for process " << p << " is " << count;
		}



		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** VALIDATION *** Calculating external process imbalance matrix...";
		ImbalanceMatrix processClusterImbMatrixValidation = util.calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray,
							this->vertexImbalance, numberOfSlaves + 1);
		bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, processClusterImbMatrixValidation.pos, 0.001, 0.1);
		bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, processClusterImbMatrixValidation.neg, 0.001, 0.1);
		BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2;
		tempProcessClusterImbMatrix = processClusterImbMatrixValidation;



		newInternalProcessImbalance = util.calculateProcessInternalImbalance(*g, currentSplitgraphClusterArray,
								globalClusterArray, numberOfProcesses);
		// TODO VALIDAR SE O NUMERO DO PRIMEIRO PROCESSO COM PARALLEL GRAPH COMECA EM 1 MESMO
		for(int px = 0; px < numberOfProcesses; px++) {
			bestClustering.setProcessImbalance(px, newInternalProcessImbalance[px]);
		}

		// Calculates the internal imbalance sum (inside each process)
		Imbalance internalImbalance = util.calculateInternalImbalanceSumOfAllProcesses(newInternalProcessImbalance);
		// Calculates the external imbalance sum (between processes)
		Imbalance externalImbalance = util.calculateExternalImbalanceSumBetweenProcesses(tempProcessClusterImbMatrix);

		// 5. Builds the clustering with the merge of each process local ILS result
		// EVITA O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
		clusteringgraph::Clustering globalClustering = clusteringgraph::Clustering(globalClusterArray, *g, problem, internalImbalance.getPositiveValue() + externalImbalance.getPositiveValue(),
				internalImbalance.getNegativeValue() + externalImbalance.getNegativeValue(), clusterProcessOrigin, newInternalProcessImbalance);

		// Validacao do calculo da FO
		// FIXME
		clusteringgraph::validation::CCImbalanceCalculator calc = CCImbalanceCalculator::instance(g->getGraphFileLocation());
		Imbalance validation = calc.objectiveFunction(globalClusterArray);
		if(fabs(validation.getValue() - globalClustering.getImbalance().getValue()) > EPS) {
			BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] I(P) DOES NOT MATCH! Correct I(P) = "
				<< std::setprecision(5) << std::fixed << validation.getValue()
				<< " vs I(P) = " << globalClustering.getImbalance().getValue();
		} else {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] move-cluster-1opt: I(P) OK!";
		}

		if(globalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 1-move-Cluster Improved solution found! I(P) = "
					<< globalClustering.getImbalance().getValue();
			BOOST_LOG_TRIVIAL(info) << "Move from process " << sourceProcess << " to process " << destinationProcess;
			// BOOST_LOG_TRIVIAL(info) << "new sourceProcess imbalance: " << leaderClustering.getImbalance().getValue();
			// BOOST_LOG_TRIVIAL(info) << "new destinationProcess imbalance: " << messageMap[workerProcess].clustering.getImbalance().getValue();

			// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
			bestClustering = globalClustering;

			bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
			util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestClustering);
			foundBetterSolution = true;
			/*  FIRST IMPROVEMENT ON PARALLEL DISTRIBUTED ILS EVALUATIONS IS DISABLED!
			break;
			*/
		} else if(initialGlobalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
			long sizeOfSourceCluster = bestSplitgraphClustering.getClusterSize(sourceProcess);
			long sizeOfDestCluster = bestSplitgraphClustering.getClusterSize(destinationProcess);
			long newSizeOfSourceCluster = sizeOfSourceCluster - listOfModifiedVertices.size();
			long newSizeOfDestCluster = sizeOfDestCluster + listOfModifiedVertices.size();

			if (globalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost move.";
			} else {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Worse solution found. Fixing to zero-cost move.";
			}
			numberOfFrustratedSolutions++;

			// WARNING: It is possible that the new solution has worse imbalance than it should, since
			//    the local ILS of one of the processes may not find the most efficient local solution.
			//    In these cases, discards the clustering solutions of each process and simply creates
			//    a zero-cost move solution (i.e. same imbalance than previous) based on the moved vertices
			//    and the previous best-known clustering.

			// => simply moves the cluster to the new process, with the same global imbalance value
			// evaluate if this zero-cost move improves the vertex balancing between processes
			if (labs(newSizeOfSourceCluster - newSizeOfDestCluster) <
					labs(sizeOfSourceCluster - sizeOfDestCluster)) {  // this zero-cost move is good

				// retrieves the initial global clustering as base for zero-cost
				bestClustering = initialGlobalClustering;

				// The imbalance matrix between processes has already been updated in the code above
				// Simply moves the cluster to the destination process -> does NOT change the imbalance
				// The internal and external imbalance values of the participating processes must be updated
				moveClusterToDestinationProcessZeroCost(g, bestClustering, bestSplitgraphClustering, clusterToMove,
						sourceProcess, destinationProcess);
				ClusterArray bestClusterArray = bestClustering.getClusterArray();
				std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, bestClusterArray,
																			globalClusterArray, numberOfProcesses);
				Imbalance imbSrc = processInternalImbalance[sourceProcess];
				bestClustering.setProcessImbalance(sourceProcess, imbSrc);
				Imbalance imbDest = processInternalImbalance[destinationProcess];
				bestClustering.setProcessImbalance(destinationProcess, imbDest);

				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost moveCluster1opt improving solution found! I(P) = "
												<< bestClustering.getImbalance().getValue();
				bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
				foundBetterSolution = true;
				util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestClustering);
			} else {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost moveCluster1opt move.";
				// TODO INCLUIR AQUI CHAMADA AO ROLLBACK DO MOVIMENTO DE VERTICES
				std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
				long N = g->getGlobalN();
				for(long i = 0; i < N; i++) {
					long k = currentSplitgraphClusterArray[i];
					verticesInProcess[k].push_back(i);
				}
				// moves the affected vertices between processes
				InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
											problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
											numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, 0, true, NULL, true,
											/*runILS=*/false, /*redistributeVertices=*/true);
				imsg.splitgraphClusterArray = currentSplitgraphClusterArray;
				for(int pr = 1; pr < np; pr++) {
					imsg.setVertexList(verticesInProcess[pr]);
					world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback message sent to process " << pr << " (runILS = false)";
				}
				imsg.setVertexList(verticesInProcess[0]);
				OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg, g);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback of master process done (runILS = false).";
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
				for(int pr = 1; pr < np; pr++) {
					OutputMessage omsg;
					mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
					int tag = stat.tag();
					int procNum = stat.source(), i = 0;
					if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
						BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
						i++;
						while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
							stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
							i++;
						}
						throw std::invalid_argument( "Error message received from slave. See error logs." );
					}
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback-done message received from process " << procNum << ".";
				}
			}
		}
		if(foundBetterSolution) {
			break;
		}
	}
	return foundBetterSolution;
}

bool ImbalanceSubgraphParallelILS::swapCluster1opt(clusteringgraph::SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
		clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses,
		ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   Swap  cluster  1 - opt  ***********************************
	// In this neighborhood the priority is swapping a bigger cluster in an overloaded process
	// 	  with a smaller cluster from a less-loaded process.
	// (prioritize vertex load balancing between processes)
	// Tries to move the top 25% of biggest clusters from all overloaded processes
	//   Vertex division between processes (has nothing to do with a solution to the CC problem)
	// *** STEP 1: INITIALIZATION AND INITIAL FULL DISTRIBUTED ILS
	std::vector<Coordinate> overloadedProcessList = util.obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** swapCluster1opt";
	BOOST_LOG_TRIVIAL(info) << "[SplitGraph swapCluster1opt] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	bool foundBetterSolution = false;
	long n = g->getGlobalN();

	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	// obtains the number of clusters from each process 'procNum'
	int np = bestSplitgraphClustering.getNumberOfClusters();

	// runs a full execution of distributed ILS, in order to find the internal imbalance value of each process
	// TODO estudar remocao da execucao do ILS distribuido inicial neste momento (requer que o zerocost move recalcule os imbalances internos dos processos)
	clusteringgraph::Clustering bestGlobalClustering = bestClustering;
	// Clustering bestGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info, splitgraphClusterArray, processClusterImbMatrix, bestClustering);
	const clusteringgraph::Clustering initialGlobalClustering = bestGlobalClustering;
	ProcessClustering initialSplitgraphClustering = bestSplitgraphClustering;
	const ClusterArray initialSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestGlobalClustering.getClusterArray();
	const ImbalanceMatrix initialProcessClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best global solution so far: I(P) = " << bestGlobalClustering.getImbalance().getValue();
	long nc = bestGlobalClustering.getNumberOfClusters();
	std::vector<unsigned int> clusterProcessOrigin = bestGlobalClustering.getClusterProcessOrigin();
	// ----------------------------------------------------------------------

	// *** STEP 2: POPULATE LIST OF MOVEMENTS
	// populates the list of movements (x, y) where x is the cluster that will be swapped to process y
	std::vector<Coordinate> movementListA, movementListB;
	// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] List of clusters to be swapped: ";
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = util.obtainListOfClustersFromProcess(*g, initialGlobalClustering, procSourceNum);
		// if there are no clusters in the process, skips the movement
		if(bigClustersList.size() == 0)  continue;
		// sorts the list of big clusters according to the number of vertices, descending order
		std::sort(bigClustersList.begin(), bigClustersList.end(), coordinate_ordering_desc());
		// TODO possible parametrization here! Parametrize the quantity of big clusters to be moved, one at a time
		long numberOfClustersToBeSwapped = bigClustersList.size();
		/*
		long numberOfClustersToBeSwapped = (int)ceil(bigClustersList.size() * PERCENTAGE_OF_MOST_IMBALANCED_CLUSTERS_TO_BE_MOVED);
		if(numberOfClustersToBeSwapped == 0) {
			numberOfClustersToBeSwapped = 1;
		} */
		for(long clusterCount = 0; clusterCount < numberOfClustersToBeSwapped; clusterCount++) {
			long clusterToSwapA = bigClustersList[clusterCount].x;
			// finds out to which process the cluster belongs to
			std::vector<long> listOfMovedVerticesFromProcessA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapA);
			int currentProcess = initialSplitgraphClusterArray[listOfMovedVerticesFromProcessA[0]];
			if(currentProcess != procSourceNum) {
				BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] Inconsistency detected between source process numbers.";
			}
			assert(currentProcess == procSourceNum);
			//BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Swap from process " << procSourceNum <<
			//		": Moving global cluster " << clusterToSwapA << " of size " << initialGlobalClustering.getClusterSize(clusterToSwapA);

			long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0);
			for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {
				if(procDestNum != procSourceNum) {
					std::vector<Coordinate> clusterListProcessB = util.obtainListOfClustersFromProcess(*g, initialGlobalClustering, procDestNum);
					// gets the cluster from process B (destination) which has less elements
					Coordinate clusterToSwapB = *std::min_element(clusterListProcessB.begin(), clusterListProcessB.end(), coordinate_ordering_asc());
					std::vector<long> listOfMovedVerticesFromProcessB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapB.x);

					// TODO only makes the swap if the max vertex threshold is respected - check if it is worth doing this
					if ( (listOfMovedVerticesFromProcessA.size() > maxVerticesAllowedInProcess)
							or ( (initialSplitgraphClustering.getClusterSize(procDestNum) + listOfMovedVerticesFromProcessA.size()
							- listOfMovedVerticesFromProcessB.size() < maxVerticesAllowedInProcess) 
							     and (initialSplitgraphClustering.getClusterSize(procSourceNum) + listOfMovedVerticesFromProcessB.size()
							- listOfMovedVerticesFromProcessA.size() < maxVerticesAllowedInProcess) ) ) {
						//BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to swap global cluster " << clusterToSwapA
						//		<< " to process " << procDestNum;
						// the cluster 'clusterToSwapA' will be swapped with the smallest cluster from process procDestNum
						Coordinate movementA(clusterToSwapA, procDestNum);
						movementListA.push_back(movementA);
						Coordinate movementB(clusterToSwapB.x, procSourceNum);
						movementListB.push_back(movementB);
					}
				}
			}
		}
	}
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of potential clusters to be swapped: " << movementListA.size();
	// -------------------------------------------

	// *** STEP 3: EVALUATE SEVERAL MOVEMENTS IN PARALLEL (DISTRIBUTED ALGORITHM)
	// ==> Uses parallelism to evaluate more than one movement at the same time, according to the number of available processes
	int numberOfParallelEvaluations = 1; // (int)std::floor(numberOfProcesses / double(2.0));
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of cluster movements evaluated in parallel: " << numberOfParallelEvaluations;
	const ClusterArray currentSplitgraphClusterArray = initialSplitgraphClustering.getClusterArray();
	while(not movementListA.empty()) {
		// generates a list with the cluster-process movements that will be executed in parallel
		std::vector<Coordinate> movementsInParallelA, movementsInParallelB;
		for(int mov = 0; mov < numberOfParallelEvaluations; mov++) {
			movementsInParallelA.push_back(movementListA.back());
			movementListA.pop_back();
			movementsInParallelB.push_back(movementListB.back());
			movementListB.pop_back();
			if(movementListA.empty())  break;
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Processing " << movementsInParallelA.size() << " movements...";

		// Executes the movements inside the list movementsInParallel
		// STEP 3.1: sends the movements to each pair of processes through MPI
		OutputMessage leaderOutputMessage;
		std::vector<InputMessageParallelILS> messagesSent;
		int leaderParticipates = 0;
		std::vector<int> participatingProcessList;
		ClusterArray tempSplitgraphClusterArray;
		std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
		for(int mov = movementsInParallelA.size() - 1; mov >= 0; mov--) {
			long clusterToSwapA = movementsInParallelA[mov].x;
			int destinationProcess = movementsInParallelA[mov].y;
			// finds out to which process the swapped cluster belongs to
			long clusterToSwapB = movementsInParallelB[mov].x;
			int sourceProcess = movementsInParallelB[mov].y;

			int workerProcess1 = sourceProcess;
			int workerProcess2 = destinationProcess;
			participatingProcessList.push_back(workerProcess1);
			participatingProcessList.push_back(workerProcess2);

			std::vector<long> listOfMovedVerticesFromProcessA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapA);
			std::vector<long> listOfMovedVerticesFromProcessB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapB);
			tempSplitgraphClusterArray = initialSplitgraphClusterArray;
			// SWAP Step 1: Move the vertices from a specific cluster from source process (cluster move from A to B)
			for(long elem = 0; elem < listOfMovedVerticesFromProcessA.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromProcessA[elem]] = destinationProcess;
			}
			// SWAP Step 2: Move the vertices from the smallest cluster in / from process B (dest), to process A (source)
			for(long elem = 0; elem < listOfMovedVerticesFromProcessB.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromProcessB[elem]] = sourceProcess;
			}
			// stores the new list of vertices from each process participating in the movement
			std::vector< long > verticesInSourceProcess;
			std::vector< long > verticesInDestinationProcess;
			for(long i = 0; i < n; i++) {
				long k = tempSplitgraphClusterArray[i];
				verticesInProcess[k].push_back(i);
				if(k == sourceProcess) {
					verticesInSourceProcess.push_back(i);
				} else if(k == destinationProcess) {
					verticesInDestinationProcess.push_back(i);
				}
			}
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to swap global clusters " << clusterToSwapA << " and " << clusterToSwapB
					<< " between processes " << sourceProcess << " and " << destinationProcess;

			if(mov >= 1) {  // 2 ILS procedures will be executed by 2 processes through MPI messages
				InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
									problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
									numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, -1, true, NULL, true,
									/*runILS=*/true, /*redistributeVertices=*/true);
				// if(k < 0) {  imsg.setClustering(&CCclustering);  }
				imsg.splitgraphClusterArray = tempSplitgraphClusterArray;
				InputMessageParallelILS imsg2 = imsg;
				// sends the modified subgraphs that will be solved by ILS
				imsg.setVertexList(verticesInSourceProcess);
				imsg2.setVertexList(verticesInDestinationProcess);
				// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
				imsg.cudaEnabled = is_overloaded_process(verticesInSourceProcess);
				imsg2.cudaEnabled = is_overloaded_process(verticesInDestinationProcess);
				int leaderProcess = 0;
				if(workerProcess1 != leaderProcess) {
					world.send(workerProcess1, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
				} else { leaderParticipates = 1; }
				if(workerProcess2 != leaderProcess) {
					world.send(workerProcess2, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg2);
				} else { leaderParticipates = 2; }
				messagesSent.push_back(imsg);
				messagesSent.push_back(imsg2);

				for(int pr = 1; pr < np; pr++) {
					if((pr != destinationProcess) and (pr != sourceProcess)) {
						InputMessageParallelILS imsg4 = imsg;
						imsg4.runILS = false;
						imsg4.setVertexList(verticesInProcess[pr]);
						world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg4);
						messagesSent.push_back(imsg4);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << pr << " (runILS = false)";
					}
				}
				participatingProcessList.push_back(workerProcess1);
				participatingProcessList.push_back(workerProcess2);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to processes " << workerProcess1 << " and " << workerProcess2;
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] verticesInSourceProcess size = " << verticesInSourceProcess.size();
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] verticesInDestinationProcess size = " << verticesInDestinationProcess.size();
				// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Size of ILS Input Message: " << (sizeof(imsg)/1024.0) << "kB.";
			} else {  // 2 ILS procedures will be executed: one by the leader process (this one) and the other through MPI message
				// STEP 3.2: Movements to be executed by workerProcess ranks 0 (leader process) and 1
				InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
									problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
									numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, -1, true, NULL, true,
									/*runILS=*/true, /*redistributeVertices=*/true);
				imsg.splitgraphClusterArray = tempSplitgraphClusterArray;
				// if(k < 0) {  imsg.setClustering(&CCclustering);  }
				imsg.setVertexList(verticesInDestinationProcess);
				// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
				imsg.cudaEnabled = is_overloaded_process(verticesInDestinationProcess);
				world.send(1, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process 1.";
				messagesSent.push_back(imsg);

				// the leader does its part of the work: runs ILS to solve the second subgraph
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking local ILS on leader process (P0).";
				VariableNeighborhoodDescent localVND = *vnd;
				localVND.setTimeLimit(LOCAL_ILS_TIME_LIMIT);
				leaderOutputMessage = runILSLocallyOnSubgraph(imsg, g);
			}
		}

		// STEP 3.3: the leader receives the processing results
		// Maps the messages according to the movement executed
		std::map<int, OutputMessage> messageMap;
		int leaderProcess = 0;
		if(leaderParticipates) {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking local ILS on leader process (P0).";
			VariableNeighborhoodDescent localVND = *vnd;
			localVND.setTimeLimit(LOCAL_ILS_TIME_LIMIT);
			InputMessageParallelILS inputMsg = messagesSent[leaderParticipates - 1];
			// moves the affected vertices between processes
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(inputMsg, g);
			  // construct, &localVND, g, iter, iterMaxILS, perturbationLevelMax,
			  //		problem, info, verticesInSourceProcess);
			numberOfTestedCombinations += leaderProcessingMessage.numberOfTestedCombinations;
			messageMap[leaderProcess] = leaderProcessingMessage;
		} else {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking redistribution on leader process (P0).";
			InputMessageParallelILS imsg4(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
									problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
									numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, -1, true, NULL, true,
									/*runILS=*/true, /*redistributeVertices=*/true);
			imsg4.runILS = false;
			imsg4.setVertexList(verticesInProcess[0]);
			imsg4.splitgraphClusterArray = tempSplitgraphClusterArray;
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg4, g);
		}
		
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
		long totalMessagesSent = messagesSent.size();
		if(leaderParticipates) {
			totalMessagesSent--;
		}
		for(int i = 0; i < totalMessagesSent; i++) {
			OutputMessage omsg;
			mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
			int tag = stat.tag();
			int procNum = stat.source();
			if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
				BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
				i++;
				while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
					stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
					i++;
				}
				throw std::invalid_argument( "Error message received from slave. See error logs." );
			}
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			// process the result of the execution of process i
			// sums the total number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// stores the time spent by each process
			// timeSpent[i+1] = omsg.timeSpent;
			messageMap[procNum] = omsg;
		}

		// STEP 3.4: For each movement, merge the partial solutions into a global solution for the whole graph
		for(int mov = movementsInParallelA.size() - 1; mov >= 0; mov--) {
			long clusterToSwapA = movementsInParallelA[mov].x;
			int destinationProcess = movementsInParallelA[mov].y;
			long clusterToSwapB = movementsInParallelB[mov].x;
			int sourceProcess = movementsInParallelB[mov].y;
			int workerProcess1 = sourceProcess;
			int workerProcess2 = destinationProcess;
			std::vector<long> listOfMovedVerticesFromProcessA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapA);
			std::vector<long> listOfMovedVerticesFromProcessB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToSwapB);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Analyzing movement pair " << mov;
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Swap global clusters " << clusterToSwapA << " and " << clusterToSwapB
								<< " between processes " << sourceProcess << " and " << destinationProcess;

			// the list containing the internal imbalance of each process local clustering
			std::vector<Imbalance> newInternalProcessImbalance = initialGlobalClustering.getInternalProcessImbalance();
			clusteringgraph::Clustering globalClustering;
			// removes all the clusters which belong to the processes participating in this cluster movement
			// and calculates the new number of clusters, to be used in the offset below
			clusteringgraph::Clustering tempClustering = initialGlobalClustering;
			tempClustering.removeAllClustersFromProcess(g, sourceProcess);
			tempClustering.removeAllClustersFromProcess(g, destinationProcess);
			std::vector<unsigned int> newClusterProcessOrigin = tempClustering.getClusterProcessOrigin();
			long clusterOffset = tempClustering.getNumberOfClusters();
			ClusterArray globalClusterArray = tempClustering.getClusterArray();
			ClusterArray previousSplitgraphClusterArray = currentSplitgraphClusterArray;
			ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;

			ImbalanceMatrix tempProcessClusterImbMatrix = initialProcessClusterImbMatrix;
			// SWAP Step 1: Move the vertices from a specific cluster from source process (cluster move from A to B)
			for(long elem = 0; elem < listOfMovedVerticesFromProcessA.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromProcessA[elem]] = destinationProcess;
			}
			// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processClusterImbMatrix
			util.updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
					listOfMovedVerticesFromProcessA, tempProcessClusterImbMatrix, numberOfSlaves + 1);
			previousSplitgraphClusterArray = tempSplitgraphClusterArray;

			// SWAP Step 2: Move the vertices from the smallest cluster in / from process B (dest), to process A (source)
			for(long elem = 0; elem < listOfMovedVerticesFromProcessB.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromProcessB[elem]] = sourceProcess;
			}
			// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processClusterImbMatrix
			util.updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
					listOfMovedVerticesFromProcessB, tempProcessClusterImbMatrix, numberOfSlaves + 1);
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from swap done.";

			// Valida se a matriz incremental e a full sao iguais
			/*
			ImbalanceMatrix tempProcessClusterImbMatrix3 = util.calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray,
					this->vertexImbalance, numberOfSlaves + 1);
			bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix3.pos, 0.001, 0.1);
			bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix3.neg, 0.001, 0.1);
			BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2;
			*/

			// Merges the partial solutions into a global solution for the whole graph
			std::vector<OutputMessage> sourceAndDestinationMovementMsg;
			sourceAndDestinationMovementMsg.push_back(messageMap[workerProcess1]);
			sourceAndDestinationMovementMsg.push_back(messageMap[workerProcess2]);
			// sets the new imbalance for cluster-move source process
			newInternalProcessImbalance[sourceProcess] = messageMap[workerProcess1].clustering.getImbalance();
			// sets the new imbalance for cluster-move destination process
			newInternalProcessImbalance[destinationProcess] = messageMap[workerProcess2].clustering.getImbalance();

			for(int proc = 0; proc <= 1; proc++) {
				OutputMessage msg = sourceAndDestinationMovementMsg[proc];
				ClusterArray localClusterArray = msg.clustering.getClusterArray();

				long msg_nc = msg.clustering.getNumberOfClusters();
				if(msg_nc > 0) {
					assert(localClusterArray.size() == msg.globalVertexId.size());
					for(long v = 0; v < localClusterArray.size(); v++) {
						// Obtains vertex v's number in the global graph
						long vglobal = msg.globalVertexId[v];
						globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
					}
					clusterOffset += msg_nc;
					// all the clusters in this interval belong either to the source or destination process in the swap movement
					for(long c = 0; c < msg_nc; c++) {
						if(proc == 0) {
							newClusterProcessOrigin.push_back(sourceProcess);
						} else {
							newClusterProcessOrigin.push_back(destinationProcess);
						}
					}
				}
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph P" << participatingProcessList[proc] << ": num_edges = " << msg.num_edges
						<< " , num_vertices = " << msg.num_vertices	<< ", I(P) = " << msg.clustering.getImbalance().getValue() << ", k = " << msg_nc;
			}

			// Calculates the internal imbalance sum (inside each process)
			Imbalance internalImbalance = util.calculateInternalImbalanceSumOfAllProcesses(newInternalProcessImbalance);
			// Calculates the external imbalance sum (between processes)
			Imbalance externalImbalance = util.calculateExternalImbalanceSumBetweenProcesses(tempProcessClusterImbMatrix);

			// 5. Builds the clustering with the merge of each process local ILS result
			// EVITA O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
			globalClustering = clusteringgraph::Clustering(globalClusterArray, *g, problem, internalImbalance.getPositiveValue() + externalImbalance.getPositiveValue(),
					internalImbalance.getNegativeValue() + externalImbalance.getNegativeValue(), newClusterProcessOrigin, newInternalProcessImbalance);

			/*
			Clustering validation(globalClusterArray, *g, problem);
			if(fabs(validation.getImbalance().getValue() - globalClustering.getImbalance().getValue()) > EPS) {
				BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] I(P) DOES NOT MATCH! Correct I(P) = "
						<< std::setprecision(5) << std::fixed << validation.getImbalance().getValue()
						<< " vs I(P) = " << globalClustering.getImbalance().getValue();

				stringstream ss;
				ss << "InternalImbalanceVector: ";
				for(int x = 0; x < numberOfProcesses; x++) {
					ss << newInternalProcessImbalance[x].getValue() << " ";
				}
				BOOST_LOG_TRIVIAL(info) << ss.str();

				std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, tempSplitgraphClusterArray,
						globalClusterArray, numberOfProcesses);
			} */

			if(globalClustering.getImbalance().getValue() < bestGlobalClustering.getImbalance().getValue()) {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 1-swap-Cluster Improved solution found! I(P) = "
						<< globalClustering.getImbalance().getValue();

				// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
				bestClustering = globalClustering;
				bestGlobalClustering = globalClustering;

				// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
				bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
				util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestGlobalClustering);
				foundBetterSolution = true;
				/*  FIRST IMPROVEMENT ON PARALLEL DISTRIBUTED ILS EVALUATIONS IS DISABLED!
				break;
				*/
			} else {
				long sizeOfSourceCluster = bestSplitgraphClustering.getClusterSize(sourceProcess);
				long sizeOfDestCluster = bestSplitgraphClustering.getClusterSize(destinationProcess);
				long newSizeOfSourceCluster = sizeOfSourceCluster - listOfMovedVerticesFromProcessA.size() + listOfMovedVerticesFromProcessB.size();
				long newSizeOfDestCluster = sizeOfDestCluster - listOfMovedVerticesFromProcessB.size() + listOfMovedVerticesFromProcessA.size();

				if (globalClustering.getImbalance().getValue() == bestGlobalClustering.getImbalance().getValue()) {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 1-swap Zero-cost move.";
				} else {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 1-swap Worse solution found (Obj = "
							<< globalClustering.getImbalance().getValue() << "). Fixing to zero-cost move.";
				}
				numberOfFrustratedSolutions++;

				// WARNING: It is possible that the new solution has worse imbalance than it should, since
				//    the local ILS of one of the processes may not find the most efficient local solution.
				//    In these cases, discards the clustering solutions of each process and simply creates
				//    a zero-cost move solution (i.e. same imbalance than previous) based on the moved vertices
				//    and the previous best-known clustering.

				// => simply moves the cluster to the new process, with the same global imbalance value
				// evaluate if this zero-cost move improves the vertex balancing between processes
				if (labs(newSizeOfSourceCluster - newSizeOfDestCluster) <
						labs(sizeOfSourceCluster - sizeOfDestCluster)) {  // this zero-cost move is good

					// retrieves the initial global clustering as base for zero-cost
					bestGlobalClustering = initialGlobalClustering;

					// Swaps the clusters between the processes, objective function value remains the same
					bestGlobalClustering.setProcessOrigin(clusterToSwapB, sourceProcess);
					bestGlobalClustering.setProcessOrigin(clusterToSwapA, destinationProcess);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost swapCluster1opt improving solution found! I(P) = "
													<< bestGlobalClustering.getImbalance().getValue();

					// But must update the internal process imbalance
					ClusterArray bestClusterArray = bestGlobalClustering.getClusterArray();
					std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, bestClusterArray,
																				globalClusterArray, numberOfProcesses);
					bestGlobalClustering.setProcessImbalance(sourceProcess, processInternalImbalance[sourceProcess]);
					bestGlobalClustering.setProcessImbalance(destinationProcess, processInternalImbalance[destinationProcess]);
					bestClustering = bestGlobalClustering;

					// Must update the imbalance matrix between processes as well
					bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
					util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestGlobalClustering);
					foundBetterSolution = true;
				} else {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost swapCluster1opt move.";
					// TODO INCLUIR AQUI CHAMADA AO ROLLBACK DO MOVIMENTO DE VERTICES
					std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
					long N = g->getGlobalN();
					for(long i = 0; i < N; i++) {
						long k = currentSplitgraphClusterArray[i];
						verticesInProcess[k].push_back(i);
					}
					// moves the affected vertices between processes
					InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
												problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
												numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, 0, true, NULL, true,
												/*runILS=*/false, /*redistributeVertices=*/true);
					imsg.splitgraphClusterArray = currentSplitgraphClusterArray;
					for(int pr = 1; pr < np; pr++) {
						imsg.setVertexList(verticesInProcess[pr]);
						world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback message sent to process " << pr << " (runILS = false)";
					}
					imsg.setVertexList(verticesInProcess[0]);
					OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg, g);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback of master process done (runILS = false).";
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
					for(int pr = 1; pr < np; pr++) {
						OutputMessage omsg;
						mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
						int tag = stat.tag();
						int procNum = stat.source(), i = 0;
						if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
							BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
							i++;
							while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
								stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
								i++;
							}
							throw std::invalid_argument( "Error message received from slave. See error logs." );
						}
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback-done message received from process " << procNum << ".";
					}
				}
			}
		}
		if(foundBetterSolution) {
			break;
		}
	}
	return foundBetterSolution;
}

bool ImbalanceSubgraphParallelILS::twoMoveCluster(clusteringgraph::SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
		clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses,
		ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   2-Move  cluster  ***********************************
	// b) The master process tries to move 2 clusters from an overloaded process to 2 less-loaded processes.
	// select a 2 clusters to move - two clusters from an overloaded process (process with too many vertices)
	// (priotitize vertex load balancing between processes)
	// for now it is the cluster c that belongs to the processPair (x, y) where x is an overloaded process
	// Tries to move the top 25% of clusters of overloaded processes

	// *** STEP 1: INITIALIZATION AND INITIAL FULL DISTRIBUTED ILS
	bool foundBetterSolution = false;
	long n = g->getGlobalN();
	long nc = bestClustering.getNumberOfClusters();

	std::vector<Coordinate> overloadedProcessList = util.obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** twoMoveCluster";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph twoMoveCluster] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());

	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	// obtains the number of clusters from each process 'procNum'
	int np = bestSplitgraphClustering.getNumberOfClusters();

	clusteringgraph::Clustering bestGlobalClustering = bestClustering;
	const clusteringgraph::Clustering initialGlobalClustering = bestGlobalClustering;
	ProcessClustering initialSplitgraphClustering = bestSplitgraphClustering;
	const ClusterArray initialSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestGlobalClustering.getClusterArray();
	const ImbalanceMatrix initialProcessClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best global solution so far: I(P) = " << bestGlobalClustering.getImbalance().getValue();
	std::vector<unsigned int> clusterProcessOrigin = bestGlobalClustering.getClusterProcessOrigin();
	// ----------------------------------------------------------------------

	// *** STEP 2: POPULATE LIST OF MOVEMENTS
	// populates the list of movements (x, y) where x is the cluster that will be swapped to process y
	std::vector<Coordinate> movementListA, movementListB;
	// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] List of 2-move-cluster movements: ";
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2-move-cluster moving cluster from process P" << procSourceNum;

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = util.obtainListOfClustersFromProcess(*g, initialGlobalClustering, procSourceNum);
		if(bigClustersList.size() <= 2)  continue;  // the origin process must remain with at least 1 cluster
		// sorts the list of big clusters according to the number of vertices, descending order
		std::sort(bigClustersList.begin(), bigClustersList.end(), coordinate_ordering_desc());
		// TODO possible parametrization here! Parametrize the quantity of big clusters to be moved, two at a time
		// Will move the 2 biggest clusters (with more vertices) from process 'procSourceNum'
		for(int a = 0; a < bigClustersList.size(); a++) {
			long clusterToMoveA = bigClustersList[a].x;
			for(int b = 0; b < bigClustersList.size(); b++) {
				long clusterToMoveB = bigClustersList[b].x;
				if(clusterToMoveA < clusterToMoveB) {  // avoid repeating combinations of A and B
					std::vector<long> listOfMovedVerticesFromClusterA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveA);
					std::vector<long> listOfMovedVerticesFromClusterB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveB);
					int currentProcess = initialSplitgraphClusterArray[listOfMovedVerticesFromClusterA[0]];
					if(currentProcess != procSourceNum) {
						BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] Inconsistency detected between source process numbers.";
					}
					assert(currentProcess == procSourceNum);
					// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] 2-Move from process " << procSourceNum << ": Moving global clusters "
					//		<< clusterToMoveA << " and " <<	clusterToMoveB << ".";
					// LOAD BALANCING: Try to move the clusters to 2 processes with less vertices than a given threshold (less loaded process)
					ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
					long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0);
					for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {  // destination process 1
						if((procDestNum != procSourceNum) and (bestSplitgraphClustering.getClusterSize(procDestNum) +
								listOfMovedVerticesFromClusterA.size() < maxVerticesAllowedInProcess)) {
							// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Move global cluster " << clusterToMoveA
							//		<< " to process " << procDestNum;
							for(int procDestNum2 = 0; procDestNum2 < numberOfProcesses; procDestNum2++) {  // destination process 2
								if((procDestNum2 != procSourceNum) and (procDestNum2 != procDestNum)
										and (bestSplitgraphClustering.getClusterSize(procDestNum2) +
										listOfMovedVerticesFromClusterB.size() < maxVerticesAllowedInProcess)) {
									// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph]     And move global cluster " << clusterToMoveB
									//			<< " to process " << procDestNum2;
									Coordinate movementA(clusterToMoveA, procDestNum);
									movementListA.push_back(movementA);
									Coordinate movementB(clusterToMoveB, procDestNum2);
									movementListB.push_back(movementB);
								}
							}
						}
					}
				}
			}
		}
	}
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of potential 2-move-cluster movements: " << movementListA.size();
	// -------------------------------------------

	// *** STEP 3: EVALUATE SEVERAL MOVEMENTS IN PARALLEL (DISTRIBUTED ALGORITHM)
	// ==> Uses parallelism to evaluate more than one movement at the same time, according to the number of available processes
	int numberOfParallelEvaluations = 1; // (int)std::floor(numberOfSlaves / double(2.0));
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of cluster movements evaluated in parallel: " << numberOfParallelEvaluations;
	const ClusterArray currentSplitgraphClusterArray = initialSplitgraphClustering.getClusterArray();
	while(not movementListA.empty()) {
		// generates a list with the cluster-process movements that will be executed in parallel
		// all movements must involve the same source process
		std::vector<Coordinate> movementsInParallelA, movementsInParallelB;
		int currentSourceProcess = clusterProcessOrigin[movementListA.back().x];
		assert(clusterProcessOrigin[movementListA.back().x] == clusterProcessOrigin[movementListB.back().x]);
		for(int mov = 0; mov < numberOfParallelEvaluations; mov++) {
			int sourceProcess = clusterProcessOrigin[movementListA.back().x];
			// all movements must involve the same source process
			if(sourceProcess != currentSourceProcess) {
				break;
			}
			movementsInParallelA.push_back(movementListA.back());
			movementListA.pop_back();
			movementsInParallelB.push_back(movementListB.back());
			movementListB.pop_back();
			if(movementListA.empty())  break;
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Processing " << movementsInParallelA.size() << " movements in parallel...";

		// The leader will run ILS (after the following loop) to solve the move on the sourceProcess
		std::vector< long > verticesInSourceProcess;

		// Executes the movements inside the list movementsInParallel
		// STEP 3.1: sends the movements to each pair of processes through MPI
		std::vector<int> participatingProcessList;
		participatingProcessList.push_back(currentSourceProcess);
		std::vector<InputMessageParallelILS> messagesSent;
		std::map<int, OutputMessage> messageMap;
		int leaderParticipates = 0;
		int leaderProcess = 0;
		ClusterArray tempSplitgraphClusterArray;
		std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
		for(int mov = 0; mov < movementsInParallelA.size(); mov++) {

			long clusterToMoveA = movementsInParallelA[mov].x;
			int destinationProcessA = movementsInParallelA[mov].y;
			long clusterToMoveB = movementsInParallelB[mov].x;
			int destinationProcessB = movementsInParallelB[mov].y;

			// workerProcess1 must begin with process of rank 1
			int workerProcess1 = destinationProcessA;// 2 * (mov + 1) - 1;
			int workerProcess2 = destinationProcessB;// 2 * (mov + 1);
			participatingProcessList.push_back(workerProcess1);
			participatingProcessList.push_back(workerProcess2);

			std::vector<long> listOfMovedVerticesFromClusterA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveA);
			std::vector<long> listOfMovedVerticesFromClusterB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveB);
			tempSplitgraphClusterArray = initialSplitgraphClusterArray;
			// MOVE Step 1: Move the vertices from a specific cluster A from source process to dest process 1
			for(long elem = 0; elem < listOfMovedVerticesFromClusterA.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromClusterA[elem]] = destinationProcessA;
			}
			// MOVE Step 2: Move the vertices from a specific cluster B from source process to dest process 2
			for(long elem = 0; elem < listOfMovedVerticesFromClusterB.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromClusterB[elem]] = destinationProcessB;
			}

			// stores the new list of vertices from each process participating in the movement
			std::vector< long > verticesInDestinationProcessA, verticesInDestinationProcessB;
			for(long i = 0; i < n; i++) {
				long k = tempSplitgraphClusterArray[i];
				verticesInProcess[k].push_back(i);
				if(k == currentSourceProcess) {
					if(mov == 0) {  // avoids duplicating elements in this list
						verticesInSourceProcess.push_back(i);
					}
				} else if(k == destinationProcessA) {
					verticesInDestinationProcessA.push_back(i);
				} else if(k == destinationProcessB) {
					verticesInDestinationProcessB.push_back(i);
				}
			}
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to move global clusters " << clusterToMoveA << " and " << clusterToMoveB
					<< " from process " << currentSourceProcess << " to processes " << destinationProcessA << " and " << destinationProcessB;

			// 3 ILS procedures will be executed by 3 processes through MPI messages
			InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
								problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
								numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, -1, true, NULL, true,
								/*runILS=*/true, /*redistributeVertices=*/true);
			imsg.splitgraphClusterArray = tempSplitgraphClusterArray;
			// if(k < 0) {  imsg.setClustering(&CCclustering);  }
			InputMessageParallelILS imsg2 = imsg;
			InputMessageParallelILS imsg3 = imsg;
			// sends the modified subgraphs that will be solved by ILS
			imsg.setVertexList(verticesInDestinationProcessA);
			imsg2.setVertexList(verticesInDestinationProcessB);
			imsg3.setVertexList(verticesInSourceProcess);
			// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
			imsg.cudaEnabled = is_overloaded_process(verticesInDestinationProcessA);
			imsg2.cudaEnabled = is_overloaded_process(verticesInDestinationProcessB);
			imsg3.cudaEnabled = true;
			if(workerProcess1 != leaderProcess) {
				world.send(workerProcess1, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
			} else { leaderParticipates = 1; }
			if(workerProcess2 != leaderProcess) {
				world.send(workerProcess2, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg2);
			} else { leaderParticipates = 2; }
			if(currentSourceProcess != leaderProcess) {
				world.send(currentSourceProcess, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg3);
			} else { leaderParticipates = 3; }
			messagesSent.push_back(imsg);
			messagesSent.push_back(imsg2);
			messagesSent.push_back(imsg3);

			for(int pr = 1; pr < np; pr++) {
				if((pr != destinationProcessA) and (pr != destinationProcessB) and (pr != currentSourceProcess)) {
					InputMessageParallelILS imsg4 = imsg;
					imsg4.runILS = false;
					imsg4.setVertexList(verticesInProcess[pr]);
					world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg4);
					messagesSent.push_back(imsg4);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << pr << " (runILS = false)";
				}
			}
			participatingProcessList.push_back(workerProcess1);
			participatingProcessList.push_back(workerProcess2);
			participatingProcessList.push_back(currentSourceProcess);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to processes " << workerProcess1 << ", " << workerProcess2 << " and " << currentSourceProcess;
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] verticesInDestinationProcessA size = " << verticesInDestinationProcessA.size();
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] verticesInDestinationProcessB size = " << verticesInDestinationProcessB.size();
			// BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Size of ILS Input Message: " << (sizeof(imsg)/1024.0) << "kB.";
		}

		if(leaderParticipates) {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking local ILS on leader process (P0).";
			VariableNeighborhoodDescent localVND = *vnd;
			localVND.setTimeLimit(LOCAL_ILS_TIME_LIMIT);
			InputMessageParallelILS inputMsg = messagesSent[leaderParticipates - 1];
			// moves the affected vertices between processes
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(inputMsg, g);
			  // construct, &localVND, g, iter, iterMaxILS, perturbationLevelMax,
			  //		problem, info, verticesInSourceProcess);
			numberOfTestedCombinations += leaderProcessingMessage.numberOfTestedCombinations;
			messageMap[leaderProcess] = leaderProcessingMessage;
		} else {
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking redistribution on leader process (P0).";
			InputMessageParallelILS imsg4(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
								problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
								numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, -1, true, NULL, true,
								/*runILS=*/true, /*redistributeVertices=*/true);
			imsg4.splitgraphClusterArray = tempSplitgraphClusterArray;
			imsg4.runILS = false;
			imsg4.setVertexList(verticesInProcess[0]);
			OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg4, g);
		}

		// STEP 3.3: the leader receives the processing results
		// Maps the messages according to the movement executed
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
		long totalMessagesSent = messagesSent.size();
		if(leaderParticipates) {
			totalMessagesSent--;
		}
		for(int i = 0; i < totalMessagesSent; i++) {
			OutputMessage omsg;
			mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
			int tag = stat.tag();
			int procNum = stat.source();
			if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
				BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
				i++;
				while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
					stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
					i++;
				}
				throw std::invalid_argument( "Error message received from slave. See error logs." );
			}
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			// process the result of the execution of process i
			// sums the total number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// stores the time spent by each process
			// timeSpent[i+1] = omsg.timeSpent;
			messageMap[procNum] = omsg;
		}

		// STEP 3.4: For each movement, merge the partial solutions into a global solution for the whole graph
		for(int mov = 0; mov < movementsInParallelA.size(); mov++) {
			long clusterToMoveA = movementsInParallelA[mov].x;
			int destinationProcessA = movementsInParallelA[mov].y;
			long clusterToMoveB = movementsInParallelB[mov].x;
			int destinationProcessB = movementsInParallelB[mov].y;
			int sourceProcess = clusterProcessOrigin[clusterToMoveA];
			// workerProcess1 must begin with process of rank 1
			int workerProcess1 = destinationProcessA;// 2 * (mov + 1) - 1;
			int workerProcess2 = destinationProcessB;// 2 * (mov + 1);
			std::vector<long> listOfMovedVerticesFromClusterA = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveA);
			std::vector<long> listOfMovedVerticesFromClusterB = util.getListOfVeticesInCluster(*g, initialGlobalClustering, clusterToMoveB);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Analyzing movement pair " << mov;
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Move global clusters " << clusterToMoveA << " and " << clusterToMoveB
						<< " from process " << sourceProcess << " to processes " << destinationProcessA << " and " << destinationProcessB;

			// the new list containing the internal imbalance of each process local clustering
			std::vector<Imbalance> newInternalProcessImbalance = initialGlobalClustering.getInternalProcessImbalance();
			clusteringgraph::Clustering globalClustering;
			// removes all the clusters which belong to the processes participating in this cluster movement
			// and calculates the new number of clusters, to be used in the offset below
			clusteringgraph::Clustering tempClustering = initialGlobalClustering;
			tempClustering.removeAllClustersFromProcess(g, sourceProcess);
			tempClustering.removeAllClustersFromProcess(g, destinationProcessA);
			tempClustering.removeAllClustersFromProcess(g, destinationProcessB);
			std::vector<unsigned int> newClusterProcessOrigin = tempClustering.getClusterProcessOrigin();
			long clusterOffset = tempClustering.getNumberOfClusters();
			ClusterArray globalClusterArray = tempClustering.getClusterArray();
			ClusterArray previousSplitgraphClusterArray = currentSplitgraphClusterArray;
			ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;

			ImbalanceMatrix tempProcessClusterImbMatrix = initialProcessClusterImbMatrix;
			// MOVE Step 1: Move the vertices from a specific cluster A from source process to dest process 1
			for(long elem = 0; elem < listOfMovedVerticesFromClusterA.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromClusterA[elem]] = destinationProcessA;
			}
			// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processC
			util.updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
					listOfMovedVerticesFromClusterA, tempProcessClusterImbMatrix, numberOfSlaves + 1);
			previousSplitgraphClusterArray = tempSplitgraphClusterArray;

			// MOVE Step 2: Move the vertices from a specific cluster B from source process to dest process 2
			for(long elem = 0; elem < listOfMovedVerticesFromClusterB.size(); elem++) {
				tempSplitgraphClusterArray[listOfMovedVerticesFromClusterB[elem]] = destinationProcessB;
			}
			// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix process
			util.updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
					listOfMovedVerticesFromClusterB, tempProcessClusterImbMatrix, numberOfSlaves + 1);
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from 2-move done.";

			// Valida se a matriz incremental e a full sao iguais
			/*
			ImbalanceMatrix tempProcessClusterImbMatrix3 = util.calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray,
					this->vertexImbalance, numberOfSlaves + 1);
			bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix3.pos, 0.001, 0.1);
			bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix3.neg, 0.001, 0.1);
			BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

			// Merges the partial solutions into a global solution for the whole graph
			std::vector<OutputMessage> sourceAndDestinationMovementMsg;
			sourceAndDestinationMovementMsg.push_back(messageMap[sourceProcess]);
			sourceAndDestinationMovementMsg.push_back(messageMap[workerProcess1]);
			sourceAndDestinationMovementMsg.push_back(messageMap[workerProcess2]);
			// sets the new imbalance for cluster-move source process (comes from leader process message - rank 0)
			newInternalProcessImbalance[sourceProcess] = messageMap[sourceProcess].clustering.getImbalance();
			long leader_nc = messageMap[sourceProcess].clustering.getNumberOfClusters();
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph P" << sourceProcess << " (source): num_edges = " << messageMap[sourceProcess].num_edges
						<< " , num_vertices = " << messageMap[sourceProcess].num_vertices	<< ", I(P) = " << messageMap[sourceProcess].clustering.getImbalance().getValue() << ", k = " << leader_nc;

			// sets the new imbalance for cluster-move destination process A
			newInternalProcessImbalance[destinationProcessA] = messageMap[workerProcess1].clustering.getImbalance();
			// sets the new imbalance for cluster-move destination process B
			newInternalProcessImbalance[destinationProcessB] = messageMap[workerProcess2].clustering.getImbalance();

			for(int proc = 0; proc < sourceAndDestinationMovementMsg.size(); proc++) {
				OutputMessage msg = sourceAndDestinationMovementMsg[proc];
				ClusterArray localClusterArray = msg.clustering.getClusterArray();

				long msg_nc = msg.clustering.getNumberOfClusters();
				if(msg_nc > 0) {
					assert(localClusterArray.size() == msg.globalVertexId.size());
					for(long v = 0; v < localClusterArray.size(); v++) {
						// Obtains vertex v's number in the global graph
						long vglobal = msg.globalVertexId[v];
						globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
					}
					clusterOffset += msg_nc;
					// all the clusters in this interval belong either to the source or destination process in the swap movement
					for(long c = 0; c < msg_nc; c++) {
						if(proc == 0) {
							newClusterProcessOrigin.push_back(sourceProcess);
						} else if(proc == 1) {
							newClusterProcessOrigin.push_back(destinationProcessA);
						} else {
							newClusterProcessOrigin.push_back(destinationProcessB);
						}
					}
				}
				if(proc > 0) {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph P" << participatingProcessList[proc + 1] << ": num_edges = " << msg.num_edges
						<< " , num_vertices = " << msg.num_vertices	<< ", I(P) = " << msg.clustering.getImbalance().getValue() << ", k = " << msg_nc;
				}
			}

			// Calculates the internal imbalance sum (inside each process)
			Imbalance internalImbalance = util.calculateInternalImbalanceSumOfAllProcesses(newInternalProcessImbalance);
			// Calculates the external imbalance sum (between processes)
			Imbalance externalImbalance = util.calculateExternalImbalanceSumBetweenProcesses(tempProcessClusterImbMatrix);

			// 5. Builds the clustering with the merge of each process local ILS result
			// EVITA O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
			globalClustering = clusteringgraph::Clustering(globalClusterArray, *g, problem, internalImbalance.getPositiveValue() + externalImbalance.getPositiveValue(),
					internalImbalance.getNegativeValue() + externalImbalance.getNegativeValue(), newClusterProcessOrigin, newInternalProcessImbalance);

			/*
			clusteringgraph::Clustering validation(globalClusterArray, *g, problem);
			if(fabs(validation.getImbalance().getValue() - globalClustering.getImbalance().getValue()) > EPS) {
				BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] I(P) DOES NOT MATCH! Correct I(P) = "
						<< std::setprecision(5) << std::fixed << validation.getImbalance().getValue()
						<< " vs I(P) = " << globalClustering.getImbalance().getValue();

				stringstream ss;
				ss << "InternalImbalanceVector: ";
				for(int x = 0; x < numberOfProcesses; x++) {
					ss << newInternalProcessImbalance[x].getValue() << " ";
				}
				BOOST_LOG_TRIVIAL(info) << ss.str();

				std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, tempSplitgraphClusterArray,
						globalClusterArray, numberOfProcesses);
			} */

			if(globalClustering.getImbalance().getValue() < bestGlobalClustering.getImbalance().getValue()) {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2-move-Cluster Improved solution found! I(P) = "
						<< globalClustering.getImbalance().getValue();

				// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
				bestClustering = globalClustering;
				bestGlobalClustering = globalClustering;

				// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
				bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
				util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestGlobalClustering);
				foundBetterSolution = true;
				/*  FIRST IMPROVEMENT ON PARALLEL DISTRIBUTED ILS EVALUATIONS IS DISABLED!
				break;
				*/
			} else {
				long sizeOfSourceProcess = bestSplitgraphClustering.getClusterSize(sourceProcess);
				long sizeOfDestProcessA = bestSplitgraphClustering.getClusterSize(destinationProcessA);
				long sizeOfDestProcessB = bestSplitgraphClustering.getClusterSize(destinationProcessB);
				long newSizeOfSourceProcess = sizeOfSourceProcess - listOfMovedVerticesFromClusterA.size() - listOfMovedVerticesFromClusterB.size();
				long newSizeOfDestProcessA = sizeOfDestProcessA + listOfMovedVerticesFromClusterA.size();
				long newSizeOfDestProcessB = sizeOfDestProcessB + listOfMovedVerticesFromClusterB.size();

				if (globalClustering.getImbalance().getValue() == bestGlobalClustering.getImbalance().getValue()) {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2-move Zero-cost move.";
				} else {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2-move Worse solution found (Obj = "
							<< globalClustering.getImbalance().getValue() << "). Fixing to zero-cost move.";
				}
				numberOfFrustratedSolutions++;

				// WARNING: It is possible that the new solution has worse imbalance than it should, since
				//    the local ILS of one of the processes may not find the most efficient local solution.
				//    In these cases, discards the clustering solutions of each process and simply creates
				//    a zero-cost move solution (i.e. same imbalance than previous) based on the moved vertices
				//    and the previous best-known clustering.

				// => simply moves the clusters to the new process, with the same global imbalance value
				// evaluate if this zero-cost move improves the vertex balancing between processes
				if (labs(newSizeOfSourceProcess - (sizeOfDestProcessA + sizeOfDestProcessB)) <
					labs(sizeOfSourceProcess - (newSizeOfDestProcessA + newSizeOfDestProcessB))) {  // this zero-cost move is good

					// retrieves the initial global clustering as base for zero-cost
					bestGlobalClustering = initialGlobalClustering;

					// Moves the clusters to the 2 destination processes
					bestGlobalClustering.setProcessOrigin(clusterToMoveA, destinationProcessA);
					bestGlobalClustering.setProcessOrigin(clusterToMoveB, destinationProcessB);

					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost twoMoveCluster improving solution found! I(P) = "
													<< bestGlobalClustering.getImbalance().getValue();

					// But must update the internal process imbalance
					ClusterArray bestClusterArray = bestClustering.getClusterArray();
					std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, bestClusterArray,
																				globalClusterArray, numberOfProcesses);
					bestGlobalClustering.setProcessImbalance(sourceProcess, processInternalImbalance[sourceProcess]);
					bestGlobalClustering.setProcessImbalance(destinationProcessA, processInternalImbalance[destinationProcessA]);
					bestGlobalClustering.setProcessImbalance(destinationProcessB, processInternalImbalance[destinationProcessB]);

					bestClustering = bestGlobalClustering;

					// Must update the imbalance matrix between processes as well
					bestSplitgraphClustering = ProcessClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
					util.validaSplitgraphArray(*g, bestSplitgraphClustering, bestGlobalClustering);
					foundBetterSolution = true;
				} else {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost twoMoveCluster Movement rejected.";
					// TODO INCLUIR AQUI CHAMADA AO ROLLBACK DO MOVIMENTO DE VERTICES
					std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
					long N = g->getGlobalN();
					for(long i = 0; i < N; i++) {
						long k = currentSplitgraphClusterArray[i];
						verticesInProcess[k].push_back(i);
					}
					// moves the affected vertices between processes
					InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
												problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, LOCAL_ILS_TIME_LIMIT,
												numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, 0, true, NULL, true,
												/*runILS=*/false, /*redistributeVertices=*/true);
					imsg.splitgraphClusterArray = currentSplitgraphClusterArray;
					for(int pr = 1; pr < np; pr++) {
						imsg.setVertexList(verticesInProcess[pr]);
						world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback message sent to process " << pr << " (runILS = false).";
					}
					imsg.setVertexList(verticesInProcess[0]);
					OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg, g);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback of master process done (runILS = false).";
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages...";
					for(int pr = 1; pr < np; pr++) {
						OutputMessage omsg;
						mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
						int tag = stat.tag();
						int procNum = stat.source(), i = 0;
						if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
							BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
							i++;
							while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
								stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
								i++;
							}
							throw std::invalid_argument( "Error message received from slave. See error logs." );
						}
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rollback-done message received from process " << procNum << ".";
					}
				}
			}
		}
		if(foundBetterSolution) {
			break;
		}
	}
	return foundBetterSolution;
}

long ImbalanceSubgraphParallelILS::variableNeighborhoodDescent(clusteringgraph::SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
		clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses,
		ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info, const double& timeSpentSoFar, int invocationNumber) {

	int r = 1, iteration = 0, l = 3;
	double timeSpentOnLocalSearch = 0.0;
	long improvedOnVND = 0;
	BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Global VND local search (Invocation " << invocationNumber << ") ...";
	std::vector<long> improvementStats(l + 1, 0);
	std::vector<long> neigExecStats(l + 1, 0);
	// stores the list of improved solutions found at each VND iteration (global clustering and splitgraph clustering)
	std::vector< std::pair<clusteringgraph::Clustering, ProcessClustering> > solutionHistory;
	solutionHistory.push_back( std::make_pair(bestClustering, bestSplitgraphClustering) );
	// stores the process balancing (number of clusters, vertices) info and solution values throughout the time
	std::vector<string> splitgraphVNDProcessBalancingHistory;
	stringstream ss;
	ss << iteration++ << ", " << timeSpentOnLocalSearch << ", " << bestClustering.getImbalance().getValue();
	splitgraphVNDProcessBalancingHistory.push_back(ss.str());
	numberOfFrustratedSolutions = 0;

	// runs a full execution of distributed ILS, in order to find the internal imbalance value of each process
	ProcessClustering bestSplitgraphClusteringVND = bestSplitgraphClustering;
	clusteringgraph::Clustering bestClusteringVND = bestClustering;
	bestClusteringVND = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem,
								info, bestSplitgraphClusteringVND, bestClusteringVND);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best global solution so far: I(P) = " << bestClustering.getImbalance().getValue();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best VND solution so far: I(P) = " << bestClusteringVND.getImbalance().getValue();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Initial solution built.";
	int global_vnd_loop_count = 1;

	while (r <= l && (timeSpentSoFar + timeSpentOnLocalSearch < vnd->getTimeLimit())) {
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();
		BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Global VND neighborhood (count = " << global_vnd_loop_count++ << ") r = " << r << " ...";

		// rebalanceClustersBetweenProcessesWithZeroCost(g, problem, bestSplitgraphClusteringVND, bestClusteringVND, numberOfProcesses);

		bool improved = false;
		if(r == 1) {
			// ************************   Move   1 - opt   cluster ***********************************
			improved = moveCluster1opt(g, bestSplitgraphClusteringVND, bestClusteringVND,
									numberOfProcesses, construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} else if(r == 4) {
			// ************************   Split  Cluster  Move ***********************************
			/*  DISABLED IN PARALLEL GRAPH VERSION!
			improved = splitClusterMove(g, bestSplitgraphClusteringVND, bestClusteringVND,
									numberOfProcesses, construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
			*/
		} else if(r == 3) {  // TODO return 2-move-cluster to r = 2
			// ************************   2-Move   cluster ***********************************
			improved = twoMoveCluster(g, bestSplitgraphClusteringVND, bestClusteringVND,
									numberOfProcesses, construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} else if(r == 2) {
			// ************************   Swap   1 - opt   cluster ***********************************
			improved = swapCluster1opt(g, bestSplitgraphClusteringVND, bestClusteringVND,
									numberOfProcesses, construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		}
		// records statistics for neighborhood structures execution
		neigExecStats[r]++;
		if(improved)  improvementStats[r]++;

		// => Finally: Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timeSpentOnLocalSearch += (end_time.wall - start_time.wall)	/ double(1000000000);

		if(bestClusteringVND.getImbalance().getValue() < 0.0) {
			BOOST_LOG_TRIVIAL(error)<< "[Splitgraph] Objective function below zero. Obj = " << bestClusteringVND.getImbalance().getValue();
			break;
		}
		if(improved) {
			// BOOST_LOG_TRIVIAL(trace) << myRank << ": New local solution found: " << setprecision(2) << il.getValue() << endl;
			r = 1;
			improvedOnVND++;
			solutionHistory.push_back( std::make_pair(bestClusteringVND, bestSplitgraphClusteringVND) );
			stringstream ss;
			ss << iteration << ", " << timeSpentOnLocalSearch << ", " << bestClusteringVND.getImbalance().getValue();
			splitgraphVNDProcessBalancingHistory.push_back(ss.str());

			// prints the new process x vertex configuration
			int np = bestSplitgraphClusteringVND.getNumberOfClusters();
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best clustering config so far:";
			for(int px = 0; px < np; px++) {
				int k = util.obtainListOfClustersFromProcess(*g, bestClusteringVND, px).size();
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << px << ": num_vertices = "
						<< bestSplitgraphClusteringVND.getClusterSize(px) << ", k = " << k;
			}

			if(bestClusteringVND.getImbalance().getValue() <= EPS) {
				BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Stopping global VND local search, since I(P) <= 0!";
				break;
			}
		} else {  // no better result found in neighborhood
			r++;
			// BOOST_LOG_TRIVIAL(debug) << "Changed to neighborhood size l = " << k;
		}
		iteration++;
	}
	BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Global VND local search done. Obj = " << bestClusteringVND.getImbalance().getValue() <<
			". Time spent: " << timeSpentOnLocalSearch << " s.";

	stringstream splitgraphVNDResults;
	BOOST_LOG_TRIVIAL(info) << "[Splitgraph] Number of frustrated ILS solutions: " << numberOfFrustratedSolutions;
	splitgraphVNDResults << "Number of frustrated ILS solutions, " << numberOfFrustratedSolutions << "\n";
	BOOST_LOG_TRIVIAL(info) << "[Splitgraph] Execution and Improvement statistics for each neighborhood: ";
	splitgraphVNDResults << "Execution and Improvement statistics for each neighborhood: \n";
	BOOST_LOG_TRIVIAL(info) << "Neighborhood size, Num exec, Num improv";
	splitgraphVNDResults << "Neighborhood size, Num exec, Num improv\n";
	for(int i = 1; i <= l; i++) {
		BOOST_LOG_TRIVIAL(info) << "r = " << i << ", " << neigExecStats[i] << ", " << improvementStats[i];
		splitgraphVNDResults << "r = " << i << ", " << neigExecStats[i] << ", " << improvementStats[i] << "\n";
	}
	splitgraphVNDResults << "Best solution found \n";

	std::vector< std::vector< long > > verticesInProcess(numberOfProcesses, std::vector< long >());
	const ClusterArray splitgraphClusterArray = bestSplitgraphClusteringVND.getClusterArray();
	long n = g->getGlobalN();
	for(long i = 0; i < n; i++) {
		long k = splitgraphClusterArray[i];
		verticesInProcess[k].push_back(i);
	}
	splitgraphVNDResults << "Subgraph, n, m, k \n";
	for(int p = 0; p < numberOfProcesses; p++) {
		// BGLSignedGraph sg(*(g->graph), verticesInProcess[p]);
		// TODO TROCAR PELO ParallelBGLSignedGraph ou equivalente contendo as informacoes da divisao de vertices entre os processos
		std::vector<Coordinate> clusterList = util.obtainListOfClustersFromProcess(*g, bestClusteringVND, p);

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_vertices = " << num_vertices(*(g->graph)) <<
				"; num_edges = " <<	num_edges(*(g->graph)) << "; k = " << clusterList.size();
		splitgraphVNDResults << p << ", " << num_vertices(*(g->graph)) << ", " << num_edges(*(g->graph)) << ", " << clusterList.size() << "\n";
	}
	splitgraphVNDResults << "I(P), " << bestClusteringVND.getImbalance().getValue() << "\n";
	splitgraphVNDResults << "Time spent, " << timeSpentOnLocalSearch << "\n";

	// Saves the splitgraph statistics to csv file
	stringstream filePrefix;
	filePrefix << "VND";
	filePrefix << (invocationNumber);
	filePrefix << "-statistics-splitgraph";
	generateOutputFile(problem, splitgraphVNDResults, info.outputFolder, info.fileId, info.executionId,
			info.processRank, filePrefix.str(), construct->getAlpha(), l, iter);

	// Exports the splitgraph solutions to csv file
	stringstream sshistory;
	sshistory << "Iteration, Time, I(P)";
	for(int i = 1; i <= numberOfProcesses; i++) {
		sshistory << ", e" << i << ", n" << i << ", k" << i;
	}
	sshistory << "\n";
	for(int iternum = 0; iternum < solutionHistory.size(); iternum++) {
		sshistory << splitgraphVNDProcessBalancingHistory[iternum];

		clusteringgraph::Clustering globalClustering = solutionHistory[iternum].first;
		ProcessClustering splitgraphClustering = solutionHistory[iternum].second;
		std::vector< std::vector< long > > verticesInProcess(numberOfProcesses, std::vector< long >());
		ClusterArray splitgraphClusterArray = splitgraphClustering.getClusterArray();
		for(long i = 0; i < n; i++) {
			long k = splitgraphClusterArray[i];
			verticesInProcess[k].push_back(i);
		}
		for(int p = 0; p < numberOfProcesses; p++) {
			// clusteringgraph::SignedGraph sg(*(g->graph), verticesInProcess[p]);
			// TODO TROCAR POR FUNCAO EQUIVALENTE DA PARALLEL BGL QUE FORNECA ESSES DADOS: NUMERO DE VERTICES E ARESTAS
			std::vector<Coordinate> clusterList = util.obtainListOfClustersFromProcess(*g, globalClustering, p);
			sshistory << ", " << num_edges(*(g->graph)) << ", " << num_vertices(*(g->graph)) << ", " << clusterList.size();
		}
		sshistory << "\n";
	}
	// includes the biggest cluster size of each process in the best solution
	clusteringgraph::Clustering globalClustering = solutionHistory.back().first;
	ProcessClustering splitgraphClustering = solutionHistory.back().second;
	for(int i = 1; i <= numberOfProcesses; i++) {
		sshistory << "verticesIn(P" << i << "), " << "biggestClusterSize(P" << i << "), ";
	}
	sshistory << "\n";
	for(int p = 0; p < numberOfProcesses; p++) {
		// BGLSignedGraph sg(*(g->graph), verticesInProcess[p]);
		// TODO TROCAR PELO ParallelBGLSignedGraph ou equivalente contendo as informacoes da divisao de vertices entre os processos
		std::vector<Coordinate> bigClustersList = util.obtainListOfClustersFromProcess(*g, globalClustering, p);
		std::vector<Coordinate>::iterator biggestCluster = std::max_element(bigClustersList.begin(), bigClustersList.end(), coordinate_ordering_asc());
		sshistory << num_vertices(*(g->graph)) << ", " << boost::math::lround(biggestCluster->value) << ", ";
	}
	sshistory << "\n";
	// includes the number of frustrated solutions
	sshistory << "Frustrated ILS solutions, " << numberOfFrustratedSolutions << "\n";

	// Saves the splitgraph solution history to csv file
	stringstream filePrefix2;
	filePrefix2 << "VND";
	filePrefix2 << (invocationNumber);
	filePrefix2 << "-solutionHistory-splitgraph";
	generateOutputFile(problem, sshistory, info.outputFolder, info.fileId, info.executionId,
			info.processRank, filePrefix2.str(), construct->getAlpha(), l, iter);

	if( (fabs(bestClusteringVND.getImbalance().getValue() - bestClustering.getImbalance().getValue()) > EPS)
			and (bestClusteringVND.getImbalance().getValue() - bestClustering.getImbalance().getValue() < EPS) ) {  // a < b
		bestClustering = bestClusteringVND;
		bestSplitgraphClustering = bestSplitgraphClusteringVND;
		return improvedOnVND;
	} else {
		return 0;
	}
}

/**
  *  Identifies a pseudo clique C+ inside an overloaded cluster, and tries to
  *  move this clique to another (not overloaded) process as a new cluster.
*/
bool ImbalanceSubgraphParallelILS::splitClusterMove(clusteringgraph::SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
		clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses,
		ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		problem::ClusteringProblem& problem, ExecutionInfo& info) {

	const ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	const ImbalanceMatrix processClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
	std::vector<Coordinate> overloadedProcessList = util.obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** splitClusterMove";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph splitClusterMove] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	bool foundBetterSolution = false;
	long n = g->getGlobalN();
	long nc = bestClustering.getNumberOfClusters();
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number
		BOOST_LOG_TRIVIAL(info) << "SplitClusterMove for overloaded process number " << procSourceNum <<
				" with " << prIter->value << " vertices...";

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = util.obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
		// sorts the list of big clusters according to the number of vertices, descending order
		std::sort(bigClustersList.begin(), bigClustersList.end(), coordinate_ordering_desc());
		// TODO possible parametrization here! Parametrize the quantity of big clusters to be moved, one at a time
		long numberOfClustersToBeInspected = (int)ceil(bigClustersList.size() * PERCENTAGE_OF_MOST_IMBALANCED_CLUSTERS_TO_BE_MOVED);
		if(numberOfClustersToBeInspected == 0) {
			numberOfClustersToBeInspected = 1;
		}
		BOOST_LOG_TRIVIAL(info) << "Will invoke splitClusterMove for " << numberOfClustersToBeInspected << " clusters...";
		for(long clusterCount = 0; clusterCount < numberOfClustersToBeInspected; clusterCount++) {
			long clusterX = bigClustersList[clusterCount].x;
			string strCl = boost::lexical_cast<std::string>(clusterX);
			BOOST_LOG_TRIVIAL(info) << "Invoking splitClusterMove for spit graph global cluster number = " << strCl << "...";

			// use a fast heuristic to find a positive clique C+ in cluster X
			std::vector<long> cliqueC = util.findPseudoCliqueC(g, bestClustering, clusterX);

			// try to move the entire clique C+ (as a new cluster) to another process B
			// as long as process B is not overloaded
			// finds out to which process the cluster belongs to
			int currentProcess = splitgraphClusterArray[cliqueC.front()];
			if(currentProcess != procSourceNum) {
				BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] Inconsistency detected between source process numbers.";
			}
			assert(currentProcess == procSourceNum);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Move pseudo clique of size " << cliqueC.size()
					<<  " from process " << procSourceNum <<
					": Moving from global cluster " << clusterX << " of size " << bestClustering.getClusterSize(clusterX);
			assert(bigClustersList[clusterCount].y == procSourceNum);

			// LOAD BALANCING: Try to move the clique to a process with less vertices than a given threshold (less loaded process)
			const ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
			// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
			// TODO possible parametrization here!
			long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0);
			for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {
				// only makes the swap if the max vertex threshold is respected
				if((procDestNum != procSourceNum) and
						(bestSplitgraphClustering.getClusterSize(procDestNum) + cliqueC.size() < bestSplitgraphClustering.getClusterSize(procSourceNum))) {

					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to move the pseudo clique to process " << procDestNum;
					ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
					// Realiza a movimentacao dos vertices de um cluster especifico (cluster move)
					for(std::vector<long>::iterator cliqueIter = cliqueC.begin(); cliqueIter != cliqueC.end(); cliqueIter++) {
						long vertex = *cliqueIter;
						tempSplitgraphClusterArray[vertex] = procDestNum;
					}
					// O RECALCULO DA MATRIZ ABAIXO EH FEITO DE FORMA INCREMENTAL, reutilizando a matrix processClusterImbMatrix
					// ImbalanceMatrix tempProcessClusterImbMatrix2 = util.calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray,
					// 		this->vertexImbalance, numberOfSlaves + 1);
					ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;
					// Realiza o calculo do delta da matriz de imbalance entre processos
					util.updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
							cliqueC, tempProcessClusterImbMatrix, numberOfSlaves + 1);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";
					ProcessClustering tempSplitgraphClustering(g, problem, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
					// Valida se a matriz incremental e a full sao iguais
					/*
					bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
					bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
					BOOST_LOG_TRIVIAL(info) << "*** As matrizes de delta sao iguais: " << igual << " e " << igual2;
					*/

					// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
					clusteringgraph::Clustering newGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem,
							info, tempSplitgraphClustering, bestClustering);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
					/*
					ClusterArray cTemp = newGlobalClustering.getClusterArray();
					clusteringgraph::Clustering validation(cTemp, *g, problem);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Full Obj Calc: I(P) = " << validation.getImbalance().getValue();
					if(fabs(validation.getImbalance().getValue() - globalClustering.getImbalance().getValue()) > EPS) {
						BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] I(P) DOES NOT MATCH! Correct I(P) = "
							<< std::setprecision(5) << std::fixed << validation.getImbalance().getValue()
							<< " vs I(P) = " << globalClustering.getImbalance().getValue();
					} */

					if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] splitClusterMove Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
						bestClustering = newGlobalClustering;

						bestSplitgraphClustering = tempSplitgraphClustering;
						foundBetterSolution = true;
						break;
					} else {
						long sizeOfSourceProcess = bestSplitgraphClustering.getClusterSize(procSourceNum);
						long sizeOfDestProcess = bestSplitgraphClustering.getClusterSize(procDestNum);
						long newSizeOfSourceProcess = sizeOfSourceProcess - cliqueC.size();
						long newSizeOfDestProcess = sizeOfDestProcess + cliqueC.size();
						// cannot count the number of frustrated ILS solutions here, since split-cluster move is not a zero-cost move
						// it can really make global imbalance value worse

						if (newGlobalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] No Zero-cost move for pseudo-clique move.";
							// WARNING: It is possible that the new solution has worse imbalance than it should, since
							//    the local ILS of one of the processes may not find the most efficient local solution.
							//    In these cases, discards the clustering solutions of each process and simply creates
							//    a zero-cost move solution (i.e. same imbalance than previous) based on the moved vertices
							//    and the previous best-known clustering.

							// => simply swaps the clusters between the processes, keeping the same global imbalance value
							// evaluate if this zero-cost move improves the vertex balancing between processes
							/* DISABLED, SINCE EXTRACTING A PSEUDOCLIQUE FROM AN EXISTING CLUSTER DOES NOT GENERATE A ZERO-COST MOVE
							if (labs(newSizeOfSourceProcess - newSizeOfDestProcess) <
									labs(sizeOfSourceProcess - sizeOfDestProcess)) {  // this zero-cost move is good

								// Update the cluster process origin, objective function value remains the same
								newGlobalClustering.setProcessOrigin(clusterToSwapB, sourceProcess);
								newGlobalClustering.setProcessOrigin(clusterToSwapA, destinationProcess);

								// But must update the internal process imbalance
								newGlobalClustering.setProcessImbalance(sourceProcess,
										util.calculateProcessInternalImbalance(g, newGlobalClustering, sourceProcess));
								newGlobalClustering.setProcessImbalance(destinationProcess,
										util.calculateProcessInternalImbalance(g, newGlobalClustering, destinationProcess));

								bestClustering = newGlobalClustering;
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost splitClusterMove improving solution found! I(P) = "
																<< bestClustering.getImbalance().getValue();

								// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
								processClusterImbMatrix = tempProcessClusterImbMatrix;

								bestSplitgraphClustering = clusteringgraph::Clustering(tempSplitgraphClusterArray, *g, problem);
								foundBetterSolution = true;
								break;
							} else {
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost splitClusterMove move.";
							} */
						} else {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Worse solution found. Nothing to do.";
						}
					}
				} else {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Clique movement rejected.";
				}
			}
			if(foundBetterSolution) {
				break;
			}
		}
		if(foundBetterSolution) {
			break;
		}
	}
	return foundBetterSolution;
}

void ImbalanceSubgraphParallelILS::rebalanceClustersBetweenProcessesWithZeroCost(clusteringgraph::SignedGraph* g, problem::ClusteringProblem& problem,
		ProcessClustering& bestSplitgraphClustering, clusteringgraph::Clustering& bestClustering, const int& numberOfProcesses) {

	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	bool foundBetterSolution = false;
	long n = g->getGlobalN();
	long nc = bestClustering.getNumberOfClusters();
	ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	RandomUtil randomUtil;

	// A vertex-overloaded process is a process with more than (2* n / numberOfProcesses) vertices.
	long numberOfMaxVertices = (long)ceil(2 * n / (double)numberOfProcesses);
	std::vector<Coordinate> overloadedProcessList = util.obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering, numberOfMaxVertices);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** rebalanceClustersBetweenProcessesWithZeroCost";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph twoMoveCluster] Found " << overloadedProcessList.size()  << " overloaded processes (with more than "
			<< numberOfMaxVertices << " vertices.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number
		long numberOfVerticesInProcess = bestSplitgraphClustering.getClusterSize(procSourceNum);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Processing overloaded process number " << procSourceNum
				<< " that has " << numberOfVerticesInProcess << " vertices.";

		// LOAD BALANCING: Try to move the clusters to other processes with less vertices than a given threshold (less loaded process)
		// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
		// TODO possible parametrization here!
		int numberOfProcesses = bestSplitgraphClustering.getNumberOfClusters();
		long maxVerticesAllowedInProcess = boost::math::lround(n / (double)numberOfProcesses);
		std::vector<long> clustersToMove;
		std::vector<Coordinate> fullClusterList = util.obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
		// sorts the list of clusters according to the number of vertices, ascending order
		std::sort(fullClusterList.begin(), fullClusterList.end(), coordinate_ordering_asc());
		long numberOfVerticesToMove = 0;

		for(std::vector<Coordinate>::iterator iter = fullClusterList.begin(); iter != fullClusterList.end(); iter++) {
			long k = iter->x;
			long clusterSize = bestClustering.getClusterSize(k);
			if((clusterSize >= maxVerticesAllowedInProcess) or (fullClusterList.size() - clustersToMove.size() <= 1) or
					(numberOfVerticesInProcess - numberOfVerticesToMove < maxVerticesAllowedInProcess)) {
				break;
			}
			clustersToMove.push_back(k);
			numberOfVerticesToMove += clusterSize;
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Moving " << numberOfVerticesToMove
									<< " vertices from process " << procSourceNum;
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] After move, process " << procSourceNum << " will remain with "
				<< (numberOfVerticesInProcess - numberOfVerticesToMove) << " vertices and "
				<< (fullClusterList.size() - clustersToMove.size()) << " cluster(s).";
		ClusterArray previousSplitgraphClusterArray = currentSplitgraphClusterArray;

		for(std::vector<long>::iterator iter = clustersToMove.begin(); iter != clustersToMove.end(); iter++) {
			// finds a cluster in procSourceNum that can be moved to another process
			long clusterA = *iter;
			std::vector<long> listOfMovedVerticesFromClusterA = util.getListOfVeticesInCluster(*g, bestClustering, clusterA);

			// moves the cluster to a potential destination process in a randomized way
			for (int procDestNum = randomUtil.next(0, numberOfProcesses - 1), cont = 0; cont < numberOfProcesses; cont++) {  // destination process
			// for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {
				if((procDestNum != procSourceNum) and (bestSplitgraphClustering.getClusterSize(procDestNum) +
						listOfMovedVerticesFromClusterA.size() < maxVerticesAllowedInProcess)) {
					 BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Moving global cluster " << clusterA << " to process " << procDestNum;

					// MOVE Step 1: Move the vertices from a specific cluster A from source process to dest process 1
					for(long elem = 0; elem < listOfMovedVerticesFromClusterA.size(); elem++) {
						currentSplitgraphClusterArray[listOfMovedVerticesFromClusterA[elem]] = procDestNum;
					}

					// Moves the clusters to the destination process -> MUST recalculate processes internal imbalance in the end of the loop
					moveClusterToDestinationProcessZeroCost(g, bestClustering, bestSplitgraphClustering, clusterA,
							procSourceNum, procDestNum);

					// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX
					const ImbalanceMatrix processClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
					ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;
					// Realiza o calculo do delta da matriz de imbalance entre processos
					util.updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, currentSplitgraphClusterArray,
							listOfMovedVerticesFromClusterA, tempProcessClusterImbMatrix, numberOfSlaves + 1);
					previousSplitgraphClusterArray = currentSplitgraphClusterArray;
					// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";
					// Valida se a matriz incremental e a full sao iguais
					/*
					ImbalanceMatrix tempProcessClusterImbMatrix2 = util.calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray,
						this->vertexImbalance, numberOfSlaves + 1);
					bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
					bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
					BOOST_LOG_TRIVIAL(info) << "Op1 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

					// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from move done.";
					// Assumes the imbalance matrix between processes has already been updated in the code above
					bestSplitgraphClustering = ProcessClustering(g, problem, currentSplitgraphClusterArray, tempProcessClusterImbMatrix);

					// the movement for this cluster is done, proceed to the next one
					break;
				}
				// increment rule
				procDestNum++;
				if(procDestNum >= numberOfProcesses) {
					procDestNum = 0;  // TODO VALIDAR SE O NUMERO DO PRIMEIRO PROCESSO COM PARALLEL GRAPH COMECA EM 1 MESMO
				}
			}  // next destination process to be tested
		}  // next cluster to be moved
	}  // next overloaded processor
	// recalculates the process-process imbalance matrix
	// processClusterImbMatrix = util.calculateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, this->vertexImbalance, numberOfSlaves + 1);




	// FIXME remover codigo de validacao
	ImbalanceMatrix tempProcessClusterImbMatrix = bestSplitgraphClustering.getInterProcessImbalanceMatrix();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] *** zero-cost move VALIDATION *** Calculating external process imbalance matrix...";

	ImbalanceMatrix processClusterImbMatrixValidation = util.calculateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray,
						this->vertexImbalance, numberOfSlaves + 1);
	bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, processClusterImbMatrixValidation.pos, 0.001, 0.1);
	bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, processClusterImbMatrixValidation.neg, 0.001, 0.1);
	BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2;





	// Recalculates internal process imbalance array
	// este calculo full eh feito apenas uma vez (nao tao custoso assim)
	if(overloadedProcessList.size() > 0) {
		std::vector<Imbalance> processInternalImbalance = util.calculateProcessInternalImbalance(*g, currentSplitgraphClusterArray,
								globalClusterArray, numberOfProcesses);
		// TODO VALIDAR SE O NUMERO DO PRIMEIRO PROCESSO COM PARALLEL GRAPH COMECA EM 1 MESMO
		for(int px = 0; px < numberOfProcesses; px++) {
			bestClustering.setProcessImbalance(px, processInternalImbalance[px]);
		}
	}

	// Validacao do calculo da FO
	// FIXME
	clusteringgraph::validation::CCImbalanceCalculator calc = CCImbalanceCalculator::instance(g->getGraphFileLocation());
	Imbalance validation = calc.objectiveFunction(globalClusterArray);
	if(fabs(validation.getValue() - bestClustering.getImbalance().getValue()) > EPS) {
		BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] rebalanceClustersBetweenProcessesWithZeroCost: I(P) DOES NOT MATCH! Correct I(P) = "
			<< std::setprecision(5) << std::fixed << validation.getValue()
			<< " vs I(P) = " << bestClustering.getImbalance().getValue();
	} else {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] rebalanceClustersBetweenProcessesWithZeroCost: I(P) OK!";
	}

	// prints the new process x vertex configuration
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] RebalanceClusters done.";
	int np = bestSplitgraphClustering.getNumberOfClusters();
	for(int px = 0; px < np; px++) {
		int k = util.obtainListOfClustersFromProcess(*g, bestClustering, px).size();
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << px << ": num_vertices = "
				<< bestSplitgraphClustering.getClusterSize(px) << ", k = " << k;
	}
}

OutputMessage ImbalanceSubgraphParallelILS::runILSLocallyOnSubgraph(InputMessageParallelILS &imsgpils, clusteringgraph::SignedGraph* g) {
	// TODO duvida: deve chamar synchronize(*pgraph); antes?
	boost::mpi::communicator world;
	int myRank = world.rank();
	long N = g->getGlobalN();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] runILSLocallyOnSubgraph()";

	// only redistributes vertices if vertexList has information available
	// otherwise this method is being called from distributed ils initialization, when no movement occurs
	if(imsgpils.redistributeVertices) {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Invoking vertex redistribution between processes...";
		// moves the affected vertices between processes
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Global graph size is " << N;
		// ClusterArray splitgraphClusterArray(N, 0);
		// for(int i = 0; i < N; i++) {
		// 	splitgraphClusterArray[i] = 0;
		// }
		/*
		for(int i = 0; i < imsgpils.vertexList.size(); i++) {
			int v = imsgpils.vertexList[i];
			// mark the vertices that belong to this process (myRank)
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Assigning rank " << myRank << " to vertex " << v;
			assert(v < N);
			splitgraphClusterArray[v] = myRank;
		}*/
		// each worker process calls vertex redistribution locally
		redistributeVerticesInProcesses(g, imsgpils.splitgraphClusterArray);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Redistribution done.";
	} else {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Vertex redistribution disabled.";
	}

	// FIXME Obtains the updated globalVertexId array
	BOOST_LOG_TRIVIAL(info) << "Invoking Parallel Graph global vertex mapping...";
	// https://github.com/boostorg/graph_parallel/blob/develop/test/adjlist_redist_test.cpp
	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g->graph));

	BOOST_LOG_TRIVIAL(info) << "Obtaining global_index...";
	mpi_process_group pg = g->graph->process_group();
	boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
	  global_index(pg, num_vertices(*(g->graph)),
				   get(vertex_index, *(g->graph)), get(vertex_global, *(g->graph)));
	BOOST_LOG_TRIVIAL(info) << "Obtaining name_map...";
	std::stringstream ss;
	BGL_FORALL_VERTICES(v, *(g->graph), ParallelGraph) {
		int idx = get(name_map, v);
		ss << idx << "; ";
		assert(owner(v) == myRank);
	}
	BOOST_LOG_TRIVIAL(debug) << "Vertices in this process (P " << myRank << "): " << ss.str() << ".";

	// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
	// global vertex id array used in split graph
	std::vector<long> globalVertexId;
	BOOST_LOG_TRIVIAL(info) << "Obtaining individual global_index...";
	ClusterArray localSplitgraphClusterArray(N, -1);
	graph_traits<ParallelGraph>::vertex_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = vertices(*(g->graph)); vi != vi_end; ++vi) {  // For each vertex v
		// BOOST_LOG_TRIVIAL(info) << "Local vertex " << vi->local << " is global vertex " << get(name_map, *vi);
		// Vertex vx = get(name_map, *vi);
		int idx = get(name_map, *vi);//vx.id;
		globalVertexId.push_back(idx);
		localSplitgraphClusterArray[idx] = myRank;
	}

	// BOOST_LOG_TRIVIAL(info) << "Obtaining processor_map...";
	// property_map<ParallelGraph, vertex_rank_t>::type to_processor_map = get(vertex_rank, *(g->graph));

	// triggers the local ILS routine
	// Chooses between the sequential or parallel search algorithm
	NeighborhoodSearch* neigborhoodSearch;
	ClusteringProblemFactory problemFactory;
	NeighborhoodSearchFactory nsFactory(MPIUtil::ALL_MASTERS_FIRST, imsgpils.numberOfMasters, imsgpils.numberOfSearchSlaves);
	neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	GainFunctionFactory functionFactory(g);
	ConstructClustering defaultConstruct(functionFactory.build(imsgpils.gainFunctionType), randomSeed, imsgpils.alpha);
	ConstructClustering noConstruct(functionFactory.build(imsgpils.gainFunctionType), randomSeed, imsgpils.alpha, &imsgpils.CCclustering);
	ConstructClustering* construct = &defaultConstruct;
	if((imsgpils.problemType == problem::ClusteringProblem::RCC_PROBLEM) and (imsgpils.k < 0)) {
			construct = &noConstruct;
	}
	VariableNeighborhoodDescent vnd(*neigborhoodSearch, randomSeed, imsgpils.l, imsgpils.firstImprovementOnOneNeig,
			imsgpils.timeLimit);
	// Additional execution info
	ExecutionInfo info(imsgpils.executionId, imsgpils.fileId, imsgpils.outputFolder, myRank);
	resolution::ils::ILS resolution;
	resolution::ils::CUDAILS CUDAILS;
	clusteringgraph::Clustering bestClustering;
	// time spent on processing
	double timeSpent = 0.0;
	// info about the processed graph
	long n = num_vertices(*(g->graph)), e = num_edges(*(g->graph));
	boost::mpi::environment  env;
	boost::mpi::communicator comm;

	BOOST_LOG_TRIVIAL(info) << "[runILSLocallyOnSubgraph] Creating subgraph";
	if((n == 0) or (not imsgpils.runILS)) {
		if(n == 0) {
			BOOST_LOG_TRIVIAL(info) << "Empty subgraph, returning zero imbalance.";
		} else {
			BOOST_LOG_TRIVIAL(info) << "Local ILS execution disabled for this process, returning zero imbalance.";
		}
		std::vector<unsigned int> clusterProcessOrigin;
		std::vector<Imbalance> internalImbalance;
		ClusterArray cArray(g->getN(), clusteringgraph::Clustering::NO_CLUSTER);
		problem::CCProblem ccProblem;
		if(g->getN() > 0) {
			bestClustering = clusteringgraph::Clustering(cArray, *g, ccProblem, 0.0, 0.0, clusterProcessOrigin, internalImbalance);
		} else {
			ClusterArray cArray2(10, clusteringgraph::Clustering::NO_CLUSTER);
			bestClustering = clusteringgraph::Clustering(cArray2, *g, ccProblem, 0.0, 0.0, clusterProcessOrigin, internalImbalance);
		}
		// n = 0;
		e = 0;
		OutputMessage omsg(bestClustering, 0, 0.0, globalVertexId, 0, 0);
		omsg.splitgraphClusterArray = localSplitgraphClusterArray;
		return omsg;
	} else {
		// TODO TROCAR POR METODO EQUIVALENTE DA PARALLEL BGL
		// clusteringgraph::SignedGraph sg(*(g->graph), vertexList);
		// std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(*(g->graph));
		// for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
			// TODO TROCAR POR METODO EQUIVALENTE DA PARALLEL BGL
			// globalVertexId.push_back(sg.graph.local_to_global(*it));
		// }

		// codigo de teste da parallel bgl
		 // https://groups.google.com/forum/#!topic/boost-list/A_IOeEGWrWY
		 boost::mpi::environment  env;
		 boost::mpi::communicator comm;
		 /*
		 BGL_FORALL_VERTICES(v, *(g->graph), ParallelGraph)
		 {
			 BOOST_LOG_TRIVIAL(trace) << "V @ P" << comm.rank() << ": " << v.local;
		 }
		 int edgecount = 0;
		 BGL_FORALL_EDGES(e, *(g->graph), ParallelGraph)
		 {
			 BOOST_LOG_TRIVIAL(debug) << "E @ P" << comm.rank() << " v_src: " << boost::source(e,*(g->graph)).local
					<< " -> v_dst: " << boost::target(e, *(g->graph)).local << "; srccpu = " <<
				e.source_processor << " dstcpu = " << e.target_processor;
			 edgecount++;
		 }
		 BOOST_LOG_TRIVIAL(debug) << "Edgecount is " << edgecount; */

		// rebuilds construct clustering objects based on partial graph 'sg'
		GainFunctionFactory functionFactory(g);
		ConstructClustering defaultConstruct(functionFactory.build(imsgpils.gainFunctionType), randomSeed, imsgpils.alpha);
		ConstructClustering noConstruct(functionFactory.build(imsgpils.gainFunctionType), randomSeed, imsgpils.alpha, &imsgpils.CCclustering);
		ConstructClustering* construct = &defaultConstruct;
		if((imsgpils.problemType == problem::ClusteringProblem::RCC_PROBLEM) and (imsgpils.k < 0)) {
				construct = &noConstruct;
		}
		problem::ClusteringProblem& problem = problemFactory.build(imsgpils.problemType, imsgpils.k);
		// Each worker process processing his own partition
		n = num_vertices(*(g->graph));
		e = num_edges(*(g->graph));
		BOOST_LOG_TRIVIAL(info) << "Processing subgraph with n =  " << n << ", " << "e =  " << e;
		clusteringgraph::Clustering bestClustering;
		double timeSpent = 0.0;
		long num_comb = 0;
		// only executes CUDA ILS on vertex overloaded processes; reason: lack of GPU memory resources
		unsigned int numberOfProcesses = numberOfSlaves + 1;
		resolution::ils::ILS resolution;
		if(is_overloaded_process(vertexList) and imsgpils.cudaEnabled) {
			resolution::ils::CUDAILS cudails;
			try {
				bestClustering = cudails.executeILS(construct, &vnd, g, imsgpils.iter,
										imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
										problemFactory.build(imsgpils.problemType, imsgpils.k), info);
				num_comb = cudails.getNumberOfTestedCombinations();
				timeSpent = cudails.getTotalTimeSpent();
			} catch(std::exception& e) {  // possibly an exception caused by lack of GPU CUDA memory
				BOOST_LOG_TRIVIAL(error) << e.what() << "\n";
				BOOST_LOG_TRIVIAL(error) << "Lack of GPU memory detected: invoking sequential ILS (non CUDA)...";
				// try to run the standard sequential ILS algorithm => FALLBACK
				bestClustering = resolution.executeILS(construct, &vnd, g, imsgpils.iter,
														imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
														problemFactory.build(imsgpils.problemType, imsgpils.k), info);
				timeSpent = resolution.getTotalTimeSpent();
	            num_comb = resolution.getNumberOfTestedCombinations();
            }
		} else {
			bestClustering = resolution.executeILS(construct, &vnd, g, imsgpils.iter,
													imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
													problemFactory.build(imsgpils.problemType, imsgpils.k), info);
											timeSpent = resolution.getTotalTimeSpent();
			timeSpent = resolution.getTotalTimeSpent();
			num_comb = resolution.getNumberOfTestedCombinations();
		}

		assert(bestClustering.getClusterArray().size() == num_vertices(*(g->graph)));
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph P" << comm.rank() << ": num_edges = " <<
								num_edges(*(g->graph)) << " , num_vertices = " << num_vertices(*(g->graph)) << ", I(P) = "
								<< bestClustering.getImbalance().getValue() << ", k = " << bestClustering.getNumberOfClusters();
		OutputMessage omsg(bestClustering, num_comb, timeSpent, globalVertexId,
						num_vertices(*(g->graph)), num_edges(*(g->graph)));
		omsg.splitgraphClusterArray = localSplitgraphClusterArray;
		return omsg;
	}
}


void ImbalanceSubgraphParallelILS::moveClusterToDestinationProcessZeroCost(clusteringgraph::SignedGraph *g, clusteringgraph::Clustering& bestClustering,
		ProcessClustering& bestSplitgraphClustering, long clusterToMove, unsigned int sourceProcess, unsigned int destinationProcess) {
	// Simply moves the cluster to the destination process, with zero-cost imbalance change
	bestClustering.setProcessOrigin(clusterToMove, destinationProcess);
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	int np = bestSplitgraphClustering.getNumberOfClusters();
	long n = g->getGlobalN();
	boost::mpi::communicator world;

	// TODO INCLUIR AQUI CHAMADA AO REDISTRIBUTE DA BOOST PARALLEL BGL!
	// All worker processes must call the vertex redistribution!
	std::vector< std::vector<long> > verticesInProcess(np, std::vector<long>());
	for(long i = 0; i < n; i++) {
		long k = splitgraphClusterArray[i];
		assert(k < np);
		verticesInProcess[k].push_back(i);
		// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] vertex " << i << " goes to process " << k;
	}
	
	InputMessageParallelILS imsg;
	imsg.l = 1;
	imsg.isSplitGraph = true;
	imsg.isParallelGraph = true;
	imsg.runILS = false;
	for(int pr = 1; pr < np; pr++) {
		imsg.setVertexList(verticesInProcess[pr]);
		imsg.splitgraphClusterArray = splitgraphClusterArray;
		world.send(pr, MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Move message sent to process " << pr << " (runILS = false)";
	}

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Invoking redistribution on leader process (P0).";
	imsg.setVertexList(verticesInProcess[0]);
	imsg.splitgraphClusterArray = splitgraphClusterArray;
	OutputMessage leaderProcessingMessage = runILSLocallyOnSubgraph(imsg, g);

	// STEP 3: the leader receives the processing results
	// Maps the messages according to the movement executed
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Waiting for slaves return messages...";
	for(int pr = 1, i = 0; pr < np; pr++) {
		OutputMessage omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfSlaves) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Message received from process " << procNum << ".";
	}
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph Zero-cost move] Redistribution of vertices done.";

	// However, the internal and external imbalance values of the participating processes must be updated
	// Assumes the imbalance matrix between processes has already been updated
	// Updates the internal imbalance of the participating processes
	/*  DISABLED, SINCE rebalanceClustersBetweenProcessesWithZeroCost() INVOKES FULL RECALCULATION IN THE END OF THE PROCEDURE
	Imbalance imbSrc = util.calculateProcessInternalImbalance(g, bestClustering, sourceProcess);
	bestClustering.setProcessImbalance(sourceProcess, imbSrc);
	Imbalance imbDest = util.calculateProcessInternalImbalance(g, bestClustering, destinationProcess);
	bestClustering.setProcessImbalance(destinationProcess, imbDest);
	*/

	// validating imbalance value FIXME remover esta validacao
	// Calculates the internal imbalance sum (inside each process)
	/*
	std::vector<Imbalance> internalProcessImbalance = bestClustering.getInternalProcessImbalance();
	Imbalance internalImbalance = utilcalculateInternalImbalanceSumOfAllProcesses(internalProcessImbalance);
	// Calculates the external imbalance sum (between processes)
	Imbalance externalImbalance = util.calculateExternalImbalanceSumBetweenProcesses(processClusterImbMatrix);
	if(internalImbalance.getPositiveValue() + externalImbalance.getPositiveValue() +
			internalImbalance.getNegativeValue() + externalImbalance.getNegativeValue() != bestClustering.getImbalance().getValue()) {
		BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] REBALANCE Zero-cost imbalance does not match!";
	}
	BOOST_LOG_TRIVIAL(info) << "AFTER => internalImbalance vector contents: ";
	stringstream ss;
	int count = 0;
	for (std::vector<Imbalance>::const_iterator i = internalProcessImbalance.begin(); i != internalProcessImbalance.end(); ++i) {
		Imbalance imb = util.calculateProcessInternalImbalance(g, bestClustering, count++);
		ss << i->getValue() << " (" << imb.getValue() << "); ";
	}
	BOOST_LOG_TRIVIAL(info) << ss.str();
	*/
}

void ImbalanceSubgraphParallelILS::redistributeVerticesInProcesses(clusteringgraph::SignedGraph *g, const ClusterArray& newSplitgraphClustering) {
	boost::mpi::communicator world;
	int myRank = world.rank();
	// TODO redistribution test
	// https://github.com/boostorg/graph_parallel/blob/develop/test/adjlist_redist_test.cpp
	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g->graph));

	mpi_process_group pg = g->graph->process_group();
	boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
	  global_index(pg, num_vertices(*(g->graph)),
				   get(vertex_index, *(g->graph)), get(vertex_global, *(g->graph)));

	property_map<ParallelGraph, vertex_rank_t>::type to_processor_map = get(vertex_rank, *(g->graph));

	// Assign a new distribution
	graph_traits<ParallelGraph>::vertex_iterator vi, vi_end;
	std::vector<int> my_vertices_before;
	std::vector<int> moved_vertices;
	for (boost::tie(vi, vi_end) = vertices(*(g->graph)); vi != vi_end; ++vi) {
		// put(to_processor_map, vi, vi->local % num_processes(pg));  // floor(vi->local / (double)g->getN())
		// BOOST_LOG_TRIVIAL(debug) << "Obtaining newSplitgraphClustering for global vertex " << get(name_map, *vi);
		// BOOST_LOG_TRIVIAL(debug) << "Going to processor " << newSplitgraphClustering[get(name_map, *vi)];
		int idx = get(name_map, *vi);
		my_vertices_before.push_back(idx);
		put(to_processor_map, *vi, newSplitgraphClustering[idx]);
		if(newSplitgraphClustering[idx] != myRank) {
			moved_vertices.push_back(idx);
		}
	}

	// if (process_id(pg) == 0) {
	BOOST_LOG_TRIVIAL(info) << "[Distributed ILS] Redistributing vertices between processes...";
	// Perform the actual redistribution
	(g->graph)->redistribute(to_processor_map);

	std::vector<int> my_vertices_after;
    {
		typedef typename property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
		typedef typename property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
		typename property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g->graph));

		boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
		global_index(pg, num_vertices(*(g->graph)),
					get(vertex_index, *(g->graph)), get(vertex_global, *(g->graph)));
		typename graph_traits<ParallelGraph>::vertex_iterator vi, vi_end;
		for (boost::tie(vi, vi_end) = vertices(*(g->graph)); vi != vi_end; ++vi) {
			int idx = get(name_map, *vi);
			my_vertices_after.push_back(idx);
			// validation(idx) = process_id(pg);
		}
    }
    BOOST_LOG_TRIVIAL(debug) << "Moved vertex list: ";
    Print(moved_vertices);
	BOOST_LOG_TRIVIAL(debug) << "*** P" << process_id(pg) << ": Before redist, my vertices are: ";
	Print(my_vertices_before);// clusterArray;
	BOOST_LOG_TRIVIAL(debug) << "*** P" << process_id(pg) << ": After redist, they are: ";
	Print(my_vertices_after);// validation;

	BOOST_LOG_TRIVIAL(info) << "[Distributed ILS] Redistribution done.";
}

void ImbalanceSubgraphParallelILS::Print(const std::vector<int>& v){
	std::stringstream ss;
    for(unsigned i = 0; i< v.size(); ++i) {
        ss << v[i] << " ";
    }
    BOOST_LOG_TRIVIAL(debug) << ss.str();
}

void ImbalanceSubgraphParallelILS::Print(const std::vector<long>& v){
	std::stringstream ss;
    for(unsigned i = 0; i< v.size(); ++i) {
        ss << v[i] << " ";
    }
    BOOST_LOG_TRIVIAL(debug) << ss.str();
}

} /* namespace ils */
} /* namespace resolution */
