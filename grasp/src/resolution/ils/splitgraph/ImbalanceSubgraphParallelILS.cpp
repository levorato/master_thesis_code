/*
 * ImbalanceSubgraphParallelILS.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#include "./include/ImbalanceSubgraphParallelILS.h"
#include "util/include/ProcessUtil.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "../include/CUDAILS.h"
#include "problem/include/RCCProblem.h"
#include "util/include/RandomUtil.h"

#include "../../construction/include/GainFunctionFactory.h"
#include "../../construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"

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

namespace ublas = boost::numeric::ublas::detail;

namespace resolution {
namespace ils {

using namespace resolution::construction;
using namespace boost::algorithm;
using namespace std;
using namespace util;
using namespace util::parallel;

ImbalanceSubgraphParallelILS::ImbalanceSubgraphParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
		const bool& split, const bool& cuda) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split), cudaEnabled(cuda), vertexImbalance(), timeSpentAtIteration() {

}

ImbalanceSubgraphParallelILS::~ImbalanceSubgraphParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering ImbalanceSubgraphParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	stringstream constructivePhaseResults;
	stringstream iterationResults;
	stringstream iterationTimeSpent;
	unsigned int numberOfProcesses = numberOfSlaves + 1;

	// *** STEP A => PREPROCESSING PHASE OF DISTRIBUTED METAHEURISTIC: New split graph partitioning here
	// Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// TODO PARAMETRIZAR DIVISAR POR VERTICE OU POR SOMA DE ARESTAS
	Clustering splitgraphClustering = preProcessSplitgraphPartitioning(g, problem, true);
	BOOST_LOG_TRIVIAL(info) << "Initial split graph processor partitioning I(P) = " << splitgraphClustering.getImbalance().getValue();

	// *** STEP B => Calls the individual ILS processing for each subgraph, invoking Parallel ILS with MPI
	// WARNING: THE SPLITGRAPHCLUSTERARRAY IS A SEPARATE DATA STRUCTURE, DIFFERENT THAN THE CURRENT CLUSTERING
	// IT CONTROLS THE PARTITIONING BETWEEN PROCESSORS (SPLIT GRAPH)
	ClusterArray initialSplitgraphClusterArray = splitgraphClustering.getClusterArray();
	ImbalanceMatrix processClusterImbMatrix = calculateProcessToProcessImbalanceMatrix(*g, initialSplitgraphClusterArray);
	Clustering Cc = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info,
			initialSplitgraphClusterArray, processClusterImbMatrix, splitgraphClustering);

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
	Cc.printClustering(g->getN());
	Clustering bestClustering(Cc);
	Clustering bestSplitgraphClustering(splitgraphClustering);

	// STEP 2: TRY TO IMPROVE THE GLOBAL SOLUTION THROUGH A DISTRIBUTED METAHEURISTIC WHICH
	// EXCHANGES VERTICES BETWEEN PROCESSES
	int improvementCount = 0;
	do {
		// 1. O processo mestre será o responsável por manter duas estruturas de dados de controle do imbalance:
		//		Um vetor com os vértices que mais contribuem para I(P) (soma de imbalance de cada vertice)
		//		Uma matriz com a soma do imbalance entre processos

		// VND: for now, using 2 neighborhood structures, first 1-opt cluster move, then 1-opt vertex move
		improvementCount = variableNeighborhoodDescent(g, bestSplitgraphClustering,
				bestClustering, numberOfProcesses,
				processClusterImbMatrix, construct, vnd,
				iter, iterMaxILS, perturbationLevelMax,
				problem, info, timeSpentInILS);

		// Last step. Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();

		// Write the results into ostream os, using csv format
		// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		// iteration results are only being shown when the solution improved
		iterationResults << (improvementCount+1) << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
				<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
				<< "," << fixed << setprecision(4) << timeSpentInILS << "\n";
		timer.resume();
		start_time = timer.elapsed();
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInILS >= vnd->getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}
	} while(improvementCount > 0);

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

	// 7. Stops the timer and stores the elapsed time
	timer.stop();
	end_time = timer.elapsed();

	// 8. Write the results into ostream os, using csv format
	// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
	timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
	iterationResults << (improvementCount+2) << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
			<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
			<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) << "\n";
	Imbalance bestValue = bestClustering.getImbalance();
	int iterationValue = 1;
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
				<< "," << bestValue.getPositiveValue()
				<< "," << bestValue.getNegativeValue()
				<< setprecision(0)
				<< "," << bestClustering.getNumberOfClusters()
				<< "," << (iterationValue+1)
				<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) // timeSpentOnBestSolution
				<< "," << iter
				<< "," << numberOfTestedCombinations << "\n";

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after 1 ILS parallel graph iteration: I(P) = " << bestClustering.getImbalance().getValue();

	// Saves the iteration results to csv file
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Save the time spent at each iteration by each process to csv file
	generateOutputFile(problem, iterationTimeSpent, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timespent-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);

	return bestClustering;
}

Clustering ImbalanceSubgraphParallelILS::preProcessSplitgraphPartitioning(SignedGraph *g, ClusteringProblem& problem, bool partitionByVertex) {
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

	// For each process pi
	for(int pi = 0; pi < numberOfProcesses; pi++) {
		BOOST_LOG_TRIVIAL(debug) << "Processing partition for processor " << pi;
		// *** STEP 1: Chooses a node with the smallest negative cardinality INSIDE the residual graph Gr
		//  (will probably result in the smallest increase in the imbalance)
		Clustering GrCluster(*g);  // uses Clustering data structure to control the elements in the set
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
		Clustering Si(*g);  // Si = empty, uses Clustering data structure to control the elements in the set
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

			boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
			UndirectedGraph::edge_descriptor e;
			UndirectedGraph::out_edge_iterator f2, l2;
			std::vector<double> deltaPosDegree(n, (double)0);
			// keeps an array containing every positive edge between ni_max and every other vertex -> O(n)
			for (boost::tie(f2, l2) = out_edges(max_ni, g->graph); f2 != l2; ++f2) {
				e = *f2;
				long j = target(*f2, g->graph);
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
	Clustering Cc(splitgraphClusterArray, *g, problem);

	return Cc;
}

Clustering ImbalanceSubgraphParallelILS::runDistributedILS(ConstructClustering *construct,
		VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
		const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
		ClusterArray& splitgraphClusterArray, ImbalanceMatrix& processClusterImbMatrix, Clustering& currentSolution) {

	Clustering c = distributeSubgraphsBetweenProcessesAndRunILS(construct, vnd, g,
			iter, iterMaxILS, perturbationLevelMax, problem, info, splitgraphClusterArray, processClusterImbMatrix);
	int retries = 1;

	// TODO return the last global solution found by distributed ILS or the less worse one?
	while((c.getImbalance().getValue() >= currentSolution.getImbalance().getValue()) and (retries <= MAX_ILS_RETRIES)) {
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Found worse global solution with I(P) = "
				<< c.getImbalance().getValue() << " (best known I(P) = " << currentSolution.getImbalance().getValue()
				<< "). Retry " << retries << "...";
		retries++;
		c = distributeSubgraphsBetweenProcessesAndRunILS(construct, vnd, g,
				iter, iterMaxILS, perturbationLevelMax, problem, info, splitgraphClusterArray, processClusterImbMatrix);
	}
	return c;
}

Clustering ImbalanceSubgraphParallelILS::distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
		VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
		const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
		ClusterArray& splitgraphClusterArray, ImbalanceMatrix& processClusterImbMatrix) {

	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	long n = g->getN();
	unsigned int numberOfProcesses = numberOfSlaves + 1;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}

	std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
	for(long i = 0; i < n; i++) {
		long k = splitgraphClusterArray[i];
		verticesInCluster[k].push_back(i);
	}

	// Creates numberOfProcesses subgraphs
	/*
	std::vector<SubGraph> subgraphList;
	// each subgraph will have a subset of the main graph's nodes and edges, based on the previous clustering
	for(int p = 0; p < numberOfProcesses; p++) {
		SubGraph sg = (g->graph).create_subgraph(); //verticesInCluster[k].begin(), verticesInCluster[k].end());
		subgraphList.push_back(sg);  // --> SUBGRAPH COPY CTOR NOT WORKING!!!
		for(std::vector<long>::iterator it = verticesInCluster[p].begin(); it != verticesInCluster[p].end(); it++) {
			add_vertex(*it, subgraphList.back());
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << p;
		}

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_edges = " <<
				num_edges(subgraphList.back()) << " , num_vertices = " << num_vertices(subgraphList.back());
	} */

	for(int i = 0; i < numberOfSlaves; i++) {
		// 2. Distribute numberOfSlaves graph parts between the ILS Slave processes
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
							problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd->getTimeLimit(),
							numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k, true);
		if(k < 0) {
			imsg.setClustering(&CCclustering);
		}
		imsg.setVertexList(verticesInCluster[i+1]);
		world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Message sent to process " << slaveList[i];
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Size of ILS Input Message: " << (sizeof(imsg)/1024.0) << "kB.";
	}
	// 2.1. the leader does its part of the work: runs ILS using the first part of the divided graph
	BOOST_LOG_TRIVIAL(info) << "[Master process] Creating subgraph";
	Clustering leaderClustering;
	// Structure containing processing time of each MH slave
	std::vector<double> timeSpent(numberOfProcesses, double(0.0));
	// structure containing the vertex id in the global graph / clustering
	std::vector<long> globalVertexId;
	if(verticesInCluster[0].size() == 0) {
		BOOST_LOG_TRIVIAL(info) << "Empty subgraph, returning zero imbalance.";
		std::vector<unsigned int> clusterProcessOrigin;
		ClusterArray cArray(g->getN(), Clustering::NO_CLUSTER);
		leaderClustering = Clustering(cArray, *g, problem, 0.0, 0.0, clusterProcessOrigin);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph 0" << ": num_edges = 0" <<
								" , num_vertices = 0" << ", I(P) = "
								<< leaderClustering.getImbalance().getValue() << ", k = " << leaderClustering.getNumberOfClusters();
	} else {
		CUDAILS cudails;
		SignedGraph sg(g->graph, verticesInCluster[0]);
		// sg.graph = (g->graph).create_subgraph(verticesInCluster[0].begin(), verticesInCluster[0].end());

		/*
		for(std::vector<long>::iterator it = verticesInCluster[0].begin(); it != verticesInCluster[0].end(); it++) {
			add_vertex(*it, sg.graph);
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << k;
		} */
		BOOST_LOG_TRIVIAL(info) << "Processing subgraph with n =  " << num_vertices(sg.graph) << ", " << "e =  " << num_edges(sg.graph);

		// rebuilds construct clustering objects based on partial graph 'sg'
		GainFunctionFactory functionFactory(&sg);
		ConstructClustering defaultConstruct(functionFactory.build(construct->getGainFunctionType()),
				construct->getRandomSeed(), construct->getAlpha());
		ConstructClustering noConstruct(functionFactory.build(construct->getGainFunctionType()),
				construct->getRandomSeed(), construct->getAlpha(), construct->getCCclustering());
		ConstructClustering* construct2 = &defaultConstruct;
		if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			RCCProblem& rp = static_cast<RCCProblem&>(problem);
			int RCCk = rp.getK();
			if(RCCk < 0) {
				construct2 = &noConstruct;
			}
		}
		// Leader processing his own partition
		leaderClustering = cudails.executeILS(construct2, vnd, &sg, iter, iterMaxILS, perturbationLevelMax, problem, info);
		timeSpent[0] = cudails.getTotalTimeSpent();
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph 0" << ": num_edges = " <<
								num_edges(sg.graph) << " , num_vertices = " << num_vertices(sg.graph) << ", I(P) = "
								<< leaderClustering.getImbalance().getValue() << ", k = " << leaderClustering.getNumberOfClusters();
	}

	// 3. the leader receives the processing results
	OutputMessage omsg;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Waiting for slaves return messages.";
	// Global cluster array
	ClusterArray globalClusterArray(g->getN(), 0);
	long clusterOffset = 0;
	double internalImbalancePosSum = 0.0;
	double internalImbalanceNegSum = 0.0;
	std::vector<unsigned int> clusterProcessOrigin;
	for(int i = 0; i < numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
		int procNum = stat.source();
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ". Obj = " <<
				omsg.clustering.getImbalance().getValue();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// stores the time spent by each process
		timeSpent[i+1] = omsg.timeSpent;
		// 4. Merge the partial solutions into a global solution for the whole graph
		ClusterArray localClusterArray = omsg.clustering.getClusterArray();
		internalImbalancePosSum += omsg.clustering.getImbalance().getPositiveValue();
		internalImbalanceNegSum += omsg.clustering.getImbalance().getNegativeValue();

		long nc = omsg.clustering.getNumberOfClusters();
		if(nc > 0) {
			assert(localClusterArray.size() == omsg.globalVertexId.size());
			for(long v = 0; v < localClusterArray.size(); v++) {
				// Obtains vertex v's number in the global graph
				long vglobal = omsg.globalVertexId[v];
				globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
				/*
				BOOST_LOG_TRIVIAL(info) << "Vertex " << v << " in local solution becomes vertex " << vglobal <<
						" in global solution and belongs to cluster " << localClusterArray[v] <<
						", that is, global cluster " << globalClusterArray[vglobal]; */
			}
			clusterOffset += nc;
			// all the clusters in this interval belong to process 'stat.source()'
			for(int clusterCount = 0; clusterCount < nc; clusterCount++) {
				clusterProcessOrigin.push_back(procNum);
			}
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << procNum << ": num_edges = " <<
						omsg.num_edges << " , num_vertices = " << omsg.num_vertices << ", I(P) = " << omsg.clustering.getImbalance().getValue()
						 << ", k = " << nc;
	}
	// includes the leader's processing result as well
	// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
	if(verticesInCluster[0].size() > 0) {
		SignedGraph sg(g->graph, verticesInCluster[0]);
		std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(sg.graph);
		for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
			globalVertexId.push_back(sg.graph.local_to_global(*it));
		}
	}
	// 4. Merges the leader's partial solution into a global solution for the whole graph
	if(leaderClustering.getNumberOfClusters() > 0) {
		ClusterArray localClusterArray = leaderClustering.getClusterArray();
		for(long v = 0; v < localClusterArray.size(); v++) {
			// Obtains vertex v's number in the global graph
			long vglobal = globalVertexId[v];
			globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
			/*
			BOOST_LOG_TRIVIAL(info) << "Vertex " << v << " in local solution becomes vertex " << vglobal <<
									" in global solution and belongs to cluster " << localClusterArray[v] <<
									", that is, global cluster " << globalClusterArray[vglobal]; */
		}
	}
	long nc = leaderClustering.getNumberOfClusters();
	clusterOffset += nc;
	// all clusters in this interval belong to process zero
	for(int clusterCount = 0; clusterCount < nc; clusterCount++) {
		clusterProcessOrigin.push_back(info.processRank);
	}
	internalImbalancePosSum += leaderClustering.getImbalance().getPositiveValue();
	internalImbalanceNegSum += leaderClustering.getImbalance().getNegativeValue();

	// Calculates the external imbalance sum (between processes)
	double externalImbalancePosSum = 0.0;
	double externalImbalanceNegSum = 0.0;
	for (unsigned i = 0; i < processClusterImbMatrix.pos.size1(); ++i) {
		for (unsigned j = 0; j < processClusterImbMatrix.pos.size2(); ++j) {
			if(i != j) {
				externalImbalancePosSum += processClusterImbMatrix.pos(i, j);
				externalImbalanceNegSum += processClusterImbMatrix.neg(i, j);
			}
		}
	}

	// 5. Builds the clustering with the merge of each process local ILS result
	// EVITA O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
	Clustering globalClustering(globalClusterArray, *g, problem, internalImbalancePosSum + externalImbalancePosSum,
			internalImbalanceNegSum + externalImbalanceNegSum, clusterProcessOrigin);

	// BOOST_LOG_TRIVIAL(info) << "Separate internal imbalance I(P) = " << (internalImbalanceSum) << " and external I(P) = " << (externalImbalanceSum);
	//BOOST_LOG_TRIVIAL(info) << "Separate clustering imbalance I(P) = " << (internalImbalancePosSum + externalImbalancePosSum +
	//		internalImbalanceNegSum + externalImbalanceNegSum)
	//		<< " versus full calculated I(P) = " << globalClustering.getImbalance().getValue();

	// 6. Stores the time spent in this iteration of distributed MH
	timeSpentAtIteration.push_back(timeSpent);

	return globalClustering;
}

ImbalanceMatrix ImbalanceSubgraphParallelILS::calculateProcessToProcessImbalanceMatrix(
		SignedGraph& g, ClusterArray& myCluster) {

	long nc = numberOfSlaves + 1;
	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	// Process to process matrix containing positive / negative contribution to imbalance
	ImbalanceMatrix clusterImbMatrix(nc);
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	// Vector containing each vertex contribution to imbalance: vertexImbalance

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			bool sameCluster = (myCluster[targ.id] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
				clusterImbMatrix.neg(ki, myCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
				clusterImbMatrix.pos(ki, myCluster[targ.id]) += fabs(weight);
			}
		}
		vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
	return clusterImbMatrix;
}

void ImbalanceSubgraphParallelILS::updateProcessToProcessImbalanceMatrix(SignedGraph& g,
			const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix) {

	long nc = numberOfSlaves + 1;
	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i in listOfModifiedVertices
	for(long item = 0; item < listOfModifiedVertices.size(); item++) {
		long i = listOfModifiedVertices[item];
		long old_ki = previousSplitgraphClusterArray[i];
		long new_ki = newSplitgraphClusterArray[i];
		UndirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			// Etapa 1: subtracao dos imbalances antigos
			long old_kj = previousSplitgraphClusterArray[targ.id];
			bool sameCluster = (old_kj == old_ki);
			// TODO VERIFICAR SE DEVE SER RECALCULADO TAMBEM O PAR OPOSTO DA MATRIZ: (KJ, KI)
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j were in the same cluster
				processClusterImbMatrix.neg(old_ki, old_kj) -= fabs(weight);
				processClusterImbMatrix.neg(old_kj, old_ki) -= fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j were NOT in the same cluster
				processClusterImbMatrix.pos(old_ki, old_kj) -= fabs(weight);
				processClusterImbMatrix.pos(old_kj, old_ki) -= fabs(weight);
			}
			// Etapa 2: acrescimo dos novos imbalances
			long new_kj = newSplitgraphClusterArray[targ.id];
			sameCluster = (new_kj == new_ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are now in the same cluster
				processClusterImbMatrix.neg(new_ki, new_kj) += fabs(weight);
				processClusterImbMatrix.neg(new_kj, new_ki) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				processClusterImbMatrix.pos(new_ki, new_kj) += fabs(weight);
				processClusterImbMatrix.pos(new_kj, new_ki) += fabs(weight);
			}
		}
		// NAO EH NECESSARIO ATUALIZAR O VERTEX IMBALANCE ABAIXO, POIS A BUSCA LOCAL (VND) EH FIRST IMPROVEMENT
		// vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
}

Coordinate ImbalanceSubgraphParallelILS::findMaximumElementInMatrix(ImbalanceMatrix &mat) {
	int x = 0, y = 0;
	double maxValue = -DBL_MAX;

	// TODO INCLUDE HANDLING FOR MATH PRECISION DOUBLE (EPSILON)
	for(int i = 0; i < mat.pos.size1(); i++) {
		for(int j = 0; j < mat.pos.size2(); j++) {
			if((mat.pos(i, j) + mat.neg(i, j)) > maxValue) {
				maxValue = mat.pos(i, j) + mat.neg(i, j);
				x = i;
				y = j;
			}
		}
	}
	return Coordinate(x, y, maxValue);
}

std::vector< Coordinate > ImbalanceSubgraphParallelILS::getMatrixElementsAsList(ImbalanceMatrix &mat) {
	std::vector< Coordinate > returnList;

	for(int i = 0; i < mat.pos.size1(); i++) {
		for(int j = 0; j < mat.pos.size2(); j++) {
			returnList.push_back(Coordinate(i, j, (mat.pos(i, j) + mat.neg(i, j))));
		}
	}
	return returnList;
}

long ImbalanceSubgraphParallelILS::findMostImbalancedVertexInProcessPair(SignedGraph& g,
		ClusterArray& splitGraphCluster, ClusterArray& globalCluster, Coordinate processPair) {

	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	double maxImbalance = -DBL_MAX;
	long maxVertex = 0;

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = globalCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// only processes vertexes belonging to the process pair
		if((splitGraphCluster[i] == processPair.x) or (splitGraphCluster[i] == processPair.y)) {
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
				e = *f;
				Vertex src = source(e, g.graph), targ = target(e, g.graph);
				double weight = ew[e].weight;
				long kj = globalCluster[targ.id];
				bool sameCluster = (kj == ki);
				if(weight < 0 and sameCluster) {  // negative edge
					// i and j are in the same cluster
					negativeSum += fabs(weight);
				} else if(weight > 0 and (not sameCluster)){ // and ((kj == processPair.x) or (kj == processPair.y))) {  // positive edge
					// only counts the imbalance between the clusters / processes in the processPair (x and y)
					// i and j are NOT in the same cluster
					positiveSum += weight;
				}
			}
			// checks if the imbalance caused by vertex i is the biggest one
			// TODO INCLUDE HANDLING FOR MATH PRECISION DOUBLE (EPSILON)
			if(positiveSum + negativeSum > maxImbalance) {
				maxImbalance = positiveSum + negativeSum;
				maxVertex = i;
			}
		}
	}
	// BOOST_LOG_TRIVIAL(info) << "Vertex " << maxVertex << " has imbalance sum of " << maxImbalance;
	return maxVertex;
}

std::vector<Coordinate> ImbalanceSubgraphParallelILS::obtainListOfImbalancedClusters(SignedGraph& g,
		ClusterArray& splitGraphCluster, Clustering& globalClustering) {

	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	// Cluster to cluster matrix containing positive / negative contribution to imbalance
	long nc = globalClustering.getNumberOfClusters();
	ClusterArray globalCluster = globalClustering.getClusterArray();
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);
	matrix<double> clusterEdgeSumMatrix = zero_matrix<double>(nc, nc);
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = globalCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			bool sameCluster = (globalCluster[targ.id] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				clusterImbMatrix(ki, globalCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				clusterImbMatrix(ki, globalCluster[targ.id]) += fabs(weight);
			}
			clusterEdgeSumMatrix(ki, globalCluster[targ.id]) += fabs(weight);
		}
	}

	std::vector<Coordinate> imbalancedClusterList;
	// For each cluster ki
	for(long ki = 0; ki < nc; ki++) {
		double percentageOfInternalNegativeEdges = 0.0;
		double percentageOfExternalPositiveEdges = 0.0;

		double sumOfInternalEdges = clusterEdgeSumMatrix(ki, ki);
		double sumOfInternalNegativeEdges = clusterImbMatrix(ki, ki);
		double sumOfExternalEdges = 0.0;
		double sumOfPositiveExternalEdges = 0.0;
		for(long kj = 0; kj < nc; kj++) {
			if(ki != kj) {
				sumOfExternalEdges += clusterEdgeSumMatrix(ki, kj) + clusterEdgeSumMatrix(kj, ki);
				sumOfPositiveExternalEdges += clusterImbMatrix(ki, kj) + clusterImbMatrix(kj, ki);
			}
		}
		if(sumOfInternalEdges > 0) {
			percentageOfInternalNegativeEdges = (100 * sumOfInternalNegativeEdges) / sumOfInternalEdges;  // item A
		}
		if(sumOfExternalEdges > 0) {
			percentageOfExternalPositiveEdges = (100 * sumOfPositiveExternalEdges) / sumOfExternalEdges;  // item B
		}
		double imbalancePercentage = max(percentageOfInternalNegativeEdges, percentageOfExternalPositiveEdges);
		double imbalance = sumOfInternalNegativeEdges + sumOfPositiveExternalEdges;
		//double imbalancePercentage = percentageOfExternalPositiveEdges;
		//imbalancedClusterList.push_back(Coordinate(ki, globalClustering.getClusterSize(ki), imbalancePercentage));
		// WHEN MOVING CLUSTERS BETWEEN PROCESSES, ORDERING BY CLUSTER IMBALANCE,
		// THE USE OF ABSOLUTE IMBALANCE VALUES (BELOW) YIELDS BETTER RESULTS THAN PERCENTUAL IMBALANCE (ABOVE)
		imbalancedClusterList.push_back(Coordinate(ki, globalClustering.getClusterSize(ki), imbalance));
	}
	return imbalancedClusterList;
}

std::vector<Coordinate> ImbalanceSubgraphParallelILS::obtainListOfOverloadedProcesses(SignedGraph& g,
		Clustering& splitGraphClustering) {

	long n = g.getN();
	int numberOfProcesses = splitGraphClustering.getNumberOfClusters();
	// A vertex-overloaded process is a process with more than (n / numberOfProcesses) vertices.
	long numberOfEquallyDividedVertices = (long)ceil(n / (double)numberOfProcesses);
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] List of overloaded processes (with more than "
			<< numberOfEquallyDividedVertices << " vertices): ";

	std::vector<Coordinate> overloadedProcessList;
	stringstream processListStr;
	// For each process px
	for(long px = 0; px < numberOfProcesses; px++) {
		if(splitGraphClustering.getClusterSize(px) > numberOfEquallyDividedVertices) {
			overloadedProcessList.push_back(Coordinate(px, 0, splitGraphClustering.getClusterSize(px)));
			processListStr << px << " ";
		}
	}
	if(processListStr.str() == "") {
		BOOST_LOG_TRIVIAL(debug) << "None.";
	} else {
		BOOST_LOG_TRIVIAL(debug) << processListStr.str() << " ";
	}
	return overloadedProcessList;
}

std::vector<Coordinate> ImbalanceSubgraphParallelILS::obtainListOfClustersFromProcess(SignedGraph& g,
		Clustering& globalClustering, int processNumber) {

	long nc = globalClustering.getNumberOfClusters();
	ClusterArray globalCluster = globalClustering.getClusterArray();
	// to which process a global cluster belongs to?
	const std::vector<unsigned int> clusterProcessOrigin = globalClustering.getClusterProcessOrigin();

	std::vector<Coordinate> clusterList;
	for(long clusterNum = 0; clusterNum < nc; clusterNum++) {
		if(clusterProcessOrigin[clusterNum] == processNumber) {
			clusterList.push_back(Coordinate(clusterNum, processNumber, globalClustering.getClusterSize(clusterNum)));
		}
	}

	return clusterList;
}

// DISABLED neighborhood!
bool ImbalanceSubgraphParallelILS::moveVertex1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	// => Calculates the imbalance info between process partitions and of each vertex individually
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	ClusterArray bestSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	// 2. O processo mestre elege dois processos condidatos a participar de uma movimentação de vertices
	// a) O processo mestre escolhe os dois processos com maior valor de imbalance entre si (usando a matriz 1b)
	// Chooses the top 20% biggest elements in the matrix
	// Coordinate processPair = findMaximumElementInMatrix(processClusterImbMatrix);
	std::vector< Coordinate > imbMatrixElementList = getMatrixElementsAsList(processClusterImbMatrix);
	int numberOfProcessCombinations = (int)floor((numberOfProcesses * (numberOfProcesses - 1)) / double(2.0));  // symmetric matrix
	int desiredNumberOfProcessPairs = (int)ceil(numberOfProcessCombinations * PERCENTAGE_OF_PROCESS_PAIRS);
	bool foundBetterSolution = false;

	// obtains and removes the next biggest element from the list
	for(int processPairCount = 0; processPairCount < desiredNumberOfProcessPairs; processPairCount++) {
		std::vector< Coordinate >::iterator pos = std::max_element(imbMatrixElementList.begin(),
				imbMatrixElementList.end(), coordinate_ordering_asc());
		Coordinate processPair = *pos;
		imbMatrixElementList.erase(pos);

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] The next process pair with biggest imbalance is (" << processPair.x
				<< ", " << processPair.y << ") with imbalance = " << processPair.value;

		// O processo mestre escolhe um vértice a ser movimentado (move) de um processo para o outro (usando o vetor 1 A)
		// select a vertex v to move - v is the vertex that causes the biggest imbalance
		// TODO the biggest imbalance in the whole graph or only between the pair of processes x and y?
		// for now it is the vertex v that belongs to the processPair and causes the biggest imbalance in the whole graph
		// Tries to move the top 25% of most imbalanced vertices
		long numberOfVerticesInProcessPair = bestSplitgraphClustering.getClusterSize(processPair.x)
				+ bestSplitgraphClustering.getClusterSize(processPair.y);
		long numberOfImbalancedVerticesToBeMoved = (int)ceil(numberOfVerticesInProcessPair * PERCENTAGE_OF_MOST_IMBALANCED_VERTICES_TO_BE_MOVED);
		ClusterArray splitgraphClusterList = bestSplitgraphClustering.getClusterArray();
		bool foundBetterSolution = false;
		for(long niv = 0; niv < numberOfImbalancedVerticesToBeMoved; niv++) {
			long v = findMostImbalancedVertexInProcessPair(*g, splitgraphClusterList, globalClusterArray, processPair);
			int currentProcess = splitgraphClusterList[v];
			// invalidates vertex v from the cluster array (== removes from the list), so that it is not chosen again
			splitgraphClusterList[v] = Clustering::NO_CLUSTER;
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Choosing vertex " << v << " as the next most imbalanced one.";

			// TODO PENSAR EM UMA HEURISTICA PARA ESCOLHER O PROCESSO DE DESTINO DO VERTICE V COM MAIOR IMBALANCE
			// OU ENTAO TENTAR MOVER O VERTICE V PARA QUALQUER OUTRO PROCESSO (ITERAR)
			ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
			ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
			for(int px = 0; px < numberOfProcesses; px++) {
				if(px != currentProcess) {
					BOOST_LOG_TRIVIAL(info) << "Trying to move vertex " << v << " to process " << px;
					// Realiza a movimentacao de um vertice (vertex 1-opt move)
					tempSplitgraphClusterArray[v] = px;
					std::vector<long> listOfModifiedVertices;
					listOfModifiedVertices.push_back(v);
					// O RECALCULO DA MATRIZ ABAIXO EH FEITO DE FORMA INCREMENTAL, reutilizando a matrix processClusterImbMatrix
					//ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
					ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;
					// Realiza o calculo do delta da matriz de imbalance entre processos
					updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
							listOfModifiedVertices, tempProcessClusterImbMatrix);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";
					// Valida se a matriz incremental e a full sao iguais
					//bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
					//bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
					//BOOST_LOG_TRIVIAL(info) << "*** As matrizes de delta sao iguais: " << igual << " e " << igual2;

					// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
					Clustering newGlobalClustering = runDistributedILS(construct, vnd, g,
								iter, iterMaxILS, perturbationLevelMax, problem, info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix,
								bestClustering);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
					if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
						bestClustering = newGlobalClustering;

						// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
						processClusterImbMatrix = tempProcessClusterImbMatrix;

						// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
						bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
						foundBetterSolution = true;
						break;
					} else {
						numberOfFrustratedSolutions++;
					}
				}
			}
			if(foundBetterSolution) {
				break;
			}
			// d) Dúvida: rodar um novo ILS inteiro para cada movimentação ou apenas uma busca local?
			// e) Cada processo envia o resultado (novo clustering) para o mestre e este calcula se vale a pena manter a alteração
			// f) O cálculo feito pelo mestre baseia-se na variação do imbalance dentro de cada processo (delta I(P) int) mais o novo imbalance externo, entre os processos (I(P) ext)
			// 		(a variação dentro de cada processo já é fornecida por cada processo, via mensagem – diferença do I(P) anterior)
			// 		(a variação do imbalance entre os processos é calculada varrendo a vizinhança do vértice movido)
			// g) Se gerar um imbalance pior, o movimento pode ser rejeitado, gerando uma volta da configuração anterior, por isso cada processo deve guardar o clustering anterior
		}
		if(foundBetterSolution)  break;
	}
	return foundBetterSolution;
}

bool ImbalanceSubgraphParallelILS::moveCluster1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   Move  cluster  1 - opt  ***********************************
	// b) The master process chooses a cluster to be moved to another process
	// select a cluster c to move - c is the cluster that causes the biggest imbalance
	// and, at the same time, giving priority to moving a cluster from an overloaded process to a less-loaded one
	// (priotitize vertex load balancing between processes)
	// TODO the biggest imbalance in the whole graph or only between the pair of processes x and y?
	// for now it is the vertex c that belongs to the processPair and causes the biggest imbalance in the whole graph
	// Tries to move the top 25% of most imbalanced clusters
	// TODO fazer os calculos corretos aqui, para clusters
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	std::vector<Coordinate> clusterImbalanceList = obtainListOfImbalancedClusters(*g, splitgraphClusterArray, bestClustering);
	// sorts the list of imbalanced clusters according to the percentage of imbalance (value)
	std::sort(clusterImbalanceList.begin(), clusterImbalanceList.end(), coordinate_ordering_desc());
	long numberOfImbalancedClustersToBeMoved = (int)ceil(clusterImbalanceList.size() * PERCENTAGE_OF_MOST_IMBALANCED_CLUSTERS_TO_BE_MOVED);
	if(numberOfImbalancedClustersToBeMoved == 0) {
		numberOfImbalancedClustersToBeMoved = 1;
	}
	bool foundBetterSolution = false;
	long n = g->getN();
	long nc = bestClustering.getNumberOfClusters();
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** moveCluster1opt";
	for(long nic = 0; nic < numberOfImbalancedClustersToBeMoved; nic++) {
		long clusterToMove = clusterImbalanceList[nic].x;
		std::vector<long> listOfModifiedVertices = getListOfVeticesInCluster(*g, bestClustering, clusterToMove);
		// finds out to which process the cluster belongs to
		int currentProcess = splitgraphClusterArray[listOfModifiedVertices[0]];
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Move " << nic << ": Moving global cluster " << clusterToMove << " of size " <<
				bestClustering.getClusterSize(clusterToMove) << " and imbalance = " << clusterImbalanceList[nc - nic - 1].value << ".";

		// LOAD BALANCING: Try to move the cluster to a process with less vertices than a given threshold (less loaded process)
		ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
		// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
		// TODO possible parametrization here!
		long maxVerticesAllowedInProcess = (2 * n) / bestSplitgraphClustering.getNumberOfClusters();
		for(int px = 0; px < numberOfProcesses; px++) {
			if((px != currentProcess)) {  // TODO the line below was disabled for narrowing the search space too much...
					// (bestSplitgraphClustering.getClusterSize(px) + listOfModifiedVertices.size() < maxVerticesAllowedInProcess)) {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to move global cluster " << clusterToMove << " to process " << px;
				ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
				// Realiza a movimentacao dos vertices de um cluster especifico (cluster move)
				for(long elem = 0; elem < listOfModifiedVertices.size(); elem++) {
					tempSplitgraphClusterArray[listOfModifiedVertices[elem]] = px;
				}
				// O RECALCULO DA MATRIZ ABAIXO EH FEITO DE FORMA INCREMENTAL, reutilizando a matrix processClusterImbMatrix
				// ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
				ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;
				// Realiza o calculo do delta da matriz de imbalance entre processos
				updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
						listOfModifiedVertices, tempProcessClusterImbMatrix);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";
				// Valida se a matriz incremental e a full sao iguais
				/*
				bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
				bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
				BOOST_LOG_TRIVIAL(info) << "*** As matrizes de delta sao iguais: " << igual << " e " << igual2;
				*/

				// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
				Clustering newGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem,
						info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix, bestClustering);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
				/*
				ClusterArray cTemp = newGlobalClustering.getClusterArray();
				Clustering validation(cTemp, *g, problem);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Full Obj Calc: I(P) = " << validation.getImbalance().getValue();
				if(validation.getImbalance().getValue() != newGlobalClustering.getImbalance().getValue()) {
					BOOST_LOG_TRIVIAL(error) << "[MoveCluster1opt] Obj functions do not match.";
				} */

				if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] moveCluster1opt Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
					bestClustering = newGlobalClustering;

					// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
					processClusterImbMatrix = tempProcessClusterImbMatrix;

					// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
					bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
					foundBetterSolution = true;
					break;
				} else {
					long sizeOfSourceCluster = bestSplitgraphClustering.getClusterSize(currentProcess);
					long sizeOfDestCluster = bestSplitgraphClustering.getClusterSize(px);
					long newSizeOfSourceCluster = sizeOfSourceCluster - listOfModifiedVertices.size();
					long newSizeOfDestCluster = sizeOfDestCluster + listOfModifiedVertices.size();

					if (newGlobalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
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

						// Realiza a movimentacao de um cluster especifico (cluster move) para outro processo
						// O valor da FO continua o mesmo
						bestClustering.setProcessOrigin(clusterToMove, px);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost moveCluster1opt improving solution found! I(P) = "
														<< bestClustering.getImbalance().getValue();

						// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
						processClusterImbMatrix = tempProcessClusterImbMatrix;

						// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
						bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
						foundBetterSolution = true;
						break;
					} else {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost moveCluster1opt move.";
					}
				}
			} else {
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected moveCluster1opt move.";
			}
		}
		if(foundBetterSolution) {
			break;
		}
		// d) Dúvida: rodar um novo ILS inteiro para cada movimentação ou apenas uma busca local?
		// e) Cada processo envia o resultado (novo clustering) para o mestre e este calcula se vale a pena manter a alteração
		// f) O cálculo feito pelo mestre baseia-se na variação do imbalance dentro de cada processo (delta I(P) int) mais o novo imbalance externo, entre os processos (I(P) ext)
		// 		(a variação dentro de cada processo já é fornecida por cada processo, via mensagem – diferença do I(P) anterior)
		// 		(a variação do imbalance entre os processos é calculada varrendo a vizinhança do vértice movido)
		// g) Se gerar um imbalance pior, o movimento pode ser rejeitado, gerando uma volta da configuração anterior, por isso cada processo deve guardar o clustering anterior
	}
	return foundBetterSolution;
}

bool ImbalanceSubgraphParallelILS::swapCluster1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   Swap  cluster  1 - opt  ***********************************
	// In this neighborhood the priority is swapping a bigger cluster in an overloaded process
	// 	  with a smaller cluster from a less-loaded process.
	// (priotitize vertex load balancing between processes)
	// Tries to move the top 25% of biggest clusters from all overloaded processes
	//   Vertex division between processes (has nothing to do with a solution to the CC problem)
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	//   Current best solution to the CC problem considering the whole graph (global solution)
	ClusterArray globalClusterArray = bestClustering.getClusterArray();

	std::vector<Coordinate> overloadedProcessList = obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** swapCluster1opt";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph swapCluster1opt] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	bool foundBetterSolution = false;
	long n = g->getN();
	long nc = bestClustering.getNumberOfClusters();
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
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
			std::vector<long> listOfMovedVerticesFromProcessA = getListOfVeticesInCluster(*g, bestClustering, clusterToSwapA);

			// finds out to which process the cluster belongs to
			int currentProcess = splitgraphClusterArray[listOfMovedVerticesFromProcessA[0]];
			assert(currentProcess == procSourceNum);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Swap from process " << procSourceNum <<
					": Moving global cluster " << clusterToSwapA << " of size " << bestClustering.getClusterSize(clusterToSwapA);

			// LOAD BALANCING: Try to move the cluster to a process with less vertices than a given threshold (less loaded process)
			ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
			// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
			// TODO possible parametrization here!
			long maxVerticesAllowedInProcess = (2 * n) / bestSplitgraphClustering.getNumberOfClusters();
			for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {
				if(procDestNum != procSourceNum) {
					std::vector<Coordinate> clusterListProcessB = obtainListOfClustersFromProcess(*g, bestClustering, procDestNum);
					// gets the cluster from process B (destination) which has less elements
					Coordinate clusterToSwapB = *std::min_element(clusterListProcessB.begin(), clusterListProcessB.end(), coordinate_ordering_asc());
					std::vector<long> listOfMovedVerticesFromProcessB = getListOfVeticesInCluster(*g, bestClustering, clusterToSwapB.x);

					// only makes the swap if the max vertex threshold is respected
					if (bestSplitgraphClustering.getClusterSize(procDestNum) + listOfMovedVerticesFromProcessA.size()
							- listOfMovedVerticesFromProcessB.size() < maxVerticesAllowedInProcess) {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to swap global cluster " << clusterToSwapA
								<< " to process " << procDestNum;
						ClusterArray previousSplitgraphClusterArray = currentSplitgraphClusterArray;
						ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
						ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;

						// SWAP Step 1: Move the vertices from a specific cluster from source process (cluster move from A to B)
						for(long elem = 0; elem < listOfMovedVerticesFromProcessA.size(); elem++) {
							tempSplitgraphClusterArray[listOfMovedVerticesFromProcessA[elem]] = procDestNum;
						}
						// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processClusterImbMatrix
						updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
								listOfMovedVerticesFromProcessA, tempProcessClusterImbMatrix);
						/*
						ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
						// Valida se a matriz incremental e a full sao iguais
						bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
						bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
						BOOST_LOG_TRIVIAL(info) << "Op1 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */
						previousSplitgraphClusterArray = tempSplitgraphClusterArray;

						// SWAP Step 2: Move the vertices from the smallest cluster in / from process B (dest), to process A (source)
						for(long elem = 0; elem < listOfMovedVerticesFromProcessB.size(); elem++) {
							tempSplitgraphClusterArray[listOfMovedVerticesFromProcessB[elem]] = procSourceNum;
						}
						// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processClusterImbMatrix
						updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
								listOfMovedVerticesFromProcessB, tempProcessClusterImbMatrix);
						/*
						ImbalanceMatrix tempProcessClusterImbMatrix3 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
						// Valida se a matriz incremental e a full sao iguais
						igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix3.pos, 0.001, 0.1);
						igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix3.neg, 0.001, 0.1);
						BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from swap done.";

						// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
						Clustering newGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem,
								info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix, bestClustering);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
						/*
						ClusterArray cTemp = newGlobalClustering.getClusterArray();
						Clustering validation(cTemp, *g, problem);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Full Obj Calc: I(P) = " << validation.getImbalance().getValue();
						if(validation.getImbalance().getValue() != newGlobalClustering.getImbalance().getValue()) {
							BOOST_LOG_TRIVIAL(error) << "[SwapCluster1opt] Obj functions do not match.";
						} */

						if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] swapCluster1opt Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
							bestClustering = newGlobalClustering;

							// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
							processClusterImbMatrix = tempProcessClusterImbMatrix;

							// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
							bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
							foundBetterSolution = true;
							break;
						} else {
							long sizeOfSourceProcess = bestSplitgraphClustering.getClusterSize(procSourceNum);
							long sizeOfDestProcess = bestSplitgraphClustering.getClusterSize(procDestNum);
							long newSizeOfSourceProcess = sizeOfSourceProcess - listOfMovedVerticesFromProcessA.size() + listOfMovedVerticesFromProcessB.size();
							long newSizeOfDestProcess = sizeOfDestProcess - listOfMovedVerticesFromProcessB.size() + listOfMovedVerticesFromProcessA.size();

							if (newGlobalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
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

							// => simply swaps the clusters between the processes, keeping the same global imbalance value
							// evaluate if this zero-cost move improves the vertex balancing between processes
							if (labs(newSizeOfSourceProcess - newSizeOfDestProcess) <
									labs(sizeOfSourceProcess - sizeOfDestProcess)) {  // this zero-cost move is good

								// Swaps the clusters between the processes
								bestClustering.setProcessOrigin(clusterToSwapB.x, procSourceNum);
								bestClustering.setProcessOrigin(clusterToSwapA, procDestNum);
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost swapCluster1opt improving solution found! I(P) = "
																<< bestClustering.getImbalance().getValue();

								// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
								processClusterImbMatrix = tempProcessClusterImbMatrix;

								// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
								bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
								foundBetterSolution = true;
								break;
							} else {
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost swapCluster1opt move.";
							}
						}
					} else {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] swapCluster1opt Movement rejected.";
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
		if(foundBetterSolution) {
			break;
		}
	}
	return foundBetterSolution;
}

bool ImbalanceSubgraphParallelILS::twoMoveCluster(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	// ************************   2-Move  cluster  ***********************************
	// b) The master process tries to move 2 clusters from an overloaded process to 2 less-loaded processes.
	// select a 2 clusters to move - two clusters from an overloaded process (process with too many vertices)
	// (priotitize vertex load balancing between processes)
	// for now it is the cluster c that belongs to the processPair (x, y) where x is an overloaded process
	// Tries to move the top 25% of clusters of overloaded processes
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	std::vector<Coordinate> clusterImbalanceList = obtainListOfImbalancedClusters(*g, splitgraphClusterArray, bestClustering);
	bool foundBetterSolution = false;
	long n = g->getN();
	long nc = bestClustering.getNumberOfClusters();

	std::vector<Coordinate> overloadedProcessList = obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** twoMoveCluster";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph twoMoveCluster] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
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
					std::vector<long> listOfMovedVerticesFromClusterA = getListOfVeticesInCluster(*g, bestClustering, clusterToMoveA);
					std::vector<long> listOfMovedVerticesFromClusterB = getListOfVeticesInCluster(*g, bestClustering, clusterToMoveB);

					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2-Move from process " << procSourceNum << ": Moving global clusters "
							<< clusterToMoveA << " and " <<	clusterToMoveB << ".";

					// LOAD BALANCING: Try to move the clusters to 2 processes with less vertices than a given threshold (less loaded process)
					ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
					// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
					// TODO possible parametrization here!
					long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0); // (2 * n) / bestSplitgraphClustering.getNumberOfClusters();
					for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {  // destination process 1
						if((procDestNum != procSourceNum) and (bestSplitgraphClustering.getClusterSize(procDestNum) +
								listOfMovedVerticesFromClusterA.size() < maxVerticesAllowedInProcess)) {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Trying to move global cluster " << clusterToMoveA
									<< " to process " << procDestNum;
							for(int procDestNum2 = 0; procDestNum2 < numberOfProcesses; procDestNum2++) {  // destination process 2
								if((procDestNum2 != procSourceNum) and (procDestNum2 != procDestNum)
										and (bestSplitgraphClustering.getClusterSize(procDestNum2) +
										listOfMovedVerticesFromClusterB.size() < maxVerticesAllowedInProcess)) {
									BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] And trying to move global cluster " << clusterToMoveB
												<< " to process " << procDestNum2;
									ClusterArray previousSplitgraphClusterArray = currentSplitgraphClusterArray;
									ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
									ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;

									// MOVE Step 1: Move the vertices from a specific cluster A from source process to dest process 1
									for(long elem = 0; elem < listOfMovedVerticesFromClusterA.size(); elem++) {
										tempSplitgraphClusterArray[listOfMovedVerticesFromClusterA[elem]] = procDestNum;
									}
									// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processC
									updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
											listOfMovedVerticesFromClusterA, tempProcessClusterImbMatrix);
									previousSplitgraphClusterArray = tempSplitgraphClusterArray;
									// Valida se a matriz incremental e a full sao iguais
									/*
									ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
									bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
									bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
									BOOST_LOG_TRIVIAL(info) << "Op1 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

									// MOVE Step 2: Move the vertices from a specific cluster B from source process to dest process 2
									for(long elem = 0; elem < listOfMovedVerticesFromClusterB.size(); elem++) {
										tempSplitgraphClusterArray[listOfMovedVerticesFromClusterB[elem]] = procDestNum2;
									}
									// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix process
									updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
											listOfMovedVerticesFromClusterB, tempProcessClusterImbMatrix);
									// Valida se a matriz incremental e a full sao iguais
									/*
									ImbalanceMatrix tempProcessClusterImbMatrix3 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
									igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix3.pos, 0.001, 0.1);
									igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix3.neg, 0.001, 0.1);
									BOOST_LOG_TRIVIAL(info) << "Op2 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

									BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from 2-move done.";

									// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
									Clustering newGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax,
											problem, info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix, bestClustering);
									BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
									/*
									ClusterArray cTemp = newGlobalClustering.getClusterArray();
									Clustering validation(cTemp, *g, problem);
									BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Full Obj Calc: I(P) = " << validation.getImbalance().getValue();
									if(validation.getImbalance().getValue() != newGlobalClustering.getImbalance().getValue()) {
										BOOST_LOG_TRIVIAL(error) << "[2-move-cluster] Obj functions do not match.";
									} */

									if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
										BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] twoMoveCluster Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
										bestClustering = newGlobalClustering;

										// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
										processClusterImbMatrix = tempProcessClusterImbMatrix;

										// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
										bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
										foundBetterSolution = true;
										break;
									} else {
										long sizeOfSourceProcess = bestSplitgraphClustering.getClusterSize(procSourceNum);
										long sizeOfDestProcessA = bestSplitgraphClustering.getClusterSize(procDestNum);
										long sizeOfDestProcessB = bestSplitgraphClustering.getClusterSize(procDestNum2);
										long newSizeOfSourceProcess = sizeOfSourceProcess - listOfMovedVerticesFromClusterA.size() - listOfMovedVerticesFromClusterB.size();
										long newSizeOfDestProcessA = sizeOfDestProcessA + listOfMovedVerticesFromClusterA.size();
										long newSizeOfDestProcessB = sizeOfDestProcessB + listOfMovedVerticesFromClusterB.size();

										if (newGlobalClustering.getImbalance().getValue() == bestClustering.getImbalance().getValue()) {
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

										// => simply swaps the clusters between the processes, keeping the same global imbalance value
										// evaluate if this zero-cost move improves the vertex balancing between processes
										if (labs(newSizeOfSourceProcess - (sizeOfDestProcessA + sizeOfDestProcessB)) <
												labs(sizeOfSourceProcess - (newSizeOfDestProcessA + newSizeOfDestProcessB))) {  // this zero-cost move is good

											// Moves the clusters to the 2 processes
											bestClustering.setProcessOrigin(clusterToMoveA, procDestNum);
											bestClustering.setProcessOrigin(clusterToMoveB, procDestNum2);
											BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost twoMoveCluster improving solution found! I(P) = "
																			<< bestClustering.getImbalance().getValue();

											// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
											processClusterImbMatrix = tempProcessClusterImbMatrix;

											// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
											bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
											foundBetterSolution = true;
											break;
										} else {
											BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost twoMoveCluster move.";
										}
									}
								} else {
									BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] twoMoveCluster Movement rejected.";
								}
								if(foundBetterSolution) {
									break;
								}
							}
							if(foundBetterSolution) {
								break;
							}
						} else {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] twoMoveCluster Movement rejected.";
						}
					}
					if(foundBetterSolution) {
						break;
					}
				}
			}
		}
	}
	return foundBetterSolution;
}

long ImbalanceSubgraphParallelILS::variableNeighborhoodDescent(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info, const double& timeSpentSoFar) {

	int r = 1, iteration = 0, l = 4;
	double timeSpentOnLocalSearch = 0.0;
	long improvedOnVND = 0;
	BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Global VND local search...";
	std::vector<long> improvementStats(l + 1, 0);
	std::vector<long> neigExecStats(l + 1, 0);
	// stores the list of improved solutions found at each VND iteration (global clustering and splitgraph clustering)
	std::vector< std::pair<Clustering, Clustering> > solutionHistory;
	solutionHistory.push_back( std::make_pair(bestClustering, bestSplitgraphClustering) );
	// stores the process balancing (number of clusters, vertices) info and solution values throughout the time
	std::vector<string> splitgraphVNDProcessBalancingHistory;
	stringstream ss;
	ss << iteration++ << ", " << timeSpentOnLocalSearch << ", " << bestClustering.getImbalance().getValue();
	splitgraphVNDProcessBalancingHistory.push_back(ss.str());
	numberOfFrustratedSolutions = 0;

	while (r <= l && (timeSpentSoFar + timeSpentOnLocalSearch < vnd->getTimeLimit())) {
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();
		BOOST_LOG_TRIVIAL(debug)<< "[Splitgraph] Global VND neighborhood " << r << " ...";

		rebalanceClustersBetweenProcessesWithZeroCost(g, problem, bestSplitgraphClustering, bestClustering,
				numberOfProcesses, processClusterImbMatrix);

		bool improved = false;
		if(r == 1) {
			// ************************   Move   1 - opt   cluster ***********************************
			improved = moveCluster1opt(g, bestSplitgraphClustering, bestClustering,
									numberOfProcesses, processClusterImbMatrix,
									construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} else if(r == 4) {
			// ************************   Split  Cluster  Move ***********************************
			improved = splitClusterMove(g, bestSplitgraphClustering, bestClustering,
									numberOfProcesses, processClusterImbMatrix,
									construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} else if(r == 2) {
			// ************************   2-Move   cluster ***********************************
			improved = twoMoveCluster(g, bestSplitgraphClustering, bestClustering,
									numberOfProcesses, processClusterImbMatrix,
									construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} else if(r == 3) {
			// ************************   Swap   1 - opt   cluster ***********************************
			improved = swapCluster1opt(g, bestSplitgraphClustering, bestClustering,
									numberOfProcesses, processClusterImbMatrix,
									construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} /* else if(r == 5) {  DISABLED NEIGHBORHOOD
			// ************************   Move   1 - opt   vertex  ***********************************
			improved = moveVertex1opt(g, bestSplitgraphClustering, bestClustering,
								numberOfProcesses, processClusterImbMatrix,
								construct, vnd, iter, iterMaxILS, perturbationLevelMax, problem, info);
		} */
		// records statistics for neighborhood structures execution
		neigExecStats[r]++;
		if(improved)  improvementStats[r]++;

		// => Finally: Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timeSpentOnLocalSearch += (end_time.wall - start_time.wall)	/ double(1000000000);

		if(bestClustering.getImbalance().getValue() < 0.0) {
			BOOST_LOG_TRIVIAL(error)<< "[Splitgraph] Objective function below zero. Obj = " << bestClustering.getImbalance().getValue();
			break;
		}
		if(improved) {
			// BOOST_LOG_TRIVIAL(trace) << myRank << ": New local solution found: " << setprecision(2) << il.getValue() << endl;
			r = 1;
			improvedOnVND++;
			solutionHistory.push_back( std::make_pair(bestClustering, bestSplitgraphClustering) );
			stringstream ss;
			ss << iteration << ", " << timeSpentOnLocalSearch << ", " << bestClustering.getImbalance().getValue();
			splitgraphVNDProcessBalancingHistory.push_back(ss.str());
			if(bestClustering.getImbalance().getValue() <= EPS) {
				BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Stopping global VND local search, since I(P) <= 0.";
				break;
			}
		} else {  // no better result found in neighborhood
			r++;
			// BOOST_LOG_TRIVIAL(debug) << "Changed to neighborhood size l = " << k;
		}
		iteration++;
	}
	BOOST_LOG_TRIVIAL(info)<< "[Splitgraph] Global VND local search done. Obj = " << bestClustering.getImbalance().getValue() <<
			". Time spent: " << timeSpentOnLocalSearch << " s.";

	stringstream splitgraphVNDResults;
	BOOST_LOG_TRIVIAL(info) << "[Splitgraph] Number of frustrated ILS solutions: " << numberOfFrustratedSolutions;
	splitgraphVNDResults << "[Splitgraph] Number of frustrated ILS solutions: " << numberOfFrustratedSolutions << "\n";
	BOOST_LOG_TRIVIAL(info) << "[Splitgraph] Execution and Improvement statistics for each neighborhood: ";
	splitgraphVNDResults << "[Splitgraph] Execution and Improvement statistics for each neighborhood: \n";
	BOOST_LOG_TRIVIAL(info) << "Neighborhood size, Num exec, Num improv";
	splitgraphVNDResults << "Neighborhood size, Num exec, Num improv\n";
	for(int i = 1; i <= l; i++) {
		BOOST_LOG_TRIVIAL(info) << "r = " << i << ", " << neigExecStats[i] << ", " << improvementStats[i];
		splitgraphVNDResults << "r = " << i << ", " << neigExecStats[i] << ", " << improvementStats[i] << "\n";
	}
	splitgraphVNDResults << "[Splitgraph] Best solution found: ";

	std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	long n = g->getN();
	for(long i = 0; i < n; i++) {
		long k = splitgraphClusterArray[i];
		verticesInCluster[k].push_back(i);
	}
	for(int p = 0; p < numberOfProcesses; p++) {
		SignedGraph sg(g->graph, verticesInCluster[p]);
		std::vector<Coordinate> clusterList = obtainListOfClustersFromProcess(*g, bestClustering, p);

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_edges = " <<
				sg.getM() << " , num_vertices = " << sg.getN()
				<< ", k = " << clusterList.size();
		splitgraphVNDResults << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_edges = " <<
				sg.getM() << " , num_vertices = " << sg.getN()
				<< ", k = " << clusterList.size() << "\n";
	}
	splitgraphVNDResults << "[Parallel ILS SplitGraph] I(P) = " << bestClustering.getImbalance().getValue() << "\n";
	splitgraphVNDResults << "[Parallel ILS SplitGraph] Time spent = " << timeSpentOnLocalSearch << " s.\n";

	// Saves the splitgraph statistics to csv file
	generateOutputFile(problem, splitgraphVNDResults, info.outputFolder, info.fileId, info.executionId,
			info.processRank, string("statistics-splitgraph"), construct->getAlpha(), l, iter);

	// Exports the splitgraph solutions to csv file
	stringstream sshistory;
	sshistory << "Iteration, Time, I(P)";
	for(int i = 1; i <= numberOfProcesses; i++) {
		sshistory << ", e" << i << ", n" << i << ", k" << i;
	}
	sshistory << "\n";
	for(int iternum = 0; iternum < solutionHistory.size(); iternum++) {
		sshistory << splitgraphVNDProcessBalancingHistory[iternum];

		Clustering globalClustering = solutionHistory[iternum].first;
		Clustering splitgraphClustering = solutionHistory[iternum].second;
		std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
		ClusterArray splitgraphClusterArray = splitgraphClustering.getClusterArray();
		for(long i = 0; i < n; i++) {
			long k = splitgraphClusterArray[i];
			verticesInCluster[k].push_back(i);
		}
		for(int p = 0; p < numberOfProcesses; p++) {
			SignedGraph sg(g->graph, verticesInCluster[p]);
			std::vector<Coordinate> clusterList = obtainListOfClustersFromProcess(*g, globalClustering, p);
			sshistory << ", " << sg.getM() << ", " << sg.getN() << ", " << clusterList.size();
		}
		sshistory << "\n";
	}
	// includes the biggest cluster size of each process in the best solution
	Clustering globalClustering = solutionHistory.back().first;
	Clustering splitgraphClustering = solutionHistory.back().second;
	for(int i = 1; i <= numberOfProcesses; i++) {
		sshistory << "verticesIn(P" << i << "), " << "biggestClusterSize(P" << i << "), ";
	}
	sshistory << "\n";
	for(int p = 0; p < numberOfProcesses; p++) {
		SignedGraph sg(g->graph, verticesInCluster[p]);
		std::vector<Coordinate> bigClustersList = obtainListOfClustersFromProcess(*g, globalClustering, p);
		std::vector<Coordinate>::iterator biggestCluster = std::max_element(bigClustersList.begin(), bigClustersList.end(), coordinate_ordering_asc());
		sshistory << sg.getN() << ", " << boost::math::lround(biggestCluster->value) << ", ";
	}
	sshistory << "\n";
	// includes the number of frustrated solutions
	sshistory << "Frustrated ILS solutions, " << numberOfFrustratedSolutions << "\n";

	// Saves the splitgraph solution history to csv file
	generateOutputFile(problem, sshistory, info.outputFolder, info.fileId, info.executionId,
			info.processRank, string("solutionHistory-splitgraph"), construct->getAlpha(), l, iter);

	// return improvedOnVND;
	return 0;  // TODO for now, we only run global VND once
}

std::vector<long> ImbalanceSubgraphParallelILS::getListOfVeticesInCluster(SignedGraph& g, Clustering& globalClustering,
		long clusterNumber) {

	long n = g.getN();
	ClusterArray globalClusterArray = globalClustering.getClusterArray();
	std::vector<long> vertexList;
	for(long vx = 0; vx < n; vx++) {
		if(globalClusterArray[vx] == clusterNumber) {
			//BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Choosing vertex " << vx << " to be moved, it belongs to cluster" << clusterToMove << ".";
			vertexList.push_back(vx);
		}
	}
	return vertexList;
}


/**
  *  Identifies a pseudo clique C+ inside an overloaded cluster, and tries to
  *  move this clique to another (not overloaded) process as a new cluster.
*/
bool ImbalanceSubgraphParallelILS::splitClusterMove(SignedGraph* g, Clustering& bestSplitgraphClustering,
		Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	std::vector<Coordinate> overloadedProcessList = obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** splitClusterMove";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph splitClusterMove] Found " << overloadedProcessList.size()  << " overloaded processes.";
	// sorts the list of overloaded processes according to the number of vertices, descending order
	std::sort(overloadedProcessList.begin(), overloadedProcessList.end(), coordinate_ordering_desc());
	bool foundBetterSolution = false;
	long n = g->getN();
	long nc = bestClustering.getNumberOfClusters();
	for(std::vector<Coordinate>::iterator prIter = overloadedProcessList.begin(); prIter != overloadedProcessList.end(); prIter++) {
		int procSourceNum = prIter->x;  // overloaded process number
		BOOST_LOG_TRIVIAL(info) << "SplitClusterMove for overloaded process number " << procSourceNum <<
				" with " << prIter->value << " vertices...";

		// obtains the list of big clusters from the overloaded process 'procNum'
		std::vector<Coordinate> bigClustersList = obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
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
			std::vector<long> cliqueC = findPseudoCliqueC(g, bestClustering, clusterX);

			// try to move the entire clique C+ (as a new cluster) to another process B
			// as long as process B is not overloaded
			// finds out to which process the cluster belongs to
			int currentProcess = splitgraphClusterArray[cliqueC.front()];
			assert(currentProcess == procSourceNum);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Move pseudo clique of size " << cliqueC.size()
					<<  " from process " << procSourceNum <<
					": Moving from global cluster " << clusterX << " of size " << bestClustering.getClusterSize(clusterX);
			assert(bigClustersList[clusterCount].y == procSourceNum);
			assert(currentProcess == procSourceNum);

			// LOAD BALANCING: Try to move the clique to a process with less vertices than a given threshold (less loaded process)
			ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
			// the maximum number of vertices allowed in a process is 2 times the initial number of vertices of a process
			// TODO possible parametrization here!
			long maxVerticesAllowedInProcess = (2 * n) / bestSplitgraphClustering.getNumberOfClusters();
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
					// ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
					ImbalanceMatrix tempProcessClusterImbMatrix = processClusterImbMatrix;
					// Realiza o calculo do delta da matriz de imbalance entre processos
					updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
							cliqueC, tempProcessClusterImbMatrix);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix done.";
					// Valida se a matriz incremental e a full sao iguais
					/*
					bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
					bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
					BOOST_LOG_TRIVIAL(info) << "*** As matrizes de delta sao iguais: " << igual << " e " << igual2;
					*/

					// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
					Clustering newGlobalClustering = runDistributedILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem,
							info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix, bestClustering);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
					/*
					ClusterArray cTemp = newGlobalClustering.getClusterArray();
					Clustering validation(cTemp, *g, problem);
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Full Obj Calc: I(P) = " << validation.getImbalance().getValue();
					if(validation.getImbalance().getValue() != newGlobalClustering.getImbalance().getValue()) {
						BOOST_LOG_TRIVIAL(error) << "[MoveCluster1opt] Obj functions do not match.";
					} */

					if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] splitClusterMove Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
						bestClustering = newGlobalClustering;

						// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
						processClusterImbMatrix = tempProcessClusterImbMatrix;

						// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
						bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
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
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost move.";
							// WARNING: It is possible that the new solution has worse imbalance than it should, since
							//    the local ILS of one of the processes may not find the most efficient local solution.
							//    In these cases, discards the clustering solutions of each process and simply creates
							//    a zero-cost move solution (i.e. same imbalance than previous) based on the moved vertices
							//    and the previous best-known clustering.

							// => simply swaps the clusters between the processes, keeping the same global imbalance value
							// evaluate if this zero-cost move improves the vertex balancing between processes
							if (labs(newSizeOfSourceProcess - newSizeOfDestProcess) <
									labs(sizeOfSourceProcess - sizeOfDestProcess)) {  // this zero-cost move is good

								bestClustering = newGlobalClustering;
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Zero-cost splitClusterMove improving solution found! I(P) = "
																<< bestClustering.getImbalance().getValue();

								// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
								processClusterImbMatrix = tempProcessClusterImbMatrix;

								// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
								bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
								foundBetterSolution = true;
								break;
							} else {
								BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Rejected zero-cost splitClusterMove move.";
							}
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

/**
  *  Returns a list containing the vertices that belong to the positive clique C+,
  *  found inside global cluster X. This is a greedy heuristic.
*/
std::vector<long> ImbalanceSubgraphParallelILS::findPseudoCliqueC(SignedGraph *g, Clustering& globalClustering,
		long clusterX) {

	string strCluster = boost::lexical_cast<std::string>(clusterX);
	BOOST_LOG_TRIVIAL(info) << "Invoking findPseudoCliqueC for global cluster number " << strCluster << "...";
	long n = g->getN();
	long m = g->getM();

	// *** A - Constructive phase: building a maximal clique
	//  Finding a maximal clique is straightforward: Starting with an arbitrary clique
	//  (for instance, a single vertex), grow the current clique one vertex at a time
	//  by iterating over the graph’s remaining vertices, adding a vertex if it is connected
	//  to each vertex in the current clique, and discarding it otherwise.
	//  This algorithm runs in linear time.
	// 1. For every vertex v in clusterX, let Dx_+(v) = number of positive neighbors inside clusterX.
	//   1.1. Obtains a list containing vertex degree info (vertex number and its positive and negative degrees)
	std::list<VertexDegree> degreesInsideClusterX = calculateDegreesInsideCluster(g, globalClustering, clusterX);
	//   1.2. Sorts the vertex list according to a weighted formula of positive and negative degrees (descending order)
	degreesInsideClusterX.sort(weighted_pos_neg_degree_ordering_desc());

	// 2. The positive clique C+ begins with the vertex v that has the biggest Dx_+(v) in clusterX
	//  Two data structures are used to control the vertices inside cliqueC: a vector and a binary cluster array
	std::list<long> cliqueC;
	ClusterArray cliqueCClusterArray(g->getN(), Clustering::NO_CLUSTER);
	long u = this->chooseRandomVertex(degreesInsideClusterX, boost::math::lround(CLUSTERING_ALPHA * degreesInsideClusterX.size()));
	cliqueC.push_back(u);
	cliqueCClusterArray[u] = 1;
	BOOST_LOG_TRIVIAL(info) << "Global cluster number " << strCluster << " currently has size " << degreesInsideClusterX.size();
	// BOOST_LOG_TRIVIAL(info) << "Building cliqueC with initial vertex " << u << "...";

	// 3. At each step, add to the clique C+ the vertex u (that also belongs to cluster X), where
	//    C+ \union {u} is also a positive clique, and the vertex u is a random vertex between
	//    the first (alpha x ||Dx_+(v)||) in clusterX => randomness factor
	long cliqueSize = 1;
	while (degreesInsideClusterX.size() > 0) { // list != empty
		// Choose vertex u randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		u = this->chooseRandomVertex(degreesInsideClusterX, boost::math::lround(CLUSTERING_ALPHA * degreesInsideClusterX.size()));
		// only adds u to the cliqueC if C+ \union {u} is also a positive clique
		if(isPseudoClique(g, cliqueC, cliqueCClusterArray, u)) {
			cliqueC.push_back(u);
			cliqueCClusterArray[u] = 1;
			cliqueSize++;
			// BOOST_LOG_TRIVIAL(info) << "Adding vertex " << it->id << " to the cliqueC.";
		}
		// limits the clique size to half of the source cluster size
		if(cliqueSize >= boost::math::lround(degreesInsideClusterX.size() / (double)2.0))  break;
	}
	BOOST_LOG_TRIVIAL(info) << "Built initial cliqueC of size " << cliqueC.size() << ".";
	// 4. The constructive procedure stops when the clique C+ is maximal.
	// That is, the list of vertices in clusterX has been fully traversed.

	// *** B - Unique local search (executed only once)
	// 1. For every vertex v in clique C+, test a 1-2-swap neighbourhood, that is, try to
	//    remove 1 vertex from C+ and insert other 2 vertices that belong to clusterX (but still not to C+),
	//    preserving the property of maximal positive clique.
	/*
	BOOST_LOG_TRIVIAL(info) << "CliqueC local search...";
	bool improvement = false;
	int iter = 0;
	std::vector<long> verticesInClusterX = getListOfVeticesInCluster(*g, globalClustering, clusterX);
	do {
		improvement = false;
		iter++;
		// BOOST_LOG_TRIVIAL(info) << "Clique local search iteration " << iter << "...";
		for(std::list<long>::iterator it = cliqueC.begin(); it != cliqueC.end(); it++) {
			long v = *it;
			// try to remove vertex v from tempCliqueC
			// BOOST_LOG_TRIVIAL(info) << "Trying to remove vertex " << v << " from the clique.";
			cliqueCClusterArray[v] = Clustering::NO_CLUSTER;

			// 1-2-swap neighbourhood:
			// try to insert 2 vertices that belong to clusterX but are still not in cliqueC (and are not the vertex v removed)
			//  for every combination of 2 vertices (a and b) in clusterX but not in cliqueC:
			for(std::vector<long>::iterator it_a = verticesInClusterX.begin(); it_a != verticesInClusterX.end(); it_a++) {
				long a = *it_a;
				if((cliqueCClusterArray[a] < 0) and (a != v)) {  // vertex a belongs to clusterX but is not in cliqueC
					for(std::vector<long>::iterator it_b = verticesInClusterX.begin(); it_b != verticesInClusterX.end(); it_b++) {
						long b = *it_b;
						if((cliqueCClusterArray[b] < 0) and (b != v) and (b < a)) {  // vertex b belongs to clusterX but is not in cliqueC
							// check if cliqueC + {a, b} is still a positive clique
							if(isPseudoClique(g, cliqueCClusterArray, cliqueC.size() - 1, a, b)) {
								// finally insert vertices a and b in cliqueC
								// BOOST_LOG_TRIVIAL(info) << "Removed vertex " << v << " from the clique.";
								// BOOST_LOG_TRIVIAL(info) << "Added vertices " << a << " and " << b  << " to the clique.";

								cliqueCClusterArray[a] = 1;
								cliqueCClusterArray[b] = 1;
								cliqueC.push_back(a);
								cliqueC.push_back(b);
								cliqueC.erase(std::remove(cliqueC.begin(), cliqueC.end(), v), cliqueC.end());
								// BOOST_LOG_TRIVIAL(info) << "CliqueC is now of size " << cliqueC.size() << ".";
								improvement = true;
								break;
							}
						}
						if(improvement)  break;
					}
				}
			}
			if(improvement){
				break;  // restarts the vertex v in cliqueC loop since cliqueC was modified
				// TODO replace this line with a time limit verification (10s?)
			} else {  // undo the movement
				cliqueCClusterArray[v] = 1;
			}
		}
	} while(improvement);
	// 2. Repeat this process until C+ has reached a local optimum related to the 1-2-swap neighbourhood.
	BOOST_LOG_TRIVIAL(info) << "Found a pseudoCliqueC of size " << cliqueC.size() << ".";
*/
	// converts the list to a vector of long
	std::vector<long> cliqueCvec(cliqueC.size());
	copy(cliqueC.begin(), cliqueC.end(), cliqueCvec.begin());
	return cliqueCvec;
}

/**
  *  Determines if the set cliqueC \union {u} is a pseudo clique in the graph g.
  *  By definition, cliqueC is already a pseudo clique.
*/
bool ImbalanceSubgraphParallelILS::isPseudoClique(SignedGraph *g, std::list<long>& cliqueC,
		ClusterArray& cliqueCClusterArray, long u) {

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
	UndirectedGraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCount = 0, positiveEdgeCount = 0, negativeEdgeCount = 0;

	// validates the edges from vertex u to cliqueC
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(u, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if(cliqueCClusterArray[j] >= 0) {
			totalEdgeCount++;
			if(ew[e].weight > 0) {
				positiveEdgeCount++;
			} else {
				negativeEdgeCount++;
			}
		}
	}
	// Must allow at least POS_EDGE_PERC_RELAX % of positive connections
	// Must allow no more than NEG_EDGE_PERC_RELAX % of negative connections
	double percOfPositiveConnections = positiveEdgeCount / (double)totalEdgeCount;
	double percOfNegativeConnections = negativeEdgeCount / (double)totalEdgeCount;
	bool isPseudoClique = ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
	return isPseudoClique;
}

/**
  *  Determines if the set cliqueC \union {u, v} is a pseudo clique in the graph g.
  *  By definition, cliqueC is already a pseudo clique.
*/
bool ImbalanceSubgraphParallelILS::isPseudoClique(SignedGraph *g, ClusterArray& cliqueCClusterArray,
		long cliqueSize, long u, long v) {

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
	UndirectedGraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCountU = 0, positiveEdgeCountU = 0, negativeEdgeCountU = 0;
	long totalEdgeCountV = 0, positiveEdgeCountV = 0, negativeEdgeCountV = 0;

	// validates the edges from vertex u to cliqueC
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(u, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if (cliqueCClusterArray[j] >= 0) {
			totalEdgeCountU++;
			if(ew[e].weight < 0) {
				negativeEdgeCountU++;
			} else {
				positiveEdgeCountU++;
			}
		}
	}
	double percOfPositiveConnections = positiveEdgeCountU / (double)totalEdgeCountU;
	double percOfNegativeConnections = negativeEdgeCountU / (double)totalEdgeCountU;
	bool isPseudoClique = ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
	if( not isPseudoClique ) {
		return false;
	}

	// validates the edges from vertex v to cliqueC
	for (boost::tie(f2, l2) = out_edges(v, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if (cliqueCClusterArray[j] >= 0) {
			totalEdgeCountV++;
			if(ew[e].weight < 0) {
				negativeEdgeCountV++;
			} else {
				positiveEdgeCountV++;
			}
		}
	}
	percOfPositiveConnections = positiveEdgeCountV / (double)totalEdgeCountV;
	percOfNegativeConnections = negativeEdgeCountV / (double)totalEdgeCountV;
	return ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
}

/**
  *  Calculates the positive and negative degrees of each vertex v in clusterX *relative to clusterX only*.
*/
std::list<VertexDegree> ImbalanceSubgraphParallelILS::calculateDegreesInsideCluster(SignedGraph *g,
		Clustering& globalClustering, long clusterX) {

	std::vector<long> verticesInClusterX = getListOfVeticesInCluster(*g, globalClustering, clusterX);
	ClusterArray myCluster(g->getN(), Clustering::NO_CLUSTER);
	for(std::vector<long>::iterator it = verticesInClusterX.begin(); it != verticesInClusterX.end(); it++) {
		long v = *it;
		myCluster[v] = 1;
	}
	std::list<VertexDegree> degreesInsideClusterX;
	for(std::vector<long>::iterator it = verticesInClusterX.begin(); it != verticesInClusterX.end(); it++) {
		long v = *it;
		// calculates the positive and negative degrees of vertex v *relative to clusterX only*
		double negDegree = fabs(g->getNegativeEdgeSumBetweenVertexAndClustering(v, myCluster));
		double posDegree = g->getPositiveEdgeSumBetweenVertexAndClustering(v, myCluster);
		degreesInsideClusterX.push_back(VertexDegree(v, posDegree, negDegree));
	}
	return degreesInsideClusterX;
}

long ImbalanceSubgraphParallelILS::chooseRandomVertex(std::list<VertexDegree>& vertexList, long x) {
	// Generates a random number between 1 and x
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	unsigned int selectedVertexSetIndex = randomUtil.next(0, x - 1);
	VertexDegree selectedVertex;

	// Finds the Vertex
	std::list<VertexDegree>::iterator pos = vertexList.begin();
	std::advance(pos, selectedVertexSetIndex);
	selectedVertex = *pos;
	vertexList.erase(pos);
	return selectedVertex.id;
}

void ImbalanceSubgraphParallelILS::rebalanceClustersBetweenProcessesWithZeroCost(SignedGraph* g, ClusteringProblem& problem,
		Clustering& bestSplitgraphClustering, Clustering& bestClustering, const int& numberOfProcesses,
		ImbalanceMatrix& processClusterImbMatrix) {

	ClusterArray splitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
	ClusterArray globalClusterArray = bestClustering.getClusterArray();
	bool foundBetterSolution = false;
	long n = g->getN();
	long nc = bestClustering.getNumberOfClusters();
	ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();

	std::vector<Coordinate> overloadedProcessList = obtainListOfOverloadedProcesses(*g, bestSplitgraphClustering);
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] ***************** rebalanceClustersBetweenProcessesWithZeroCost";
	BOOST_LOG_TRIVIAL(debug) << "[SplitGraph twoMoveCluster] Found " << overloadedProcessList.size()  << " overloaded processes.";
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
		long maxVerticesAllowedInProcess = boost::math::lround(n / (double)2.0); // (2 * n) / bestSplitgraphClustering.getNumberOfClusters();
		std::vector<long> clustersToMove;
		std::vector<Coordinate> fullClusterList = obtainListOfClustersFromProcess(*g, bestClustering, procSourceNum);
		// sorts the list of clusters according to the number of vertices, ascending order
		std::sort(fullClusterList.begin(), fullClusterList.end(), coordinate_ordering_asc());
		long numberOfVerticesToMove = 0;

		for(std::vector<Coordinate>::iterator iter = fullClusterList.begin(); iter != fullClusterList.end(); iter++) {
			long k = iter->x;
			long clusterSize = bestClustering.getClusterSize(k);
			if((clusterSize >= maxVerticesAllowedInProcess) or
					(numberOfVerticesInProcess - numberOfVerticesToMove < maxVerticesAllowedInProcess)) {
				break;
			}
			clustersToMove.push_back(k);
			numberOfVerticesToMove += clusterSize;
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Moving " << numberOfVerticesToMove
									<< " vertices from process " << procSourceNum;
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] After move, process " << procSourceNum << " will remain with "
				<< (numberOfVerticesInProcess - numberOfVerticesToMove) << " vertices.";

		for(std::vector<long>::iterator iter = clustersToMove.begin(); iter != clustersToMove.end(); iter++) {
			// finds a cluster in procSourceNum that can be moved to another process
			long clusterA = *iter;
			std::vector<long> listOfMovedVerticesFromClusterA = getListOfVeticesInCluster(*g, bestClustering, clusterA);

			// moves the cluster to a potential destination process
			for(int procDestNum = 0; procDestNum < numberOfProcesses; procDestNum++) {  // destination process 1
				if((procDestNum != procSourceNum) and (bestSplitgraphClustering.getClusterSize(procDestNum) +
						listOfMovedVerticesFromClusterA.size() < maxVerticesAllowedInProcess)) {
					BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Moving global cluster " << clusterA
							<< " to process " << procDestNum;

					// MOVE Step 1: Move the vertices from a specific cluster A from source process to dest process 1
					for(long elem = 0; elem < listOfMovedVerticesFromClusterA.size(); elem++) {
						currentSplitgraphClusterArray[listOfMovedVerticesFromClusterA[elem]] = procDestNum;
					}
					// INCREMENTAL UPDATE OF THE CLUSTER IMBALANCE MATRIX, reusing the matrix processC
					/*
					updateProcessToProcessImbalanceMatrix(*g, previousSplitgraphClusterArray, tempSplitgraphClusterArray,
							listOfMovedVerticesFromClusterA, tempProcessClusterImbMatrix);
					previousSplitgraphClusterArray = tempSplitgraphClusterArray; */
					// Valida se a matriz incremental e a full sao iguais
					/*
					ImbalanceMatrix tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
					bool igual = ublas::equals(tempProcessClusterImbMatrix.pos, tempProcessClusterImbMatrix2.pos, 0.001, 0.1);
					bool igual2 = ublas::equals(tempProcessClusterImbMatrix.neg, tempProcessClusterImbMatrix2.neg, 0.001, 0.1);
					BOOST_LOG_TRIVIAL(info) << "Op1 *** As matrizes de delta sao iguais: " << igual << " e " << igual2; */

					// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix from move done.";

					// Moves the clusters to the destination process -> does NOT change the imbalance
					bestClustering.setProcessOrigin(clusterA, procDestNum);

					// TODO validar o calculo da FO AQUI! e estudar calculo incremental aqui tb
					bestSplitgraphClustering = Clustering(currentSplitgraphClusterArray, *g, problem);
					// the movement for this cluster is done, proceed to the next one
					break;
				}
			}  // next destination process to be tested
		}  // next cluster to be moved
	}  // next overloaded processor
	// recalculates the process-process imbalance matrix
	processClusterImbMatrix = calculateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray);
}

} /* namespace ils */
} /* namespace resolution */
