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
		splitGraph(split), cudaEnabled(cuda) {

}

ImbalanceSubgraphParallelILS::~ImbalanceSubgraphParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering ImbalanceSubgraphParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}

	stringstream constructivePhaseResults;
	stringstream iterationResults;

	// 1. Divide the graph into (numberOfSlaves + 1) non-overlapping parts
	// Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned int numberOfProcesses = numberOfSlaves + 1;
	string strNP = boost::lexical_cast<std::string>(numberOfProcesses);
	BOOST_LOG_TRIVIAL(info) << "Invoking partitioning for spit graph (numberOfProcesses = " << strNP << ")...";

	// New split graph partitioning here
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
	// initializes a vector containing vertex degree info (vertex number and its pos and neg degrees)
	std::vector<VertexDegree> Gr;
	for(long i = 0; i < n; i++) {
		Gr.push_back(VertexDegree(i, 0.0, 0.0));
	}
	long desiredCardinality = long(floor(n / (double)numberOfProcesses));
	// a list containing the partition number each vertex belongs to
	ClusterArray clusterArray(n);
	BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " vertices.";

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

		while(Si.getClusterSize(0) < desiredCardinality) {
			// updates the positive edge sum of each vertex ni in Gr (pos edge sum between ni and Si)
			for(int i = 0; i < Gr.size(); i++) {
				long ni = Gr[i].id;
				Gr[i].positiveDegree = g->getPositiveEdgeSumBetweenVertexAndClustering(ni, Si.getClusterArray());
			}
			std::vector<VertexDegree>::iterator it_max = std::max_element(Gr.begin(), Gr.end(), pos_degree_ordering_asc());  // O(n)
			BOOST_LOG_TRIVIAL(debug) << "ClusterSize = " << Si.getClusterSize(0) << ": The vertex in (Gr - Si) with the biggest pos edge sum is " << it_max->id << " with sum = " << it_max->positiveDegree;

			// Chooses a node ni belonging to (Gr \ Si) with the biggest cardinality of positive edges between ni and Si;
			long ni = it_max->id;
			// Si = Si + ni
			Si.addNodeToCluster(*g, problem, ni, 0, false);
			// Gr = Gr - ni
			Gr.erase(it_max);
		}
		BOOST_LOG_TRIVIAL(debug) << "Partition ready.";
		// Gr = Gr – Si  TODO: subtrair a lista partialGraph da lista residualGraph -> ja feito acima
		// Envia Si para o processo pi+1
		ClusterArray SiArray = Si.getClusterArray();
		for(int i = 0; i < n; i++) {
			if(SiArray[i] >= 0) {
				clusterArray[i] = pi;
			}
		}
	}

	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timer.resume();
	start_time = timer.elapsed();

	// processes the initial solution generated by the split graph partitioning
	BOOST_LOG_TRIVIAL(info) << "Generated a cluster array of " << clusterArray.size() << " vertices.";
	Clustering Cc(clusterArray, *g, problem);
	BOOST_LOG_TRIVIAL(info) << "Initial split graph clustering I(P) = " << Cc.getImbalance().getValue();
	constructivePhaseResults << 1 << "," << Cc.getImbalance().getValue() << ","
					<< Cc.getImbalance().getPositiveValue()
					<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << Cc.getImbalance().getValue()
						<< "\n";
	std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
	for(long i = 0; i < n; i++) {
		long k = clusterArray[i];
		verticesInCluster[k].push_back(i);
	}

	// Creates numberOfProcesses subgraphs
	std::vector<SubGraph> subgraphList;
	// each subgraph will have a subset of the main graph's nodes and edges, based on the previous clustering
	for(int k = 0; k < numberOfProcesses; k++) {
		SubGraph sg = (g->graph).create_subgraph(); //verticesInCluster[k].begin(), verticesInCluster[k].end());
		subgraphList.push_back(sg);  // --> SUBGRAPH COPY CTOR NOT WORKING!!!
		for(std::vector<long>::iterator it = verticesInCluster[k].begin(); it != verticesInCluster[k].end(); it++) {
			add_vertex(*it, subgraphList.back());
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << k;
		}

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << k << ": num_edges = " << num_edges(subgraphList.back()) << " , num_vertices = " << num_vertices(subgraphList.back());
	}

	for(int i = 0; i < numberOfSlaves; i++) {
		// 2. Distribute numberOfSlaves graph parts between the ILS Slave processes
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
							problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd->getTimeLimit(),
							numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k);
		if(k < 0) {
			imsg.setClustering(&CCclustering);
		}
		imsg.setVertexList(verticesInCluster[i+1]);
		world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << slaveList[i];
	}
	// 2.1. the leader does its part of the work: runs ILS using the first part of the divided graph
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
	Clustering leaderClustering = cudails.executeILS(construct2, vnd, &sg, iter, iterMaxILS, perturbationLevelMax, problem, info);

	// 3. the leader receives the processing results
	OutputMessage omsg;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages.";
	// Global cluster array
	ClusterArray globalClusterArray(g->getN(), 0);
	long clusterOffset = 0;
	for(int i = 0; i < numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source() << ". Obj = " <<
				omsg.clustering.getImbalance().getValue();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// 4. Merge the partial solutions into a global solution for the whole graph
		ClusterArray localClusterArray = omsg.clustering.getClusterArray();
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
		long nc = omsg.clustering.getNumberOfClusters();
		clusterOffset += nc;
	}
	// includes the leader's processing result as well
	// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
	std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(sg.graph);
	std::vector<long> globalVertexId;
	for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
		globalVertexId.push_back(sg.graph.local_to_global(*it));
	}
	// 4. Merges the leader's partial solution into a global solution for the whole graph
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
	long nc = leaderClustering.getNumberOfClusters();
	clusterOffset += nc;

	// 5. Builds the clustering with the merge of each process results
	Clustering globalClustering(globalClusterArray, *g, problem);

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Initial solution found: I(P) = " << globalClustering.getImbalance().getValue();
	globalClustering.printClustering(g->getN());

	// 6. Runs a local search over the merged global solution, trying to improve it
	// TODO Run a full ILS iteration over the split graph initial solution
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Applying local search over the initial solution.";
	NeighborhoodSearch* neigborhoodSearch;
	NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, this->numberOfSlaves, this->numberOfSearchSlaves);
	neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	VariableNeighborhoodDescent vnd2(*neigborhoodSearch, vnd->getRandomSeed(), vnd->getNeighborhoodSize(),
			false, vnd->getTimeLimit());
	globalClustering = vnd2.localSearch(g, globalClustering, 0, problem, timeSpentInILS, info.processRank);
	CBest = globalClustering;

	// runs a ILS iteration over the merged global solution
	ConstructClustering NoConstruct(functionFactory.build(construct->getGainFunctionType()), vnd->getRandomSeed(),
			construct->getAlpha(), &globalClustering);
	globalClustering = cudails.executeILS(&NoConstruct, vnd, g, iter, 1, perturbationLevelMax, problem, info);

	// 7. Stops the timer and stores the elapsed time
	timer.stop();
	end_time = timer.elapsed();

	// 8. Write the results into ostream os, using csv format
	// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
	timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
	iterationResults << (1) << "," << globalClustering.getImbalance().getValue() << "," << globalClustering.getImbalance().getPositiveValue()
			<< "," << globalClustering.getImbalance().getNegativeValue() << "," << globalClustering.getNumberOfClusters()
			<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) << "\n";
	Imbalance bestValue = globalClustering.getImbalance();
	int iterationValue = 1;
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
				<< "," << bestValue.getPositiveValue()
				<< "," << bestValue.getNegativeValue()
				<< setprecision(0)
				<< "," << globalClustering.getNumberOfClusters()
				<< "," << (iterationValue+1)
				<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) // timeSpentOnBestSolution
				<< "," << iter
				<< "," << numberOfTestedCombinations << "\n";

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after 1 ILS parallel graph iteration: I(P) = " << globalClustering.getImbalance().getValue();

	// Saves the iteration results to csv file
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);

	return globalClustering;
}

} /* namespace ils */
} /* namespace resolution */
