/*
 * GraclusParallelILS.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#include "./include/GraclusParallelILS.h"
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

GraclusParallelILS::GraclusParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
		const bool& split, const bool& cuda) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split), cudaEnabled(cuda) {

}

GraclusParallelILS::~GraclusParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering GraclusParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {
	/**
	  *
		1- Paralelizar o grafo: supondo n maquinas
			a. definir uma medida de centralidade na rede para os vertices (caminho minimo eh uma)
			b. expandir o grafo a partir de n vertices mais centrais:
					Busca em largura a partir de cada vertice.
					Rodar uma iteracao da busca para cada grafo expandido registrando as intersecoes
					Parar qdo a intersecao chegar a um determinado valor (parametro a definir)
			c. Rodar a melhor versao do CC e RCC sobre os grafos expandidos
			d. Resolver de forma simples as intersecoes e obter uma solucao
			e. Rodar busca local na solucao entrada: idealmente a busca local nao deveria mudar muito
														a solucao encontrada. Somente nas intersecoes.

			http://liuweipingblog.cn/cpp/an-example-of-boost-betweenness-centrality/
	 */

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
	// 1.1. Convert the graph to Graclus input format
	string s = g->convertToGraclusInputFormat();
	// obtains only the filename (without path)
	string filename = g->getGraphFileLocation();
	filename = filename.substr(filename.find_last_of("/") + 1);
	generateGraclusOutputFile(filename, s);

	// 1.2. Invoke Graclus to process the graph
	// Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned int numberOfProcesses = numberOfSlaves + 1;
	string strNP = boost::lexical_cast<std::string>(numberOfProcesses);
	BOOST_LOG_TRIVIAL(info) << "Invoking Graclus partitioning for spit graph (numberOfProcesses = " << strNP << ")...";
	if(ProcessUtil::exec((string("./graclus ") + filename + string(" ") + strNP).c_str())) {
		BOOST_LOG_TRIVIAL(error) << "FATAL: Error invoking Glacus. Please check the logs. Program will now exit!";
		throw std::invalid_argument("FATAL: Error invoking Glacus. Please check the logs.");
	} else {
		BOOST_LOG_TRIVIAL(info) << "Successful. Processing partition...";
	}

	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timer.resume();
	start_time = timer.elapsed();

	// 1.3. Process Graclus output
	std::vector<long> clusterArray = readGraclusResultFile(filename + string(".part.") + strNP);
	BOOST_LOG_TRIVIAL(info) << "Generated a cluster array of " << clusterArray.size() << " vertices.";
	Clustering Cc(clusterArray, *g, problem);
	BOOST_LOG_TRIVIAL(info) << "Initial Graclus clustering I(P) = " << Cc.getImbalance().getValue();
	constructivePhaseResults << 1 << "," << Cc.getImbalance().getValue() << ","
					<< Cc.getImbalance().getPositiveValue()
					<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << Cc.getImbalance().getValue()
						<< "\n";
	std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
	int n = g->getN();
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
		if(not (localClusterArray.size() == omsg.globalVertexId.size())) {
			BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] FAILED to assert: localClusterArray.size() == omsg.globalVertexId.size()";
		}
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
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Applying local search over the initial solution.";
	NeighborhoodSearch* neigborhoodSearch;
	NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, this->numberOfSlaves, this->numberOfSearchSlaves);
	neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	VariableNeighborhoodDescent vnd2(*neigborhoodSearch, vnd->getRandomSeed(), vnd->getNeighborhoodSize(),
			false, vnd->getTimeLimit());
	globalClustering = vnd2.localSearch(g, globalClustering, 0, problem, timeSpentInILS, info.processRank);
	CBest = globalClustering;

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

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after VND: I(P) = " << globalClustering.getImbalance().getValue();

	// Saves the iteration results to csv file
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);

	return globalClustering;
}

void GraclusParallelILS::generateGraclusOutputFile(string filename, string fileContents) {
	namespace fs = boost::filesystem;
	// Generates the output file in the current dir (runtime dir)
	fs::path newFile(filename);
	BOOST_LOG_TRIVIAL(info) << "Creating graclus output file in " << newFile.c_str();
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		BOOST_LOG_TRIVIAL(fatal) << "Can't open graclus output file! Path: " << newFile.c_str();
		throw "Cannot open graclus output file.";
	}
	// Writes file contents to the output file
	os << fileContents;
	os.close();
}

std::vector<long> GraclusParallelILS::readGraclusResultFile(string filename) {
	namespace fs = boost::filesystem;
	fs::path exFile(filename);
	BOOST_LOG_TRIVIAL(info) << "Reading graclus input file in " << exFile.c_str();
	ifstream in;
	in.open(exFile.c_str(), std::ios::in | std::ios::binary);
	if (!in) {
		BOOST_LOG_TRIVIAL(fatal) << "Can't open graclus result file! Path: " << exFile.c_str();
		throw "Cannot open graclus result file.";
	}
	// Reads file contents
	std::ostringstream contents;
	contents << in.rdbuf();
	in.close();

	// parse each line, reading vertexes clusters
	int clusterNumber = 0;
	std::vector<long> clusters;
	stringstream ss(contents.str());
	while(not ss.eof()) {
		string sv;
		ss >> sv;
		trim(sv);
		if(sv.size() == 0)  continue;

		istringstream iss(sv);
		iss >> clusterNumber;
		clusters.push_back(clusterNumber);
	}
	return clusters;
}

} /* namespace ils */
} /* namespace resolution */
