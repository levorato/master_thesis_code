/*
 * ParallelILS.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Mario Levorato
 */

#include "include/ParallelILS.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "problem/include/RCCProblem.h"
#include "./include/CUDAILS.h"
#include "util/include/ProcessUtil.h"
#include "../construction/include/GainFunctionFactory.h"
#include "../construction/include/GainFunction.h"

#include <cstring>
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

namespace resolution {
namespace ils {

using namespace std;
using namespace util;
using namespace util::parallel;
using namespace resolution::construction;
using namespace boost::algorithm;
namespace mpi = boost::mpi;


ParallelILS::ParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves, const bool& split) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split) {

}

ParallelILS::~ParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {
	mpi::communicator world;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Initiating parallel ILS...";
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	if(not splitGraph) {  // traditional parallel ILS approach, dividing the number of multistart iterations
		for(int i = 0; i < numberOfSlaves; i++) {
			InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
								problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd->getTimeLimit(),
								numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k);
			if(k < 0) {
				imsg.setClustering(&CCclustering);
			}
			world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Message sent to process " << slaveList[i];
		}
		// the leader does its part of the work
		CUDAILS cudails;
		Clustering bestClustering = cudails.executeILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info);

		// the leader receives the processing results
		OutputMessage omsg;
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Waiting for slaves return messages.";
		for(int i = 0; i < numberOfSlaves; i++) {
			mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Message received from process " << stat.source() << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			// process the result of the execution of process i
			// sums the total number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// check if solution value has improved
			if(omsg.clustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
				Clustering clustering = omsg.clustering;
				bestClustering = clustering;
				BOOST_LOG_TRIVIAL(trace) << "*** [Parallel ILS] Better value found for objective function in node " << stat.source() << ": " <<
						omsg.clustering.getImbalance().getValue();
			}
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Best solution found: I(P) = " << bestClustering.getImbalance().getValue();
		bestClustering.printClustering(g->getN());
		return bestClustering;
	} else {
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

		// 1. Divide the graph into (numberOfSlaves + 1) non-overlapping parts
		// 1.1. Convert the graph to Graclus input format
		string s = g->convertToGraclusInputFormat();
		// obtains only the filename (without path)
		string filename = g->getGraphFileLocation();
		filename = filename.substr(filename.find_last_of("/") + 1);
		generateGraclusOutputFile(filename, s);

		// 1.2. Invoke Graclus to process the graph
		//string s = ProcessUtil::exec("./graclus1.2/graclus --help");
		unsigned int numberOfProcesses = numberOfSlaves + 1;
		string strNP = boost::lexical_cast<std::string>(numberOfProcesses);
		BOOST_LOG_TRIVIAL(info) << "Invoking Graclus partitioning for spit graph (numberOfProcesses = " << strNP << ")...";
		if(ProcessUtil::exec((string("./graclus ") + filename + string(" ") + strNP).c_str())) {
			BOOST_LOG_TRIVIAL(error) << "FATAL: Error invoking Glacus. Please check the logs. Program will now exit!";
			throw std::invalid_argument("FATAL: Error invoking Glacus. Please check the logs.");
		} else {
			BOOST_LOG_TRIVIAL(info) << "Successful. Processing partition...";
		}

		// 1.3. Process Graclus output
		std::vector<long> clusterArray = readGraclusResultFile(filename + string(".part.") + strNP);
		BOOST_LOG_TRIVIAL(info) << "Generated a cluster array of " << clusterArray.size() << " vertices.";
		Clustering c(clusterArray, *g, problem);
		BOOST_LOG_TRIVIAL(info) << "Initial Graclus clustering I(P) = " << c.getImbalance().getValue();
		std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
		int n = g->getN();
		for(long i = 0; i < n; i++) {
			long k = clusterArray[i];
			verticesInCluster[k].push_back(i);
		}

		std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(g->graph);
		// Creates numberOfProcesses subgraphs
		std::vector<SubGraph> subgraphList;
		// each subgraph will have a subset of the main graph's nodes and edges, based on the previous clustering
		for(int k = 0; k < numberOfProcesses; k++) {
			// SubGraph sg = (g->graph).create_subgraph(verticesInCluster[k].begin(), verticesInCluster[k].end());
			// SubGraph sg = (g->graph).create_subgraph(v_it.first, v_it.second);
			SubGraph sg = (g->graph).create_subgraph();
			for(std::vector<long>::iterator it = verticesInCluster[k].begin(); it != verticesInCluster[k].end(); it++) {
				add_vertex(*it, sg);
				// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << k;
			}

			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] vertices_num = " << verticesInCluster[k].size();
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 1. num_edges = " << num_edges(sg) << " , num_vertices = " << num_vertices(sg);
			subgraphList.push_back(sg);

			std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(sg);
			stringstream ss("Vertices in subgraph: ");
			for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
				ss << *it << " (" << sg.local_to_global(*it) << "), ";
			}
			BOOST_LOG_TRIVIAL(info) << ss.str();
			stringstream ss2("Edges in subgraph: ");
			boost::property_map<SubGraph, edge_properties_t>::type ew = boost::get(edge_properties, sg);
			for(int i = 0; i < num_vertices(sg); i++) {
				// iterates over out-edges of vertex i
				//typedef graph_traits<Graph> GraphTraits;
				DirectedGraph::out_edge_iterator out_i, out_end;
				DirectedGraph::edge_descriptor e;

				// std::cout << "out-edges of " << i << ": ";
				for (tie(out_i, out_end) = out_edges(i, sg); out_i != out_end; ++out_i) {
					e = *out_i;
					Vertex src = source(e, sg), targ = target(e, sg);
					double weight = ew[e].weight;
					ss2 << "*(" << src.id << "," << targ.id << ") = " << weight  << ", ";
				}
			}
			BOOST_LOG_TRIVIAL(info) << "Edges: " << ss2.str();
		}

		for(int i = 0; i < numberOfSlaves; i++) {
			// SubGraph sg = subgraphList[i];
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] 2. num_edges = " << num_edges(sg) << " , num_vertices = " << num_vertices(sg);

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
		SignedGraph sg(verticesInCluster[0].size());
		sg.graph = (g->graph).create_subgraph();

		for(std::vector<long>::iterator it = verticesInCluster[0].begin(); it != verticesInCluster[0].end(); it++) {
			add_vertex(*it, sg.graph);
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << k;
		}
		BOOST_LOG_TRIVIAL(info) << "n =  " << num_vertices(sg.graph);
		BOOST_LOG_TRIVIAL(info) << "e =  " << num_edges(sg.graph);

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

		Clustering bestClustering = cudails.executeILS(construct2, vnd, &sg, iter, iterMaxILS, perturbationLevelMax, problem, info);

		// 3. the leader receives the processing results
		OutputMessage omsg;
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Waiting for slaves return messages.";
		for(int i = 0; i < numberOfSlaves; i++) {
			mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source() << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			// process the result of the execution of process i
			// sums the total number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// 4. Merge the partial solutions into a global solution for the whole graph

			/*
			if(omsg.clustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
				Clustering clustering = omsg.clustering;
				bestClustering = clustering;
				BOOST_LOG_TRIVIAL(trace) << "*** [Parallel ILS] Better value found for objective function in node " << stat.source() << ": "  <<
						omsg.clustering.getImbalance().getValue();
			} */
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Best solution found: I(P) = " << bestClustering.getImbalance().getValue();
		bestClustering.printClustering(g->getN());

		// 5. Run a local search over the merged global solution, trying to improve it


		BOOST_LOG_TRIVIAL(error) << "[Parallel ILS SplitGraph] TODO: ParallelILS with split graph.";
		cout << "[Parallel ILS SplitGraph] TODO: ParallelILS with split graph.\n";
		return Clustering();
	}
}

void ParallelILS::generateGraclusOutputFile(string filename, string fileContents) {
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

std::vector<long> ParallelILS::readGraclusResultFile(string filename) {
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

} /* namespace grasp */
} /* namespace resolution */
