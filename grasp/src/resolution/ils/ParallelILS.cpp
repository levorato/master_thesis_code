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

#include <cstring>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

namespace resolution {
namespace ils {

using namespace std;
using namespace util;
using namespace util::parallel;
using namespace resolution::construction;
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
		//string s = ProcessUtil::exec("./graclus1.2/graclus --help");
		BOOST_LOG_TRIVIAL(info) << "Invoking Graclus partitioning for spit graph...";
		if(ProcessUtil::exec("./graclus test.g 2")) {
			BOOST_LOG_TRIVIAL(error) << "FATAL: Error invoking Glacus. Please check the logs. Program will now exit!";
			throw std::invalid_argument("FATAL: Error invoking Glacus. Please check the logs.");
		} else {
			BOOST_LOG_TRIVIAL(info) << "Successful. Processing partition...";
		}

		for(int i = 0; i < numberOfSlaves; i++) {
			// 2. Distribute numberOfSlaves graph parts between the ILS Slave processes
			InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
								problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd->getTimeLimit(),
								numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k);
			if(k < 0) {
				imsg.setClustering(&CCclustering);
			}
			world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << slaveList[i];
		}
		// 2.1. the leader does its part of the work: runs ILS using the first part of the divided graph
		CUDAILS cudails;
		Clustering bestClustering = cudails.executeILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info);

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

} /* namespace grasp */
} /* namespace resolution */
