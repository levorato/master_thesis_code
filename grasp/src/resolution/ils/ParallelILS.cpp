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

#include "../construction/include/GainFunctionFactory.h"
#include "../construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "./splitgraph/include/GraclusParallelILS.h"
#include "splitgraph/include/ImbalanceSubgraphParallelILS.h"

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

using namespace std;
using namespace util;
using namespace util::parallel;
using namespace resolution::construction;
using namespace boost::algorithm;
namespace mpi = boost::mpi;


ParallelILS::ParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
		const bool& split, const bool& cuda) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split), cudaEnabled(cuda) {

}

ParallelILS::ParallelILS(ParallelILS &parILS) : machineProcessAllocationStrategy(parILS.machineProcessAllocationStrategy),
		numberOfSlaves(parILS.numberOfSlaves), numberOfSearchSlaves(parILS.numberOfSearchSlaves),
		splitGraph(parILS.splitGraph), cudaEnabled(parILS.cudaEnabled) {

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
								numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k, cudaEnabled);
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
		// GraclusParallelILS gpils(machineProcessAllocationStrategy, numberOfSlaves, numberOfSearchSlaves, splitGraph, cudaEnabled);
		ImbalanceSubgraphParallelILS gpils(machineProcessAllocationStrategy, numberOfSlaves, numberOfSearchSlaves, splitGraph, cudaEnabled);
		Clustering Gc = gpils.executeILS(construct, vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info);
		return Gc;
	}
}



} /* namespace grasp */
} /* namespace resolution */
