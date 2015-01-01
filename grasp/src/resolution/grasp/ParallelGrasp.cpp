/*
 * ParallelGrasp.cpp
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#include "include/ParallelGrasp.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include <cstring>
#include <vector>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

namespace resolution {
namespace grasp {

using namespace std;
using namespace util;
using namespace util::parallel;
using namespace resolution::construction;
namespace mpi = boost::mpi;

/**
 * @param numberOfSlaves number of slaves used for parallel ILS processing
 * @param numberOfSearchSlaves number of slaves used for parallel VND processing
 */
ParallelGrasp::ParallelGrasp(const int& allocationStrategy, const int& slaves, const int& searchSlaves) :
		Grasp(), machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves) {

}

ParallelGrasp::~ParallelGrasp() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelGrasp::executeGRASP(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
		SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info) {
	mpi::communicator world;
	BOOST_LOG_TRIVIAL(info) << "[Parallel GRASP] Initiating parallel GRASP...";
	// max number of clusters (RCC Problem Only)
	long k = 0;
	Clustering* CCclustering = NULL;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = construct.getCCclustering();
		}
	}
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);
	for(int i = 0; i < numberOfSlaves; i++) {
		InputMessageParallelGrasp imsg(g->getId(), g->getGraphFileLocation(), iter, construct.getAlpha(), vnd.getNeighborhoodSize(),
				problem.getType(), construct.getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd.getTimeLimit(),
				numberOfSlaves, numberOfSearchSlaves, vnd.isFirstImprovementOnOneNeig(), k, CCclustering);
		world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_GRASP_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel GRASP] Message sent to process " << slaveList[i];
	}
	// the leader does its part of the work
	Clustering bestClustering = Grasp::executeGRASP(construct, vnd, g, iter, problem, info);

	// the leader receives the processing results
	OutputMessage omsg;
	for(int i = 0; i < numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel GRASP] Message received from process " << stat.source() << ". Obj = " <<
				omsg.clustering.getImbalance().getValue();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// check if solution value has improved
		if(omsg.clustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
			Clustering clustering = omsg.clustering;
			bestClustering = clustering;
			BOOST_LOG_TRIVIAL(trace) << "*** [Parallel GRASP] Better value found for objective function in node " << stat.source() << ": " <<
					omsg.clustering.getImbalance().getValue();
		}
	}
	BOOST_LOG_TRIVIAL(info) << "[Parallel GRASP] Best solution found: I(P) = " << bestClustering.getImbalance().getValue();
	bestClustering.printClustering(g->getN());
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
