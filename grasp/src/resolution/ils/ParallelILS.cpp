/*
 * ParallelILS.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Mario Levorato
 */

#include "include/ParallelILS.h"
#include "../../util/include/MPIMessage.h"
#include "../../util/parallel/include/MPIUtil.h"
#include <cstring>
#include <vector>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

namespace resolution {
namespace ils {

using namespace std;
using namespace util;
using namespace util::parallel;
using namespace resolution::construction;
namespace mpi = boost::mpi;

ParallelILS::ParallelILS(const int& slaves, const int& searchSlaves) : ILS(),
		numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves) {

}

ParallelILS::~ParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelILS::executeILS(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
		SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info) {
	mpi::communicator world;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Initiating parallel ILS...";
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);
	for(int i = 0; i < numberOfSlaves; i++) {
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct.getAlpha(), vnd.getNeighborhoodSize(),
				problem.getType(), construct.getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd.getTimeLimit(),
				numberOfSlaves, numberOfSearchSlaves, vnd.isFirstImprovementOnOneNeig());
		world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS] Message sent to process " << slaveList[i];
	}
	// the leader does its part of the work
	Clustering bestClustering = ILS::executeILS(construct, vnd, g, iter, problem, info);

	// the leader receives the processing results
	OutputMessage omsg;
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
	bestClustering.printClustering();
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
