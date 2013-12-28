/*
 * ParallelGrasp.cpp
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#include "include/ParallelGrasp.h"
#include "../../util/include/MPIMessage.h"
#include <cstring>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

namespace resolution {
namespace grasp {

using namespace util;
namespace mpi = boost::mpi;

ParallelGrasp::ParallelGrasp(GainFunction* f, unsigned long seed) : Grasp(f, seed) {

}

ParallelGrasp::~ParallelGrasp() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelGrasp::executeGRASP(SignedGraph *g, const int& iter,
		const double& alpha, const int& l, const bool& firstImprovementOnOneNeig,
		ClusteringProblem& problem, string& executionId, string& fileId, string& outputFolder,
		const long& timeLimit, const int& numberOfSlaves, const int& myRank,
		const int& numberOfSearchSlaves) {
	mpi::communicator world;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel GRASP] Initiating parallel GRASP...";
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i <= numberOfSlaves; i++) {
		InputMessageParallelGrasp imsg(g->getId(), g->getGraphAsText(), iter, alpha, l,
				problem.getType(), gainFunction->getType(), executionId, fileId, outputFolder, timeLimit,
				numberOfSlaves, numberOfSearchSlaves, firstImprovementOnOneNeig);
		world.send(i * (numberOfSearchSlaves + 1), MPIMessage::INPUT_MSG_PARALLEL_GRASP_TAG, imsg);
		BOOST_LOG_TRIVIAL(trace) << "[Parallel GRASP] Message sent to process " << i;
	}
	// the leader does its part of the work
	ClusteringPtr bestClustering = Grasp::executeGRASP(g, iter, alpha,
			l, firstImprovementOnOneNeig, problem, executionId, fileId,
			outputFolder, timeLimit, numberOfSlaves,
			myRank, numberOfSearchSlaves);

	// the leader receives the processing results
	OutputMessage omsg;
	for(int i = 1; i <= numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
		BOOST_LOG_TRIVIAL(trace) << "[Parallel GRASP] Message received from process " << stat.source() << ": " <<
				omsg.clustering.getImbalance().getValue() << endl << omsg.clustering.toString();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// check if solution value has improved
		if(omsg.clustering.getImbalance().getValue() < bestClustering->getImbalance().getValue()) {
			ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
			bestClustering.reset();
			bestClustering = clustering;
			BOOST_LOG_TRIVIAL(trace) << "*** [Parallel GRASP] Better value found for objective function in node " << stat.source() << ": " <<
					omsg.clustering.getImbalance().getValue();
		}
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel GRASP] Best solution found: I(P) = " << bestClustering->getImbalance().getValue();
	bestClustering->printClustering();
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
