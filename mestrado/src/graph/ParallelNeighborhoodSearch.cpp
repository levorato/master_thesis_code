/*
 * ParallelNeighborhoodSearch.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#include "include/ParallelNeighborhoodSearch.h"
#include "include/SequentialNeighborhoodSearch.h"
#include "../util/include/MPIMessage.h"
#include "../util/parallel/include/MPIUtil.h"
#include "../resolution/grasp/include/ParallelGrasp.h"
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

#include <limits>

namespace resolution {
namespace grasp {

using namespace boost::mpi;
using namespace util;
using namespace util::parallel;
using namespace std;

ParallelNeighborhoodSearch::ParallelNeighborhoodSearch(unsigned int offset, unsigned int numproc) :
		numberOfSlaves(offset), numberOfSearchSlaves(numproc) {

}

ParallelNeighborhoodSearch::~ParallelNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, bool firstImprovementOnOneNeig,
		unsigned long k) {

	// Resets the number of combinations tested on neighborhood search
	numberOfTestedCombinations = 0;

	// Splits the processing in (numberOfClusters / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = clustering->getNumberOfClusters() / numberOfSearchSlaves;
	unsigned long remainingClusters = clustering->getNumberOfClusters() % numberOfSearchSlaves;
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	int i = 0;
	std::vector<int> slaveList;
	if(sizeOfChunk > 0) {
		MPIUtil::populateListOfVNSSlaves(slaveList, myRank, numberOfSlaves, numberOfSearchSlaves);
		// Sends the parallel search (VNS) message to the search slaves via MPI
		for(i = 0; i < numberOfSearchSlaves; i++) {
			InputMessageParallelVNS imsgpvns(g->getId(), l, g->getGraphFileLocation(), *clustering, problem.getType(),
					timeSpentSoFar, timeLimit, i * sizeOfChunk, (i + 1) * sizeOfChunk - 1, numberOfSlaves,
					numberOfSearchSlaves, k);
			world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_VNS_TAG, imsgpvns);
			BOOST_LOG_TRIVIAL(trace) << "VNS Message sent to process " << slaveList[i] << "; [" << i * sizeOfChunk
					<< ", " << (i + 1) * sizeOfChunk - 1 << "]" << endl;
		}
	}
	// the leader (me) does its part of the work too
	BOOST_LOG_TRIVIAL(trace) << "Total number of clusters is " << clustering->getNumberOfClusters() << endl;
	BOOST_LOG_TRIVIAL(trace) << "VNS Parallelization status: " << numberOfSearchSlaves <<
			" search slaves will process " << sizeOfChunk << " clusters each one." << endl;
	BOOST_LOG_TRIVIAL(trace) << "RemainingClusters is " << remainingClusters << endl;

	double bestValue = numeric_limits<double>::infinity();
	ClusteringPtr bestClustering;
	if(remainingClusters > 0) {
		bestClustering = this->searchNeighborhood(l, g, clustering,
				problem, timeSpentSoFar, timeLimit, randomSeed, myRank,
				i * sizeOfChunk, i * sizeOfChunk + remainingClusters - 1, false, k);
		bestValue = bestClustering->getImbalance().getValue();
	}
	BOOST_LOG_TRIVIAL(debug) << "Waiting for slaves return messages...\n";
	// the leader receives the processing results, if that is the case
	// TODO implement global first improvement (l = 2-opt) in parallel VNS here!
	if(sizeOfChunk > 0) {
		OutputMessage omsg;
		for(i = 0; i < numberOfSearchSlaves; i++) {
			mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
			BOOST_LOG_TRIVIAL(debug) << "Message received from process " << stat.source() << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			BOOST_LOG_TRIVIAL(trace) << omsg.clustering.toString();
			// processes the result of the execution of process p(i)
			// sums the number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// checks if the obj value improved
			if(omsg.clustering.getImbalance().getValue() < bestValue) {
				ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
				bestClustering = clustering;
				bestValue = omsg.clustering.getImbalance().getValue();
				BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VNS] Better value found for objective function in node "
						<< stat.source() << ": " <<
						omsg.clustering.getImbalance().getValue() << endl;
				// IMPORTANT: Terminate other slaves' VNS search if 2-opt and CC Problem
				if( (l == 2) and (problem.getType() == ClusteringProblem::CC_PROBLEM) ) {
					InputMessageParallelVNS imsg;
					BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VNS] First improvement on 2-opt: interrupting other VNS slaves.";
					for(int j = 0; j < numberOfSearchSlaves; j++) {
						world.send(slaveList[j], MPIMessage::INTERRUPT_MSG_PARALLEL_VNS_TAG, imsg);
					}
				}
			}
		}
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel VNS] Best solution found: Obj = " << bestClustering->getImbalance().getValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}

ClusteringPtr ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex, bool firstImprovementOnOneNeig, unsigned long k) {

	assert(initialClusterIndex < clustering->getNumberOfClusters());
	assert(finalClusterIndex < clustering->getNumberOfClusters());

	if (l == 1) {  // 1-opt
		// Parallel search always does best improvement in 1-opt
		// Therefore, parameter firstImprovementOnOneNeig is ignored
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, false, k);
	} else {  // 2-opt
		// TODO implement global first improvement in parallel 2-opt
		// first improvement found in one process must break all other processes' loop
		// IMPORTANT: Parallel VNS does first improvement on 2-opt if CC Problem
		// Or best improvement on 2-opt if RCC Problem
		bool firstImprovementOn2Opt = true;
		if(problem.getType() == ClusteringProblem::CC_PROBLEM) {
			firstImprovementOn2Opt = true;
		} else if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			firstImprovementOn2Opt = false;
		}
		return this->search2opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, firstImprovementOn2Opt, k);
	}
}


} /* namespace grasp */
} /* namespace resolution */
