/*
 * ParallelNeighborhoodSearch.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#include "include/ParallelNeighborhoodSearch.h"
#include "include/SequentialNeighborhoodSearch.h"
#include "../util/include/MPIMessage.h"
#include "../resolution/grasp/include/ParallelGrasp.h"
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

#include <limits>

namespace resolution {
namespace grasp {

using namespace boost::mpi;
using namespace util;
using namespace std;

ParallelNeighborhoodSearch::ParallelNeighborhoodSearch(unsigned int offset, unsigned int numproc) :
		numberOfSlaves(offset), numberOfSearchSlaves(numproc) {

}

ParallelNeighborhoodSearch::~ParallelNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, bool firstImprovementOnOneNeig) {

	// Splits the processing in (numberOfClusters / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = clustering->getNumberOfClusters() / numberOfSearchSlaves;
	unsigned long remainingClusters = clustering->getNumberOfClusters() % numberOfSearchSlaves;
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	// The formulas below determine the first and last search slaves of this grasp slave process
	// IMPORTANT! Process mapping:
	// * 1..numberOfSlaves => p(i) GRASP slave processes (execute parellel GRASP iterations with MPI)
	// * numberOfSlaves+1..numberOfSearchSlaves => VNS slave processes for p(0) (execute parallel VNS for p(0))
	// * numberOfSlaves+1+i*numberOfSearchSlaves..numberOfSlaves+(i+1)*numberOfSearchSlaves => VNS slave processes for p(i)
	// Here, p(i) is represented by myRank.
	int i = 0, cont = 0;
	int firstSlave = numberOfSlaves + 1 + myRank * numberOfSearchSlaves;
	int lastSlave = numberOfSlaves + 1 + (myRank + 1) * numberOfSearchSlaves;
	if(sizeOfChunk > 0) {
		i = firstSlave;
		// Sends the parallel search (VNS) message to the slaves via MPI
		for(cont = 0; i < lastSlave; i++, cont++) {
			InputMessageParallelVNS imsgpvns(g->getId(), l, g->getGraphAsText(), *clustering, problem.getType(),
					timeSpentSoFar, timeLimit, cont * sizeOfChunk, (cont + 1) * sizeOfChunk - 1, numberOfSlaves,
					numberOfSearchSlaves);
			world.send(i, ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG, imsgpvns);
			BOOST_LOG_TRIVIAL(trace) << "VNS Message sent to process " << i << "; [" << cont * sizeOfChunk
					<< ", " << (cont + 1) * sizeOfChunk - 1 << "]" << endl;
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
				cont * sizeOfChunk, cont * sizeOfChunk + remainingClusters - 1, false);
		bestValue = bestClustering->getImbalance().getValue();
	}
	BOOST_LOG_TRIVIAL(debug) << "Waiting for slaves return messages...\n";
	// the leader receives the processing results, if that is the case
	// TODO implement global first improvement (l = 2-opt) in parallel VNS here!
	if(sizeOfChunk > 0) {
		OutputMessage omsg;
		for(i = firstSlave; i < lastSlave; i++) {
			mpi::status stat = world.recv(mpi::any_source, ParallelGrasp::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
			BOOST_LOG_TRIVIAL(debug) << "Message received from process " << stat.source() << ": " <<
					omsg.clustering.getImbalance().getValue() << endl << omsg.clustering.toString()  << "\n";
			// processes the result of the execution of process p(i)
			if(omsg.clustering.getImbalance().getValue() < bestValue) {
				ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
				bestClustering = clustering;
				bestValue = omsg.clustering.getImbalance().getValue();
				BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VNS] Better value found for objective function in node "
						<< stat.source() << ": " <<
						omsg.clustering.getImbalance().getValue() << endl;
				// terminate other slaves' VNS search
				for(int pj = firstSlave; pj < lastSlave; pj++) {
					world.send(pj, ParallelGrasp::INTERRUPT_MSG_PARALLEL_VNS_TAG);
				}
			}
		}
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel VNS] Best solution found: I(P) = " << bestClustering->getImbalance().getValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}

ClusteringPtr ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex, bool firstImprovementOnOneNeig) {

	assert(initialClusterIndex < clustering->getNumberOfClusters());
	assert(finalClusterIndex < clustering->getNumberOfClusters());

	if (l == 1) {  // 1-opt
		// Parallel search always does best improvement in 1-opt
		// Therefore, parameter firstImprovementOnOneNeig is ignored
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, false);
	} else {  // 2-opt
		// TODO implement global first improvement in parallel 2-opt
		// first improvement found in one process must break all other processes' loop
		return this->search2opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, true);
	}
}


} /* namespace grasp */
} /* namespace resolution */
