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
		double timeLimit, unsigned long randomSeed, int myRank) {

	// Splits the processing in (numberOfClusters / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = clustering->getNumberOfClusters() / numberOfSearchSlaves;
	unsigned long remainingClusters = clustering->getNumberOfClusters() % numberOfSearchSlaves;
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	// The formulas below determine the first and last search slaves of this grasp slave process
	// IMPORTANT! Process mapping:
	// * 0..numberOfSlaves-1 => p(i) GRASP slave processes (execute parellel GRASP iterations with MPI)
	// * numberOfSlaves..numberOfSearchSlaves-1 => VNS slave processes for p(0) (execute parallel VNS for p(0))
	// * numberOfSlaves+i*numberOfSearchSlaves..numberOfSlaves+(i+1)*numberOfSearchSlaves-1 => VNS slave processes for p(i)
	// Here, p(i) is represented by myRank.
	int i = 0, cont = 0;
	int firstSlave = numberOfSlaves + myRank * numberOfSearchSlaves;
	int lastSlave = numberOfSlaves + (myRank + 1) * numberOfSearchSlaves;
	if(sizeOfChunk > 0) {
		i = firstSlave;
		// Sends the parallel search (VNS) message to the slaves via MPI
		for(cont = 0; i < lastSlave; i++, cont++) {
			InputMessageParallelVNS imsgpvns(l, g->getGraphAsText(), *clustering, problem.getType(),
					timeSpentSoFar, timeLimit, cont * sizeOfChunk, (cont + 1) * sizeOfChunk - 1, numberOfSlaves,
					numberOfSearchSlaves);
			world.send(i, ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG, imsgpvns);
			BOOST_LOG_TRIVIAL(trace) << "VNS Message sent to process " << i << "; [" << cont * sizeOfChunk << ", " << (cont + 1) * sizeOfChunk - 1 << "]" << endl;
		}
	}
	// the leader (me) does its part of the work too
	SequentialNeighborhoodSearch neig;

	BOOST_LOG_TRIVIAL(trace) << "Total number of clusters is " << clustering->getNumberOfClusters() << endl;
	BOOST_LOG_TRIVIAL(trace) << "VNS Parallelization status: " << numberOfSearchSlaves <<
			" search slaves will process " << sizeOfChunk << " clusters each one." << endl;
	BOOST_LOG_TRIVIAL(trace) << "RemainingClusters is " << remainingClusters << endl;

	double bestValue = numeric_limits<double>::infinity();
	ClusteringPtr bestClustering;
	if(remainingClusters > 0) {
		bestClustering = neig.searchNeighborhood(l, g, clustering,
				problem, timeSpentSoFar, timeLimit, randomSeed,
				cont * sizeOfChunk, cont * sizeOfChunk + remainingClusters - 1);
		bestValue = bestClustering->getImbalance().getValue();
	}
	BOOST_LOG_TRIVIAL(debug) << "Waiting for slaves return messages...\n";
	// the leader receives the processing results, if any
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
			}
		}
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel VNS] Best solution found: I(P) = " << bestClustering->getImbalance().getValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}


} /* namespace grasp */
} /* namespace resolution */
