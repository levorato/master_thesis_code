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

namespace resolution {
namespace grasp {

using namespace boost::mpi;
using namespace util;

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
	int firstSlave = numberOfSlaves + myRank * numberOfSearchSlaves;
	int lastSlave = numberOfSlaves + (myRank + 1) * numberOfSearchSlaves;
	int i = firstSlave;
	// Sends the parallel search (VNS) message to the slaves via MPI
	for(; i < lastSlave; i++) {
		InputMessageParallelVNS imsgpvns(l, g->getGraphAsText(), *clustering, problem.getType(),
				timeSpentSoFar, timeLimit, i * sizeOfChunk, (i + 1) * sizeOfChunk - 1, numberOfSlaves,
				numberOfSearchSlaves);
		world.send(i, ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG, imsgpvns);
		cout << "VNS Message sent to process " << i << endl;
	}
	// the leader (me) does its part of the work too
	SequentialNeighborhoodSearch neig;
	ClusteringPtr bestClustering = neig.searchNeighborhood(l, g, clustering,
			problem, timeSpentSoFar, timeLimit, randomSeed,
			i * sizeOfChunk, i * sizeOfChunk + remainingClusters);
	// the leader receives the processing results
	OutputMessage omsg;
	for(int i = firstSlave; i < lastSlave; i++) {
		mpi::status stat = world.recv(mpi::any_source, ParallelGrasp::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
		cout << "Message received from process " << stat.source() << ": " <<
				omsg.clustering.getImbalance().getValue() << endl << omsg.clustering.toString()  << "\n";
		// processes the result of the execution of process p(i)
		if(omsg.clustering.getImbalance().getValue() < bestClustering->getImbalance().getValue()) {
			ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
			bestClustering.reset();
			bestClustering = clustering;
			cout << "*** [Parallel VNS] Better value found for objective function in node " << stat.source() << ": " <<
					omsg.clustering.getImbalance().getValue() << endl;
		}
	}
	cout << "[Parallel VNS] Best solution found: I(P) = " << bestClustering->getImbalance().getValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}


} /* namespace grasp */
} /* namespace resolution */
