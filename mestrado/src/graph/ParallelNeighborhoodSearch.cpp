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
		offset(offset), numberOfProcesses(numproc) {
	// TODO Auto-generated constructor stub

}

ParallelNeighborhoodSearch::~ParallelNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, unsigned long numberOfSlaves, int myRank,
		unsigned long numberOfSearchSlaves) {

	// Splits the processing in (numberOfClusters / numberOfProcesses) chunks,
	// to be consumed by numberOfProcesses processes
	unsigned long sizeOfChunk = clustering->getNumberOfClusters() / numberOfProcesses;
	unsigned long remainingClusters = clustering->getNumberOfClusters() % numberOfProcesses;
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	int firstSlave = offset + myRank * numberOfProcesses;
	int lastSlave = offset + (myRank + 1) * numberOfProcesses;
	int i = firstSlave;
	for(; i < lastSlave; i++) {
		InputMessageParallelVNS imsgpvns(l, g->getGraphAsText(), *clustering, problem.getType(),
				timeSpentSoFar, timeLimit, i * sizeOfChunk, (i + 1) * sizeOfChunk - 1, numberOfSlaves, numberOfSearchSlaves);
		world.send(i, ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG, imsgpvns);
		cout << "VNS Message sent to process " << i << endl;
	}
	// the leader does its part of the work
	// ClusteringProblemFactory problemFactory;
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
		// process the result of the execution of process i
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
