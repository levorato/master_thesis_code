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
		const double& alpha, const int& l, ClusteringProblem& problem,
		string& timestamp, string& fileId, string& outputFolder, 
		const long& timeLimit, const int& numberOfSlaves,
		const int& myRank, const int& numberOfSearchSlaves) {
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i < numberOfSlaves; i++) {
		InputMessageParallelGrasp imsg(g->getId(), g->getGraphAsText(), iter, alpha, l,
				problem.getType(), gainFunction->getType(), fileId, outputFolder, timeLimit,
				numberOfSlaves, numberOfSearchSlaves);
		world.send(i, INPUT_MSG_PARALLEL_GRASP_TAG, imsg);
		cout << "Message sent to process " << i << endl;
	}
	// the leader does its part of the work
	ClusteringPtr bestClustering = Grasp::executeGRASP(g, iter, alpha,
			l, problem, timestamp, fileId, outputFolder, timeLimit, numberOfSlaves,
			myRank, numberOfSearchSlaves);

	// the leader receives the processing results
	OutputMessage omsg;
	for(int i = 1; i < numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
		cout << "Message received from process " << stat.source() << ": " <<
				omsg.clustering.getImbalance().getValue() << endl << omsg.clustering.toString()  << "\n";
		// process the result of the execution of process i
		if(omsg.clustering.getImbalance().getValue() < bestClustering->getImbalance().getValue()) {
			ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
			bestClustering.reset();
			bestClustering = clustering;
			cout << "*** [Parallel GRASP] Better value found for objective function in node " << stat.source() << ": " <<
					omsg.clustering.getImbalance().getValue() << endl;
		}
	}
	cout << "Best solution found: I(P) = " << bestClustering->getImbalance().getValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
