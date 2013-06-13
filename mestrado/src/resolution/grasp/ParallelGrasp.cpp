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

ParallelGrasp::ParallelGrasp() {
	// TODO Auto-generated constructor stub

}

ParallelGrasp::~ParallelGrasp() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelGrasp::executeGRASP(SignedGraph *g, int iter, float alpha, int l,
			ClusteringProblem& problem, string& timestamp, string& fileId, int &np, int &myRank) {
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i < np; i++) {
		InputMessage imsg(g->getGraphAsText(), iter, alpha, l, problem.getType(), fileId);
		world.send(i, INPUT_MSG_TAG, imsg);
		cout << "Message sent to process " << i << endl;
	}
	// the leader does its part of the work
	ClusteringPtr bestClustering = Grasp::executeGRASP(g, iter, alpha,
			l, problem, timestamp, fileId, myRank);

	// the leader receives the processing results
	OutputMessage omsg;
	for(int i = 1; i < np; i++) {
		mpi::status stat = world.recv(mpi::any_source, OUTPUT_MSG_TAG, omsg);
		cout << "Message received from process " << stat.source() << ": " <<
				omsg.clustering.getObjectiveFunctionValue() << endl << omsg.clustering.toString()  << "\n";
		// process the result of the execution of process i
		if(omsg.clustering.getObjectiveFunctionValue() < bestClustering->getObjectiveFunctionValue()) {
			ClusteringPtr clustering = make_shared<Clustering>(omsg.clustering);
			bestClustering.reset();
			bestClustering = clustering;
			cout << "*** Better value found for object function in node " << stat.source() << ": " <<
					omsg.clustering.getObjectiveFunctionValue() << endl;
		}
	}
	cout << "Best solution found: I(P) = " << bestClustering->getObjectiveFunctionValue() << endl;
	bestClustering->printClustering();
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
