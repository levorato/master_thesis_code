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

// TODO complete the MPI code
ClusteringPtr ParallelGrasp::executeGRASP(SignedGraph *g, int iter, float alpha, int l,
			ClusteringProblem& problem, string& fileId, int &np, int &myRank) {
	int myIter = iter / np;
	int resto = iter - (myIter * np);
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i < np; i++) {
		InputMessage imsg(g->getGraphAsText(), myIter, alpha, l, problem.getType(), fileId);
		world.send(i, INPUT_MSG_TAG, imsg);
		cout << "Message sent to process " << i << endl;
	}
	// o resto da divisao vai pro processo 0
	ClusteringPtr bestClustering = Grasp::executeGRASP(g, myIter + resto, alpha,
			l, &problem, fileId, myRank);

	// receives the processing results
	for(int i = 1; i < np; i++) {
		OutputMessage omsg;
		mpi::status stat = world.recv(mpi::any_source, OUTPUT_MSG_TAG, omsg);
		cout << "Message received from process " << stat.source() << ": " << omsg.clusteringAsText << "\n";
		// processa o resultado da execucao do processo i
		if(omsg.objectiveFunctionValue < bestClustering->getObjectiveFunctionValue()) {
			// bestClustering.reset();
			// bestClustering = clustering;
			cout << "*** Encontrado melhor valor para a FO: " << omsg.objectiveFunctionValue << endl;
		}

	}
	// retorna o melhor clustering dentre todos os processos
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
