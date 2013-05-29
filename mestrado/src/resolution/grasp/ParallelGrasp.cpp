/*
 * ParallelGrasp.cpp
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#include "include/ParallelGrasp.h"

namespace resolution {
namespace grasp {

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
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i < np; i++) {
		InputMessage imsg(g->getGraphAsText(), iter, alpha, l, problem.getType(), fileId);
		// processMessage(matrixPtr.get(), &msg, i, partitionSize, TO_MESSAGE);
		MPI_Send(&imsg, sizeof(InputMessage), MPI_BYTE, i, InputMessage::TAG, MPI_COMM_WORLD);
		cout << "Message sent to process " << i << endl;
	}
	// o resto da divisao vai pro processo 0
	Grasp::executeGRASP(g, myIter + resto, alpha, l, problem, fileId, myRank);
	// receives the processing results
	for(int i = 1; i < np; i++) {
		OutputMessage omsg;
		MPI_Status status;
		MPI_Recv(&omsg, sizeof(OutputMessage), MPI_BYTE, MPI_ANY_SOURCE, OutputMessage::TAG, MPI_COMM_WORLD, &status);
		cout << "Message received from process " << i << "\n";
		// reune os resultados da execucao do processo i
		// processMessage(newmatrix, &msg, origem - 1, partitionSize, TO_MATRIX);
	}

	// return ;
}

} /* namespace grasp */
} /* namespace resolution */
