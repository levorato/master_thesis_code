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
			ClusteringProblem& problem, std::ostream& os, int &np) {
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	for(int i = 1; i < np; i++) {
		// populateMessage(matrixPtr.get(), &msg, j - 1, partitionSize, TO_MESSAGE);
		// MPI_Isend(&msg, sizeof(Message), MPI_BYTE, j, tag, MPI_COMM_WORLD, &flag[j]);
	}
	Grasp::executeGRASP(g, iter / np, alpha, l, problem, os);
	// receives the processing results
	for(int i = 1; i < np; i++) {
		// MPI_Recv(&msg, sizeof(Message), MPI_BYTE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
		// reune os resultados da execucao do processo i
		// populateMessage(newmatrix, &msg, origem - 1, partitionSize, TO_MATRIX);
	}

	// return ;
}

} /* namespace grasp */
} /* namespace resolution */
