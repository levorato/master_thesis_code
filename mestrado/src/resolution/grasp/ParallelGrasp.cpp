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
		// TODO realizar marshalling da mensagem atraves da boost
		InputMessage imsg; // (g->getGraphAsText(), myIter, alpha, l, problem.getType(), fileId);
		imsg.iter = myIter;
		// processMessage(matrixPtr.get(), &msg, i, partitionSize, TO_MESSAGE);
		MPI_Send(&imsg, sizeof(InputMessage), MPI_BYTE, i, 50, MPI_COMM_WORLD);
		cout << "Message sent to process " << i << endl;
	}
	// o resto da divisao vai pro processo 0
	ClusteringPtr bestClustering = Grasp::executeGRASP(g, myIter + resto, alpha, l, problem, fileId, myRank);

	// receives the processing results
	for(int i = 1; i < np; i++) {
		OutputMessage omsg;
		MPI_Status status;
		// TODO realizar unmarshalling da mensagem pela boost
		MPI_Recv(&omsg, sizeof(OutputMessage), MPI_BYTE, MPI_ANY_SOURCE, OutputMessage::TAG, MPI_COMM_WORLD, &status);
		// TODO colocar numero do processo (obter do MPI)
		cout << "Message received from process " << i << ": " << omsg.clusteringAsText << "\n";
		// reune os resultados da execucao do processo i
		// processMessage(newmatrix, &msg, origem - 1, partitionSize, TO_MATRIX);
		if(clustering->getObjectiveFunctionValue() < bestClustering->getObjectiveFunctionValue()) {
			bestClustering = clustering;
		}
	}
	// retorna o melhor clustering dentre todos os processos
	return bestClustering;
}

} /* namespace grasp */
} /* namespace resolution */
