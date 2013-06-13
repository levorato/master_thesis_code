/*
 * ParallelGrasp.h
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#ifndef PARALLELGRASP_H_
#define PARALLELGRASP_H_

#include "Grasp.h"
#include <mpi.h>

using namespace problem;

namespace resolution {
namespace grasp {

class ParallelGrasp : resolution::grasp::Grasp {
public:
	static const int INPUT_MSG_TAG = 50;
	static const int OUTPUT_MSG_TAG = 60;
	static const int TERMINATE_MSG_TAG = 70;
	static const int LEADER_ID = 0;

	ParallelGrasp();
	virtual ~ParallelGrasp();

	/**
	 * Triggers the parallel execution of the GRASP algorithm using MPI.
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, const int& iter, const double& alpha,
			const int& l, ClusteringProblem& problem, string& timestamp, string& fileId,
			const long& timeLimit, const int& np, const int& myRank);
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELGRASP_H_ */
