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
	// Parallel Grasp messa ta
	static const int INPUT_MSG_PARALLEL_GRASP_TAG = 50;
	static const int OUTPUT_MSG_PARALLEL_GRASP_TAG = 60;
	// Parallel VNS messa
	static const int INPUT_MSG_PARALLEL_VNS_TAG = 70;
	static const int OUTPUT_MSG_PARALLEL_VNS_TAG = 80;
	static const int TERMINATE_MSG_TAG = 90;
	static const int LEADER_ID = 0;

	ParallelGrasp(GainFunction* f, unsigned long seed);
	virtual ~ParallelGrasp();

	/**
	 * Triggers the parallel execution of the GRASP algorithm using MPI.
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, const int& iter, const double& alpha,
			const int& l, ClusteringProblem& problem, string& timestamp, string& fileId,
			string& outputFolder, const long& timeLimit, const int& np, const int& myRank);
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELGRASP_H_ */
