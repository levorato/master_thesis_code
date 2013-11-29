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
	ParallelGrasp(GainFunction* f, unsigned long seed);
	virtual ~ParallelGrasp();

	/**
	 * Triggers the parallel execution of the GRASP algorithm using MPI.
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, const int& iter, const double& alpha,
			const int& l, const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const long& timeLimit,
			const int& numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELGRASP_H_ */
