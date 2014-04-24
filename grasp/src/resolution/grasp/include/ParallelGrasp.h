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
using namespace resolution::construction;

namespace resolution {
namespace grasp {

class ParallelGrasp : resolution::grasp::Grasp {
public:
	ParallelGrasp(const int& slaves, const int& searchSlaves);
	virtual ~ParallelGrasp();

	/**
	 * Triggers the parallel execution of the GRASP algorithm using MPI.
	 */
	Clustering executeGRASP(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
			SignedGraph *g, const int& iter, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const int& myRank);

private:
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELGRASP_H_ */
