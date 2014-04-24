/*
 * ParallelILS.h
 *
 *  Created on: Mar 31, 2014
 *      Author: Mario Levorato
 */

#ifndef PARALLELILS_H_
#define PARALLELILS_H_

#include "ILS.h"
#include <mpi.h>

using namespace problem;
using namespace resolution::construction;

namespace resolution {
namespace ils {

class ParallelILS : resolution::ils::ILS {
public:
	ParallelILS();
	virtual ~ParallelILS();

	/**
	 * Triggers the parallel execution of the ILS algorithm using MPI.
	 */
	Clustering executeILS(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
			SignedGraph *g, const int& iter, const int& l, const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const long& timeLimit,
			const int& numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);
};

} /* namespace ils */
} /* namespace resolution */
#endif /* PARALLELILS_H_ */
