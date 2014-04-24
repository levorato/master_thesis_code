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
	/**
	 * @param numberOfSlaves number of slaves used for parallel ILS processing
	 * @param numberOfSearchSlaves number of slaves used for parallel VND processing
	 */
	ParallelILS(const int& slaves, const int& searchSlaves);
	virtual ~ParallelILS();

	/**
	 * Triggers the parallel execution of the ILS algorithm using MPI.
	 */
	Clustering executeILS(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
			SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info);

private:
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
};

} /* namespace ils */
} /* namespace resolution */
#endif /* PARALLELILS_H_ */
