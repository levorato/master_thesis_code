/*
 * ParallelNeighborhoodSearch.h
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#ifndef PARALLELNEIGHBORHOODSEARCH_H_
#define PARALLELNEIGHBORHOODSEARCH_H_

#include "Neighborhood.h"

namespace resolution {
namespace grasp {

class ParallelNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	ParallelNeighborhoodSearch(unsigned int offset, unsigned int numproc);
	virtual ~ParallelNeighborhoodSearch();

	/**
	 * This method is the parallelized MPI version: neighborhood generation
	 * is done across several processes.
	 */
	virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, const ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					unsigned long numberOfSlaves, int myRank, unsigned long numberOfSearchSlaves);

private:
	/**
	 * Initial process index.
	 */
	unsigned int offset;
	/**
	 * Number of processes that will execute the parallel search.
	 */
	unsigned int numberOfProcesses;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELNEIGHBORHOODSEARCH_H_ */
