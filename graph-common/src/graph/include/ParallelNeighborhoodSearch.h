/*
 * ParallelNeighborhoodSearch.h
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#ifndef PARALLELNEIGHBORHOODSEARCH_H_
#define PARALLELNEIGHBORHOODSEARCH_H_

#include "NeighborhoodSearch.h"

namespace clusteringgraph {

class ParallelNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	ParallelNeighborhoodSearch(int allocStrategy, unsigned int offset, unsigned int numproc);
	virtual ~ParallelNeighborhoodSearch();

	/**
	 * This method is the parallelized MPI version: neighborhood generation
	 * is done across several processes.
	 */
	virtual Clustering searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					int myRank, bool firstImprovementOnOneNeig);

	/**
	 * Searches the l-neighborhood of the given clustering, by removing a vertex
	 * from a cluster (vertex index is between initialSearchIndex and finalSearchIndex)
	 * and adding it to any other cluster available (and also a new cluster). In the case
	 * of 2-opt, the initial and final search indices only apply to the first vertex moved
	 * in the search (not the second one).
	 * This parallel version does best-improvement for 1-opt e first improvement for 2-opt.
	 */
	virtual Clustering searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					int myRank,	unsigned long initialSearchIndex, unsigned long finalSearchIndex,
					bool firstImprovementOnOneNeig, unsigned long k);

private:
	/**
	 * Machine-process allocation strategy for VND slave processes.
	 */
	int machineProcessAllocationStrategy;
	/**
	 * Number of slave processes that will execute the grasp iterations in parallel.
	 */
	unsigned int numberOfSlaves;
	/**
	 * Number of slave processes that will execute the parallel search (parallel VND).
	 */
	unsigned int numberOfSearchSlaves;

};

} /* namespace graph */
#endif /* PARALLELNEIGHBORHOODSEARCH_H_ */
