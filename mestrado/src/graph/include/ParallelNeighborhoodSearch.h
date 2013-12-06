/*
 * ParallelNeighborhoodSearch.h
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#ifndef PARALLELNEIGHBORHOODSEARCH_H_
#define PARALLELNEIGHBORHOODSEARCH_H_

#include "NeighborhoodSearch.h"

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
					Clustering* clustering, ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					int myRank, bool firstImprovementOnOneNeig, unsigned long k);

	/**
	 * Searches the l-neighborhood of the given clustering, by removing a vertex
	 * from a cluster which index is between initialClusterIndex and finalClusterIndex
	 * and adding it to any other cluster available (and also a new cluster). In the case
	 * of 2-opt, the initial and final cluster indices only apply to the first vertex in
	 * the search (not the second one).
	 * This parallel version does best-improvement for 1-opt e first improvement for 2-opt.
	 */
	virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					int myRank,	unsigned long initialClusterIndex, unsigned long finalClusterIndex,
					bool firstImprovementOnOneNeig, unsigned long k);

private:
	/**
	 * Number of slave processes that will execute the grasp iterations in parallel.
	 */
	unsigned int numberOfSlaves;
	/**
	 * Number of slave processes that will execute the parallel search (parallel VNS).
	 */
	unsigned int numberOfSearchSlaves;

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELNEIGHBORHOODSEARCH_H_ */
