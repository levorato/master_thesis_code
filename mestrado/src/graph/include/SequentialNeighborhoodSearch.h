/*
 * SequentialNeighborhoodGen.h
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#ifndef SEQUENTIALNEIGHBORHOODGEN_H_
#define SEQUENTIALNEIGHBORHOODGEN_H_

#include "Neighborhood.h"

namespace clusteringgraph {

class SequentialNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	SequentialNeighborhoodSearch();
	virtual ~SequentialNeighborhoodSearch();

	virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
				Clustering* clustering, const ClusteringProblem& problem,
				double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
				int myRank);

	/**
	 * Searches the l-neighborhood of the given clustering, by removing a vertex
	 * from a cluster which index is between initialClusterIndex and finalClusterIndex
	 * and adding it to any other cluster available (and also a new cluster). In the case
	 * of 2-opt, the initial e final cluster indices only apply to the first vertex in
	 * the search (not the second one).
	 */
	ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, const ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					unsigned long initialClusterIndex, unsigned long finalClusterIndex);
};

} /* namespace clusteringgraph */
#endif /* SEQUENTIALNEIGHBORHOODGEN_H_ */
