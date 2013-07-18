/*
 * SequentialNeighborhoodSearch.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/SequentialNeighborhoodSearch.h"

namespace clusteringgraph {

SequentialNeighborhoodSearch::SequentialNeighborhoodSearch() {

}

ClusteringPtr SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank) {
	unsigned long nc = clustering->getNumberOfClusters();
	return SequentialNeighborhoodSearch::searchNeighborhood(l, g, clustering, problem,
			timeSpentSoFar, timeLimit, randomSeed, myRank, 0, nc - 1);
}

ClusteringPtr SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex) {

	assert(initialClusterIndex < clustering->getNumberOfClusters());
	assert(finalClusterIndex < clustering->getNumberOfClusters());

	if (l == 1) {  // 1-opt
		// Sequential search does first improvement in 1-opt
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, true);
	} else {  // 2-opt
		return this->search2opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, true);
	}
}

} /* namespace clusteringgraph */

namespace clusteringgraph {

SequentialNeighborhoodSearch::~SequentialNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

} /* namespace clusteringgraph */
