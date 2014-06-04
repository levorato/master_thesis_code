/*
 * SequentialNeighborhoodGen.h
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#ifndef SEQUENTIALNEIGHBORHOODGEN_H_
#define SEQUENTIALNEIGHBORHOODGEN_H_

#include "NeighborhoodSearch.h"

namespace clusteringgraph {

class SequentialNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	SequentialNeighborhoodSearch();
	virtual ~SequentialNeighborhoodSearch();

	virtual Clustering searchNeighborhood(int l, SignedGraph* g,
				Clustering* clustering, ClusteringProblem& problem,
				double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
				int myRank, bool firstImprovementOnOneNeig);

private:
	/**
	 * Searches the l-neighborhood of the given clustering, by removing a vertex
	 * whose index is between initialSearchIndex and finalSearchIndex
	 * and adding it to any other cluster available (and also a new cluster). In the case
	 * of 2-opt, the initial and final vertex indices only apply to the first vertex in
	 * the search (not the second one).
	 * The sequential search does first improvement for 1-opt and 2-opt searches.
	 */
	virtual Clustering searchNeighborhood(int l, SignedGraph* g,
					Clustering* clustering, ClusteringProblem& problem,
					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
					int myRank,	unsigned long initialSearchIndex, unsigned long finalSearchIndex,
					bool firstImprovementOnOneNeig, unsigned long k);

};

} /* namespace clusteringgraph */
#endif /* SEQUENTIALNEIGHBORHOODGEN_H_ */
