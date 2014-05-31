/*
 * SequentialNeighborhoodSearch.h
 *
 *  Created on: May 25, 2014
 *      Author: mlevorato
 */

#ifndef SEQNEIGHBORHOODSEARCH_H_
#define SEQNEIGHBORHOODSEARCH_H_

#include <graph/include/NeighborhoodSearch.h>

namespace clusteringgraph {

class SequentialNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	SequentialNeighborhoodSearch();
	virtual ~SequentialNeighborhoodSearch();

	/**
	 * This method is the parallelized CUDA version: neighborhood generation
	 * is done across several cuda cores.
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

	Clustering search1opt(SignedGraph* g,
	                Clustering* clustering, ClusteringProblem& problem,
	                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
	                int myRank, unsigned long initialSearchIndex,
	        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k);

	Clustering search2opt(SignedGraph* g,
		                Clustering* clustering, ClusteringProblem& problem,
		                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
		                int myRank, unsigned long initialSearchIndex,
		        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k);
};

} /* namespace clusteringgraph */

#endif /* SEQNEIGHBORHOODSEARCH_H_ */
