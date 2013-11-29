/*
 * Neighborhood.h
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include <vector>
#include <boost/config.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <cassert>

#include "Graph.h"
#include "Clustering.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

// Defines the neighborhood list
class NeighborhoodSearch {
public:
        NeighborhoodSearch() : numberOfTestedCombinations(0) {

        }

        virtual ~NeighborhoodSearch() {

        }

        /**
         * Traverses the l-neighborhood for the given clustering, searching
         * for the best value of objective function.
         * @return NeighborhoodList*
         */
        virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
                        Clustering* clustering, const ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        int myRank, bool firstImprovementOnOneNeig, unsigned long k) = 0;

	long getNumberOfTestedCombinations() {
		return numberOfTestedCombinations;
	}

protected:

    	/**
    	 * Searches the l-neighborhood of the given clustering, by removing a vertex
    	 * from a cluster which index is between initialClusterIndex and finalClusterIndex
    	 * and adding it to any other cluster available (and also a new cluster). In the case
    	 * of 2-opt, the initial and final cluster indices only apply to the first vertex in
    	 * the search (not the second one).
    	 * @param k the maximum number of clusters of RCC Problem (optional).
    	 */
    	virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
    					Clustering* clustering, const ClusteringProblem& problem,
    					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
    					int myRank,	unsigned long initialClusterIndex, unsigned long finalClusterIndex,
    					bool firstImprovementOnOneNeig, unsigned long k) = 0;

        virtual ClusteringPtr search1opt(SignedGraph* g,
                        Clustering* clustering, const ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        int myRank, unsigned long initialClusterIndex,
                		unsigned long finalClusterIndex, bool firstImprovement, unsigned long k);

		virtual ClusteringPtr search2opt(SignedGraph* g,
						Clustering* clustering, const ClusteringProblem& problem,
						double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
						int myRank, unsigned long initialClusterIndex,
						unsigned long finalClusterIndex, bool firstImprovement, unsigned long k);

        /**
         * Generates a new cluster by swithcing the node i from cluster k1 to k3, and
         * node j from cluster k2 to k4. Includes the cases of cluster deletion and
         * new cluster creation, when the clusters contain only one element.
         * The parameter n represents the number of nodes in the graph / clustering.
         * If parameter k3 == -1, inserts node i in a new cluster (alone).
         * If parameter k4 == -1, inserts node j in a new cluster (anlone).
         */
        ClusteringPtr process2optCombination(SignedGraph& g, Clustering* clustering,
        		const ClusteringProblem& problem, int k1, int k2, int k3, int k4,
        		int n, int i, int j);

	/* Number of tested combinations during neighborhood search */
	long numberOfTestedCombinations;
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
