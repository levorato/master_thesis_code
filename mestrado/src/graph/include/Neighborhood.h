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

#include "Graph.h"
#include "Clustering.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

// Defines the neighborhood list
class NeighborhoodListGenerator {
public:
	NeighborhoodListGenerator(int n) : numberOfNodes(n) {

	}

	/**
	 * Generates a l-neighborhood for this clustering.
	 * @return NeighborhoodList*
	 */
	virtual ClusteringPtr generateNeighborhood(int l, SignedGraph* g,
			Clustering* clustering, const ClusteringProblem& problem) = 0;

protected:
	int numberOfNodes;

	/**
	 * Generates a new cluster by swithcing the node i from cluster k1 to k3, and
	 * node j from cluster k2 to k4. Includes the cases of cluster deletion and
	 * new cluster creation, when the clusters contain only one element.
	 * The parameter n represents the number of nodes in the graph / clustering.
	 * If parameter k3 == -1, inserts node i in a new cluster (alone).
	 * If parameter k4 == -1, inserts node j in a new cluster (anlone).
	 */
	ClusteringPtr process2optCombination(Clustering* clustering, int k1, int k2,
			int k3, int k4, int n, int i, int j) {

		ClusteringPtr cTemp = make_shared < Clustering > (*clustering, n);
		// removes node i from cluster1 and inserts in cluster3
		// TODO check if the removal of node i destroys cluster1
		cTemp->removeNodeFromCluster(i, k1);
		if (k3 >= 0) {
			// inserts i in existing cluster k3
			cTemp->addNodeToCluster(i, k3);
		} else {
			// inserts i in a new cluster (alone)
			int nodeArray[1] = { i };
			cTemp->addCluster(nodeArray, 1);
		}
		// removes node j from cluster2 and inserts in cluster4
		cTemp->removeNodeFromCluster(j, k2);
		if (k4 >= 0) {
			// inserts j in existing cluster k4
			cTemp->addNodeToCluster(j, k4);
		} else {
			// inserts j in a new cluster (alone)
			int nodeArray[1] = { j };
			cTemp->addCluster(nodeArray, 1);
		}

		return cTemp;
	}
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
