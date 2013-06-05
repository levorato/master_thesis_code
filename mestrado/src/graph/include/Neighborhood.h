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
	ClusteringPtr process2optCombination(SignedGraph& g, Clustering* clustering, int k1, int k2,
			int k3, int k4, int n, int i, int j) {

		// cout << "2-opt-comb: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << ", " << i << ", " << j << endl;
		// clustering->printClustering();
		ClusteringPtr cTemp = make_shared < Clustering > (*clustering, n);
		int nc = cTemp->getNumberOfClusters();

		// 1st step: insertion of nodes i and j into other clusters
		if (k3 >= 0) {
			// inserts i in existing cluster k3
			cTemp->addNodeToCluster(g, i, k3);
			assert(cTemp->getCluster(k3)[i]);
		} else {
			// inserts i in a new cluster (alone)
			cTemp->addCluster(g, i);
		}
		if (k4 >= 0) {
			// inserts j in existing cluster k4
			cTemp->addNodeToCluster(g, j, k4);
		} else {
			// inserts j in a new cluster (alone)
			cTemp->addCluster(g, j);
		}
		// 2nd step: removal of nodes i and j from previous clusters (k1 and k2)
		// attention: in reverse order, to avoid out of range!
		if(k2 > k1) {
			cTemp->removeNodeFromCluster(g, j, k2);
			cTemp->removeNodeFromCluster(g, i, k1);
		} else {
			cTemp->removeNodeFromCluster(g, j, k1);
			cTemp->removeNodeFromCluster(g, i, k2);
		}

		return cTemp;
	}
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
