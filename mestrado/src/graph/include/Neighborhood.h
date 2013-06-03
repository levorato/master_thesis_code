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
	ClusteringPtr process2optCombination(SignedGraph& g, Clustering* clustering, int k1, int k2,
			int k3, int k4, int n, int i, int j) {

		ClusteringPtr cTemp = make_shared < Clustering > (*clustering, n);
		int nc = cTemp->getNumberOfClusters();
		// removes node i from cluster1 and inserts in cluster3
		// TODO check if the removal of node i destroys cluster1
		cTemp->removeNodeFromCluster(g, i, k1);
		// recalculates the number of clusters, as one of them may have been removed
		int newnc1 = cTemp->getNumberOfClusters();
		if (k3 >= 0) {
			// inserts i in existing cluster k3
			if(newnc1 < nc && k3 > k1) {
				// cluster k1 has been removed
				// cluster k3 is in fact (k3 - 1)
				cTemp->addNodeToCluster(g, i, k3 - 1);
			} else {
				cTemp->addNodeToCluster(g, i, k3);
			}
		} else {
			// inserts i in a new cluster (alone)
			cTemp->addCluster(g, i);
		}
		// removes node j from cluster2 and inserts in cluster4
		cTemp->removeNodeFromCluster(g, j, k2);
		if (k4 >= 0) {
			// inserts j in existing cluster k4
			// recalculates the number of clusters, as one or two of them may have been removed
			int newnc2 = cTemp->getNumberOfClusters();
			if((newnc1 < nc) && (k4 > k1)) {
				if((newnc2 < newnc1) && (k4 > k2)) {
					// clusters k1 and k2 have been removed
					// cluster k4 is in fact (k4 - 2)
					cTemp->addNodeToCluster(g, i, k4 - 2);
				} else {
					// only cluster k1 has been removed
					// cluster k4 is in fact (k4 - 1)
					cTemp->addNodeToCluster(g, i, k4 - 1);
				}
			} else {
				if((newnc2 < nc) && (k4 > k2)) {
					// only cluster k2 has been removed
					// cluster k4 is in fact (k4 - 1)
					cTemp->addNodeToCluster(g, i, k4 - 1);
				} else {
					// no cluster has been removed
					// cluster k4 is in fact k4
					cTemp->addNodeToCluster(g, i, k4);
				}
			}
		} else {
			// inserts j in a new cluster (alone)
			cTemp->addCluster(g, j);
		}

		return cTemp;
	}
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
