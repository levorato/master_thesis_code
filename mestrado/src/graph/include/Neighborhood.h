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

#include "Graph.h"
#include "Clustering.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

// Defines the neighborhood list (of cluster lists) and its pointer
class NeighborhoodList : vector<ClusterListPtr> {
public:
	ClusteringPtr candidateClustering;

	NeighborhoodList(Clustering *c, int n);

	/**
	 * Generates a l-neighborhood for this clustering.
	 * @return NeighborhoodList*
	 */
	ClusteringPtr generateNeighborhood(int l, SignedGraph* g, ClusteringProblem* problem);

private:
	int numberOfNodes;
	ClusteringPtr clusteringPtr;
	GraphPtr graphPtr;

	/**
	 * Generates a new cluster by swithcing the node i from cluster k1 to k3, and
	 * node j from cluster k2 to k4. Includes the cases of cluster deletion and
	 * new cluster creation, when the clusters contain only one element.
	 * The parameter n represents the number of nodes in the graph / clustering.
	 * If parameter k3 == -1, inserts node i in a new cluster (alone).
	 * If parameter k4 == -1, inserts node j in a new cluster (anlone).
	 */
	ClusteringPtr process2optCombination(int k1, int k2, int k3, int k4, int n,
			int i, int j);
};
typedef shared_ptr<NeighborhoodList> NeighborhoodListPtr;

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
