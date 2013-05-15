/*
 * Clustering.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <vector>
#include <boost/config.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph.h"

using namespace boost;
using namespace std;

namespace clusteringgraph {

// uses dynamic_bitset for bool array, a high performance and space saving structure
// based on real bits
// the following array is initially empty and needs to be dynamically intialized
typedef boost::dynamic_bitset<> BoolArray;
// Defines the cluster list and its pointer
// the list is made of boolean arrays, indicating that node i is in the cluster
typedef vector<BoolArray> ClusterList;
typedef shared_ptr<ClusterList> ClusterListPtr;
// Defines the neighborhood list (of cluster lists) and its pointer
typedef vector<ClusterListPtr> NeighborhoodList;
typedef shared_ptr<NeighborhoodList> NeighborhoodListPtr;

class Clustering {
public:
	/**
	 * Creates a Clustering object with n nodes and no clusters.
	 */
	Clustering(int n);
	/**
	 * Creates a Clustering object with n nodes based on the clusterList.
	 */
	Clustering(ClusterList* clusterList, int numberOfNodes);
	/**
	 * Creates a Clustering object based on the clusterList.
	 */
	Clustering(ClusterList* clusterList);
	virtual ~Clustering();

	/**
	 * Adds a new cluster to the clustering configuration.
	 */
	void addCluster(int vertexList[], unsigned int arraySize);
	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering();

	/**
	 * Returns the number of nodes in the graph.
	 */
	int getNumberOfNodes();

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	int getNumberOfClusters();

	/**
	 * Generates a l-neighborhood for this clustering.
	 * @return NeighborhoodList*
	 */
	NeighborhoodList* generateNeighborhood(int l);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering *c);

private:
	/** number of nodes in the graph (n) */
	int numberOfNodes;
	/** the cluster list, with dimensions k x n */
	ClusterListPtr clusterListPtr;
	/** the l-neighborhood list of clusters */
	NeighborhoodListPtr neighborhoodListPtr;

	void print(std::ostream& os, ClusterList *l);
};

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
