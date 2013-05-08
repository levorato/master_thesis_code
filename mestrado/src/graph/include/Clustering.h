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
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;

namespace clusteringgraph {

// Defines the boolean matrix and its pointer
typedef multi_array<bool, 2> BoolMatrix;
typedef shared_ptr<BoolMatrix> ClusterMatrixPtr;
// Defines the neighborhood list (of cluster matrices) and its pointer
typedef vector<ClusterMatrixPtr> NeighborhoodList;
typedef shared_ptr<NeighborhoodList> NeighborhoodListPtr;

class Clustering {
public:
	/**
	 * Creates a Clustering object with n nodes and no clusters.
	 */
	Clustering(int n);
	/**
	 * Creates a Clustering object based on the number of nodes (n)
	 * and numbers of clusters (k).
	 */
	Clustering(int n, int k);
	/**
	 * Creates a Clustering object based on the boolean clustering
	 * matrix (clusterMatrix).
	 */
	Clustering(BoolMatrix* clusterMatrix);
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
	/** number of clusters (k) */
	int numberOfClusters;
	/** the clustering matrix, with dimensions k x n */
	ClusterMatrixPtr clusterMatrixPtr;
	/** the l-neighborhood list of clusters */
	NeighborhoodListPtr neighborhoodListPtr;

	void print(std::ostream& os, BoolMatrix *m);
};

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
