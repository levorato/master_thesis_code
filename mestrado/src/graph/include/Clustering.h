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
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;

namespace clusteringgraph {

// Defines the boolean matrix and its pointer
typedef multi_array<bool, 2> BoolMatrix;
typedef scoped_ptr<BoolMatrix> ClusterMatrixPtr;
// Defines the neighborhood list (of cluster matrices) and its pointer
typedef vector<ClusterMatrixPtr> NeighborhoodList;
typedef shared_ptr<NeighborhoodList> NeighborhoodListPtr;

class Clustering {
public:
	/**
	 * Creates a Clustering object based on the number of nodes (n)
	 * and numbers of clusters (k).
	 */
	Clustering(int n, int k);
	/**
	 * Creates a Clustering object based on the boolean clustering
	 * matrix (clusterMatrix).
	 */
	Clustering(BoolMatrix *clusterMatrix);
	virtual ~Clustering();

	void addCluster();
	void printClustering();
	int getN();

	/**
	 * Generates a l-neighborhood for this clustering.
	 * @return NeighborhoodList*
	 */
	NeighborhoodList* generateNeighborhood(int l);

private:
	/** number of nodes in the graph (n) */
	int numberOfNodes;
	/** number of clusters (k) */
	int numberOfClusters;
	/** the clustering matrix, with dimensions k x n */
	ClusterMatrixPtr clusterMatrixPtr;
	/** the l-neighborhood list of clusters */
	NeighborhoodListPtr neighborhoodListPtr;
};

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
