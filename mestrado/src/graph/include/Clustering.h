/*
 * Clustering.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <boost/config.hpp>
#include <boost/multi_array.hpp>
#include <boost/scoped_ptr.hpp>

namespace clusteringgraph {

typedef boost::multi_array<bool, 2> BoolMatrix;
typedef boost::scoped_ptr<BoolMatrix> ClusterMatrixPtr;

class Clustering {
public:
	Clustering(int n, int k);
	virtual ~Clustering();

private:
	/** number of nodes in the graph (n) */
	int numberOfNodes;
	/** number of clusters (k) */
	int numberOfClusters;
	/** Create the clustering matrix, with dimensions k x n */
	ClusterMatrixPtr clusterMatrixPtr;
};

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
