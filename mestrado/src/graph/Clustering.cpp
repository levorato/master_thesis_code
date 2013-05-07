/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"

namespace clusteringgraph {

Clustering::Clustering(int n, int k) {
	numberOfNodes = n;
	numberOfClusters = k;

	// Create the cluster array: n x k
	clusterMatrixPtr = new 2DBoolArray(boost::extents[k][n]);
}

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}

} /* namespace clusteringgraph */
