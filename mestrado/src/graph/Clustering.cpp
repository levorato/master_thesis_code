/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"

namespace clusteringgraph {

Clustering::Clustering(int n, int k) : numberOfNodes(n), numberOfClusters(k),
		clusterMatrixPtr(new BoolMatrix(boost::extents[k][n])) {
	// Create the cluster array: n x k
}

// TODO test dimension attribution
Clustering::Clustering(BoolMatrix* clusterMatrix) : clusterMatrixPtr(clusterMatrix) {
	numberOfClusters = clusterMatrixPtr->shape()[0];
	numberOfNodes = clusterMatrixPtr->shape()[1];
}

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}

int Clustering::getN() {
	return numberOfNodes;
}

void Clustering::addCluster() {

}

void Clustering::printClustering() {
	std::cout << "Clustering configuration: TODO" << std::endl;
}

// TODO: Implement this
NeighborhoodList* Clustering::generateNeighborhood(int l) {
	return NULL;
}

} /* namespace clusteringgraph */
