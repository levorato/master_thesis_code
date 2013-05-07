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

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}

int Clustering::getN() {
	return numberOfNodes;
}

void Clustering::printClustering() {
	std::cout << "Clustering configuration: TODO" << std::endl;
}

std::vector<Clustering*> Clustering::generateNeighborhood(Clustering* c, int l) {

}

} /* namespace clusteringgraph */
