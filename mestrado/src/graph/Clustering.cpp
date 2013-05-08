/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"

namespace clusteringgraph {

Clustering::Clustering(int n) : numberOfNodes(n), numberOfClusters(0),
		clusterMatrixPtr(), neighborhoodListPtr() {

}

Clustering::Clustering(int n, int k) : numberOfNodes(n), numberOfClusters(k),
		clusterMatrixPtr(new BoolMatrix(boost::extents[k][n])), neighborhoodListPtr() {
	// Create the cluster array: n x k, filled with zeroes
	std::fill(clusterMatrixPtr->origin(), clusterMatrixPtr->origin() +
			clusterMatrixPtr->size(), false);
}

// TODO test dimension attribution
Clustering::Clustering(BoolMatrix* clusterMatrix) : clusterMatrixPtr(clusterMatrix),
		neighborhoodListPtr() {
	numberOfClusters = clusterMatrixPtr->shape()[0];
	numberOfNodes = clusterMatrixPtr->shape()[1];
}

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}

int Clustering::getNumberOfNodes() {
	return numberOfNodes;
}

int Clustering::getNumberOfClusters() {
	return numberOfClusters;
}

void Clustering::addCluster(int vertexList[], unsigned int arraySize) {
	BoolMatrix m = *(clusterMatrixPtr.get());

	// 1. Create a new line for the new cluster in the matrix
	// the method below keeps existing array data
	if(numberOfClusters == 0) {
		clusterMatrixPtr.reset(new BoolMatrix(boost::extents[numberOfClusters++][numberOfNodes]));
	} else {
		clusterMatrixPtr->resize(boost::extents[numberOfClusters++][numberOfNodes]);
		for(int i = 0; i < numberOfNodes; i++) {
			m[numberOfClusters - 1][i] = false;
		}
	}

	// 2. For every vertex in the list, remove the vertex from
	// any other cluster and add it to the newly created cluster
	for(unsigned int i = 0; i < arraySize; i++) {
		int vertex = vertexList[i];
		for(int k = 0; k < numberOfClusters; k++) {
			m[k][vertex] = false;
		}
		m[numberOfClusters - 1][vertex] = true;
	}
}

void Clustering::printClustering() {
	std::cout << "Clustering configuration: " << std::endl;
	print(std::cout, *this->clusterMatrixPtr.get());
}

void Clustering::print(std::ostream& os, const BoolMatrix& A)
{
	int numberOfClusters = A.shape()[0];
	int numberOfNodes = A.shape()[1];
    for(int k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << ": [";
    	for(int i = 0; i < numberOfNodes; i++) {
    		if(A[k][i]) {
    			os << i << " ";
    		}
    	}
    	os << "] \n";
    }
}

// TODO: Implement this
NeighborhoodList* Clustering::generateNeighborhood(int l) {
	return NULL;
}

bool Clustering::equals(Clustering* c) {
	if(std::equal(this->clusterMatrixPtr.get(), this->clusterMatrixPtr.get() +
			this->clusterMatrixPtr->size(), c->clusterMatrixPtr.get()))
		return true;
	else
		return false;
}

} /* namespace clusteringgraph */
