/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"

namespace clusteringgraph {

Clustering::Clustering(int n) : numberOfNodes(n),
		clusterListPtr(new ClusterList()) {

}

// TODO test dimension attribution
Clustering::Clustering(ClusterList* clusterList, int numberOfNodes) : numberOfNodes(numberOfNodes),
		clusterListPtr(clusterList) {

}

Clustering::~Clustering() {
	// TODO Auto-generated destructor stub
}

int Clustering::getNumberOfNodes() {
	return numberOfNodes;
}

int Clustering::getNumberOfClusters() {
	return this->clusterListPtr->size();
}

void Clustering::addCluster(int vertexList[], unsigned int arraySize) {
	// 1. Create a new cluster in the list
	BoolArray array(MAX_NODES);

	// 2. For every vertex in the list, remove the vertex from
	// any other cluster and add it to the newly created cluster
	int numberOfClusters = this->getNumberOfClusters();
	for(unsigned int i = 0; i < arraySize; i++) {
		std::cout << "Adding vertex " << vertexList[i] << " to cluster " << (numberOfClusters - 1) << std::endl;
		int vertex = vertexList[i];
		for(int k = 0; k < numberOfClusters; k++) {
			clusterListPtr->at(k)[vertex] = false;
		}
		array[vertex] = true;
	}
	this->clusterListPtr->push_back(array);
}

BoolArray Clustering::getCluster(int clusterNumber) {
	return clusterListPtr->at(clusterNumber);
}

void Clustering::printClustering() {
	std::cout << "Clustering configuration: " << std::endl;
	print(std::cout, clusterListPtr.get());
}

void Clustering::print(std::ostream& os, ClusterList* l)
{
	int numberOfClusters = l->size();
	for(int k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << ": [ ";
    	for(int i = 0; i < numberOfNodes; i++) {
    		if(clusterListPtr->at(k)[i]) {
    			os << i << " ";
    		}
    	}
    	os << "] \n";
    }
}

// TODO verificar se essa igualdade funciona
bool Clustering::equals(Clustering* c) {
	if(std::equal(this->clusterListPtr.get(), this->clusterListPtr.get() +
			this->clusterListPtr->size(), c->clusterListPtr.get()))
		return true;
	else
		return false;
}

} /* namespace clusteringgraph */
