/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include <limits>

namespace clusteringgraph {

Clustering::Clustering(int n) : numberOfNodes(n), clusterListPtr(new ClusterList()) {

}

// TODO test dimension attribution
Clustering::Clustering(Clustering* clustering, int numberOfNodes) : numberOfNodes(numberOfNodes),
		// deep copy of the clusterlist data
		clusterListPtr(new ClusterList(*(clustering->clusterListPtr.get()))) {

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
		std::cout << "Adding vertex " << vertexList[i] << " to cluster " << numberOfClusters << std::endl;
		int vertex = vertexList[i];
		for(int k = 0; k < numberOfClusters; k++) {
			clusterListPtr->at(k)[vertex] = false;
		}
		array[vertex] = true;
	}
	this->clusterListPtr->push_back(array);
}

const BoolArray& Clustering::getCluster(int clusterNumber) {
	return clusterListPtr->at(clusterNumber);
}

void Clustering::addNodeToCluster(int i, int k) {
	BoolArray cluster = this->getCluster(k);
	cluster[i] = true;
}

// TODO tratar o caso em que o cluster k desaparece
void Clustering::removeNodeFromCluster(int i, int k) {
	BoolArray cluster = this->getCluster(k);
	cluster[i] = false;
}

// TODO test this method
float Clustering::gain(SignedGraph* graph, const int &a) {
	ModularityMatrix* modularityMatrix = graph->getModularityMatrix();
	float max = (*modularityMatrix)[a][a];
	// For each cluster k...
	int nc = this->getNumberOfClusters();
	for(int k = 0; k < nc; k++) {
		int sum = 0;
		// Cluster(k)
		BoolArray cluster = this->clusterListPtr->at(k);
		// j in Cluster(k)
		for(int j = 0; j < this->numberOfNodes; j++) {
			if(cluster[j]) {
				sum += 2 * (*modularityMatrix)[a][j];
			}
		}
		sum += (*modularityMatrix)[a][a];
		if(sum > max) {
			max = sum ;
		}
	}
	return max;
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
