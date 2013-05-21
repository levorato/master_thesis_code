/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include <limits>

using namespace std;

namespace clusteringgraph {

Clustering::Clustering(int n) : numberOfNodes(n), clusterListPtr(new ClusterList()) {

}

// TODO test dimension attribution
Clustering::Clustering(Clustering* clustering, int numberOfNodes) : numberOfNodes(numberOfNodes),
		// deep copy of the clusterlist data
		clusterListPtr(new ClusterList(*(clustering->clusterListPtr.get()))) {

}

Clustering::~Clustering() {
	cout << "Freeing memory of Clustering object..." << endl;
}

int Clustering::getNumberOfNodes() {
	return numberOfNodes;
}

int Clustering::getNumberOfClusters() {
	return this->clusterListPtr->size();
}

void Clustering::addCluster(int vertexList[], unsigned int arraySize) {
	// 1. Create a new cluster in the list
	BoolArrayPtr array(new BoolArray(MAX_NODES));

	// 2. For every vertex in the list, remove the vertex from
	// any other cluster and add it to the newly created cluster
	int numberOfClusters = this->getNumberOfClusters();
	for(unsigned int i = 0; i < arraySize; i++) {
		std::cout << "Adding vertex " << vertexList[i] << " to cluster " << numberOfClusters << std::endl;
		int vertex = vertexList[i];
		for(int k = 0; k < numberOfClusters; k++) {
			(*clusterListPtr->at(k))[vertex] = false;
		}
		(*array)[vertex] = true;
	}
	this->clusterListPtr->push_back(array);
	this->numberOfNodes += arraySize;
}

BoolArray* Clustering::getCluster(int clusterNumber) {
	return clusterListPtr->at(clusterNumber).get();
}

void Clustering::addNodeToCluster(int i, int k) {
	std::cout << "Adding vertex " << i << " to cluster " << k << std::endl;
	BoolArray* cluster = this->getCluster(k);
	(*cluster)[i] = true;
	this->numberOfNodes++;
}

template <typename T>
void remove(vector<T>* vec, size_t pos) {
    typename vector<T>::iterator it = vec->begin();
    advance(it, pos);
    vec->erase(it);
}

void Clustering::removeCluster(int k) {
	ClusterList* ptr = clusterListPtr.get();
	remove <BoolArrayPtr> (ptr, k);
}

// TODO tratar o caso em que o cluster k desaparece
// Ainda esta mantendo o cluster vazio na lista de clusters
void Clustering::removeNodeFromCluster(int i, int k) {
	BoolArray* cluster = this->getCluster(k);
	// verifica se o cluster eh unitario
	if(cluster->size() == 1) {
		this->removeCluster(k);
	} else {
		(*cluster)[i] = false;
	}
	this->numberOfNodes--;
	std::cout << "Removing vertex " << i << " from cluster " << k << std::endl;
}

// TODO test this method
GainCalculation Clustering::gain(SignedGraph* graph, const int &a) {
	GainCalculation gainCalculation;
	ModularityMatrix* modularityMatrix = graph->getModularityMatrix();
	float max = (*modularityMatrix)[a][a];
	gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

	// For each cluster k...
	int nc = this->getNumberOfClusters();
	for(int k = 0; k < nc; k++) {
		int sum = 0;
		// Cluster(k)
		BoolArray cluster = *(this->clusterListPtr->at(k));
		// j in Cluster(k)
		for(int j = 0; j < this->numberOfNodes; j++) {
			if(cluster[j]) {
				sum += 2 * (*modularityMatrix)[a][j];
			}
		}
		sum += (*modularityMatrix)[a][a];
		if(sum > max) {
			max = sum;
			gainCalculation.clusterNumber = k;
		}
	}

	gainCalculation.value = max;
	return gainCalculation;
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
    	BoolArrayPtr arrayPtr = clusterListPtr->at(k);
    	for(int i = 0; i < MAX_NODES; i++) {
    		if((*arrayPtr)[i]) {
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
