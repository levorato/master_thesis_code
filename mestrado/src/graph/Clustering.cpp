/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include <limits>

namespace clusteringgraph {

Clustering::Clustering(int n) : numberOfNodes(n), clusterListPtr(new ClusterList()),
		modularityMatrix(boost::extents[n][n]) {

}

// TODO test dimension attribution
Clustering::Clustering(ClusterList* clusterList, int numberOfNodes) : numberOfNodes(numberOfNodes),
		clusterListPtr(clusterList), modularityMatrix(boost::extents[numberOfNodes][numberOfNodes]) {

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

// TODO calculate the modularity matrix of weighed graphs
void Clustering::calculateModularityMatrix(SignedGraph* g) {
	int m = g->getM();
	int degree[g->getN()];
	// Prestore the degrees for optimezed lookup
	for(int i = 0; i < this->numberOfNodes; i++) {
		degree[i] = g->getDegree(i);
	}

	for(int i = 0; i < this->numberOfNodes; i++) {
		for(int j = 0; j < this->numberOfNodes; j++) {
			int a = (g->getEdge(i, j) != 0) ? 1 : 0;
			modularityMatrix[i][j] = a - ( (degree[i] * degree[j]) / (2 * m) );
		}
	}
}

// TODO test this method
float Clustering::gain(int a) {
	float max = modularityMatrix[a][a];
	// For each cluster k...
	int nc = this->getNumberOfClusters();
	for(int k = 0; k < nc; k++) {
		int sum = 0;
		// Cluster(k)
		BoolArray cluster = this->clusterListPtr->at(k);
		// j in Cluster(k)
		for(int j = 0; j < this->numberOfNodes; j++) {
			if(cluster[j]) {
				sum += 2 * modularityMatrix[a][j];
			}
		}
		sum += modularityMatrix[a][a];
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
