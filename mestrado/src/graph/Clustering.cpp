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

Clustering::Clustering() : numberOfNodes(0), clusterList(),
		objectiveFunctionValue(0.0F) {

}

// TODO test dimension attribution
Clustering::Clustering(const Clustering& clustering, int n) : numberOfNodes(n),
		clusterList(), objectiveFunctionValue(0.0F) {
	// deep copy of the clusterlist data
	for(unsigned int i = 0; i < clustering.clusterList.size(); i++) {
		BoolArray boolArray = clustering.clusterList.at(i);
		this->clusterList.push_back(BoolArray(boolArray));
	}
}

Clustering::~Clustering() {
	// cout << "Freeing memory of Clustering object..." << endl;
}

int Clustering::getNumberOfNodes() {
	return numberOfNodes;
}

int Clustering::getNumberOfClusters() {
	return this->clusterList.size();
}

void Clustering::addCluster(int vertexList[], unsigned int arraySize) {
	// 1. Create a new cluster in the list
	BoolArray array(MAX_NODES);

	// 2. For every vertex in the list, add it to the newly created cluster
	for(unsigned int i = 0; i < arraySize; i++) {
		// std::cout << "Adding vertex " << vertexList[i] << " to cluster " << numberOfClusters << std::endl;
		int vertex = vertexList[i];
		array[vertex] = true;
	}
	this->clusterList.push_back(array);
	this->numberOfNodes += arraySize;
}

BoolArray& Clustering::getCluster(int clusterNumber) {
	return clusterList.at(clusterNumber);
}

void Clustering::addNodeToCluster(int i, int k) {
	// std::cout << "Adding vertex " << i << " to cluster " << k << std::endl;
	BoolArray cluster = this->getCluster(k);
	cluster[i] = true;
	this->removeCluster(k);
	this->clusterList.push_back(cluster);
	this->numberOfNodes++;
}

template <typename T>
void remove(vector<T>* vec, size_t pos) {
    typename vector<T>::iterator it = vec->begin();
    advance(it, pos);
    vec->erase(it);
}

void Clustering::removeCluster(int k) {
	// Swaps the k-th and the last element, to avoid linear-time removal
	//swap (clusterList[k], clusterList[clusterList.size() - 1]);
	//clusterList.erase(clusterList.end() - 1);
	remove <BoolArray> (&clusterList, k);
}

int Clustering::clusterSize(int k) {
	BoolArray cluster = this->getCluster(k);
	return cluster.count();
}

void Clustering::removeNodeFromCluster(int i, int k) {
	BoolArray cluster = this->getCluster(k);
	// verifica se o cluster eh unitario
	// TODO possivel otimizacao: verificar se pelo menos 2 bits estao setados
	if(clusterSize(k) == 1) {
		// cout << "Deleting cluster " << k << endl;
		this->removeCluster(k);
	} else {
		cluster[i] = false;
	}
	this->numberOfNodes--;
	// std::cout << "Removing vertex " << i << " from cluster " << k << std::endl;
}

// TODO test this method
GainCalculation Clustering::gain(SignedGraph& graph, const int &a) {
	GainCalculation gainCalculation;
	ModularityMatrix& modularityMatrix = graph.getModularityMatrix();
	float max = modularityMatrix[a][a];
	gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

	// For each cluster k...
	int nc = this->getNumberOfClusters();
	for(int k = 0; k < nc; k++) {
		int sum = 0;
		// Cluster(k)
		BoolArray cluster = this->clusterList.at(k);
		// j in Cluster(k)
		for(int j = 0; j < this->numberOfNodes; j++) {
			if(cluster[j]) {
				sum += 2 * modularityMatrix[a][j];
			}
		}
		sum += modularityMatrix[a][a];
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
	print(std::cout, clusterList);
}

void Clustering::printClustering(ostream& os) {
	os << "Clustering configuration: " << std::endl;
	print(os, clusterList);
}

void Clustering::print(std::ostream& os, ClusterList& l)
{
	int numberOfClusters = l.size();
	for(int k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << " (" << clusterSize(k) <<  "): [ ";
    	BoolArray array = l.at(k);
    	for(int i = 0; i < MAX_NODES; i++) {
    		if(array[i]) {
    			os << i << " ";
    		}
    	}
    	os << "] \n";
    }
}

string Clustering::toString() {
	stringstream ss;
	printClustering(ss);
	return ss.str();
}

float Clustering::getObjectiveFunctionValue() {
	return objectiveFunctionValue;
}

void Clustering::setObjectiveFunctionValue(float f) {
	objectiveFunctionValue = f;
}

// TODO verificar se essa igualdade funciona
bool Clustering::equals(Clustering& c) {
	if(std::equal(&this->clusterList, &this->clusterList +
			this->clusterList.size(), &c.clusterList))
		return true;
	else
		return false;
}

} /* namespace clusteringgraph */
