/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include <limits>
#include <iostream>

using namespace std;

namespace clusteringgraph {

Clustering::Clustering() : numberOfNodes(0), clusterList(),
		objectiveFunctionValue(0.0F) {

}

// TODO test dimension attribution
Clustering::Clustering(const Clustering& clustering, int n) : numberOfNodes(n),
		clusterList(), objectiveFunctionValue(clustering.objectiveFunctionValue) {
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

void Clustering::addCluster(SignedGraph& g, const int& i) {
	// 1. Create a new cluster in the list
	BoolArray array(MAX_NODES);

	// Add i to the newly created cluster
	// std::cout << "Adding vertex " << i << " to a new cluster."<< std::endl;
	array[i] = true;
	this->clusterList.push_back(array);
	this->numberOfNodes++;
	int k = this->clusterList.size() - 1;
	this->objectiveFunctionValue += calculateDeltaObjectiveFunction(g, k, i);
}

BoolArray& Clustering::getCluster(int clusterNumber) {
	return clusterList.at(clusterNumber);
}

void Clustering::addNodeToCluster(SignedGraph& g, int i, int k) {
	// std::cout << "Adding vertex " << i << " to cluster " << k << std::endl;
	BoolArray cluster = this->getCluster(k);
	cluster[i] = true;
	this->removeCluster(g, k);
	this->clusterList.push_back(cluster);
	this->numberOfNodes++;
	k = this->clusterList.size() - 1;
	this->objectiveFunctionValue += calculateDeltaObjectiveFunction(g, k, i);
}

template <typename T>
void remove(vector<T>* vec, size_t pos) {
    typename vector<T>::iterator it = vec->begin();
    advance(it, pos);
    vec->erase(it);
}

void Clustering::removeCluster(SignedGraph& g, int k) {
	// Swaps the k-th and the last element, to avoid linear-time removal
	//swap (clusterList[k], clusterList[clusterList.size() - 1]);
	//clusterList.erase(clusterList.end() - 1);
	remove <BoolArray> (&clusterList, k);
}

int Clustering::clusterSize(int k) {
	BoolArray cluster = this->getCluster(k);
	return cluster.count();
}

void Clustering::removeNodeFromCluster(SignedGraph& g, int i, int k) {
	BoolArray cluster = this->getCluster(k);
	int n = getNumberOfNodes();
	// verifica se o cluster eh unitario
	// TODO possivel otimizacao: verificar se pelo menos 2 bits estao setados
	this->objectiveFunctionValue -= calculateDeltaObjectiveFunction(g, k, i);
	if(clusterSize(k) == 1) {
		// cout << "Deleting cluster " << k << endl;
		this->removeCluster(g, k);
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
	std::cout << "Clustering configuration: I(P) = " <<
			this->objectiveFunctionValue << "\n";
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

// Calculates the delta of the objective function
int Clustering::calculateDeltaObjectiveFunction(SignedGraph& g, const int& k, const int& i) {
	int negativeSum = 0, positiveSum = 0;
	int n = g.getN();
	BoolArray cluster = getCluster(k);
	for(int b = 0; b < n; b++) {
		if(b != i) {
			if(cluster[b]) {
				// nodes i and b are in the same cluster
				// 1. calculates the change in the sum of internal
				//    negative edges (within the same cluster)
				if(g.getEdge(i, b) < 0) {
					negativeSum += (g.getEdge(i, b) * (-1));
				}
			} else {
				// nodes i and b are in different clusters
				// 2. calculates the change in the sum of external
				//    positive edges (within different clusters)
				if(g.getEdge(i, b) > 0) {
					positiveSum += g.getEdge(i, b);
				}
			}
		}
	}
	return negativeSum + positiveSum;
}

} /* namespace clusteringgraph */
