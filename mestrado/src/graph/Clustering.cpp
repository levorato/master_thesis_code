/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

namespace clusteringgraph {

Clustering::Clustering() : clusterList(),
		imbalance(0.0, 0.0) {

}

// TODO test dimension attribution
Clustering::Clustering(const Clustering& clustering) :
		clusterList(), imbalance(clustering.imbalance) {
	// deep copy of the clusterlist data
	for(unsigned long i = 0; i < clustering.clusterList.size(); i++) {
		BoolArray boolArray = clustering.clusterList.at(i);
		this->clusterList.push_back(BoolArray(boolArray));
	}
}

Clustering::~Clustering() {
	// cout << "Freeing memory of Clustering object..." << endl;
}

unsigned long Clustering::getNumberOfClusters() {
	return this->clusterList.size();
}

BoolArray Clustering::addCluster(SignedGraph& g, const unsigned long& i) {
	// 1. Create a new cluster in the list
	BoolArray array(g.getN());

	// Add i to the newly created cluster
	// std::cout << "Adding vertex " << i << " to a new cluster."<< std::endl;
	array[i] = true;
	this->clusterList.push_back(array);
	this->imbalance += calculateDeltaObjectiveFunction(g, array, i);
	return array;
}

BoolArray& Clustering::getCluster(unsigned long clusterNumber) {
	return clusterList.at(clusterNumber);
}

void Clustering::addNodeToCluster(SignedGraph& g, unsigned long i, unsigned long k) {
	// std::cout << "Adding vertex " << i << " to cluster " << k << std::endl;
	this->getCluster(k)[i] = true;
	this->imbalance += calculateDeltaObjectiveFunction(g, this->getCluster(k), i);
}

void Clustering::removeCluster(SignedGraph& g, unsigned long k) {
	// Swaps the k-th and the last element, to avoid linear-time removal
	//swap (clusterList[k], clusterList[clusterList.size() - 1]);
	//clusterList.erase(clusterList.end() - 1);
	clusterList.erase(clusterList.begin()+k);
}

unsigned long Clustering::clusterSize(unsigned long k) {
	return this->getCluster(k).count();
}

void Clustering::removeNodeFromCluster(SignedGraph& g, unsigned long i, unsigned long k) {
	// verifica se o cluster eh unitario
	// TODO possivel otimizacao: verificar se pelo menos 2 bits estao setados
	// std::cout << "Removing vertex " << i << " from cluster " << k << std::endl;
	this->imbalance -= calculateDeltaObjectiveFunction(g, this->getCluster(k), i);
	if(clusterSize(k) == 1) {
		// cout << "Deleting cluster " << k << endl;
		this->removeCluster(g, k);
	} else {
		this->getCluster(k)[i] = false;
	}
}

void Clustering::printClustering() {
	std::cout << "Clustering configuration: I(P) = " << fixed << setprecision(2)
			<< this->imbalance.getValue() << "\n";
	print(std::cout, clusterList);
}

void Clustering::printClustering(ostream& os) {
	os << "Clustering configuration: " << std::endl;
	print(os, clusterList);
}

void Clustering::print(std::ostream& os, ClusterList& l)
{
	unsigned long numberOfClusters = l.size();
	for(unsigned long k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << " (" << clusterSize(k) <<  "): [ ";
    	BoolArray array = l.at(k);
    	for(unsigned long i = 0; i < array.size(); i++) {
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

// TODO verificar se essa igualdade funciona
bool Clustering::equals(Clustering& c) {
	if(std::equal(&this->clusterList, &this->clusterList +
			this->clusterList.size(), &c.clusterList))
		return true;
	else
		return false;
}

// Calculates the delta of the objective function
Imbalance Clustering::calculateDeltaObjectiveFunction(SignedGraph& g, BoolArray& cluster, const unsigned long& i) {
	double negativeSum = 0, positiveSum = 0;
	unsigned long n = g.getN();
	for(unsigned long b = 0; b < n; b++) {
		if(b != i) {
			if(cluster[b]) {
				// nodes i and b are in the same cluster
				// 1. calculates the change in the sum of internal
				//    negative edges (within the same cluster)
				if(g.isNegativeEdge(i, b)) {
					negativeSum += fabs(g.getEdge(i, b));
				}
				if(g.isNegativeEdge(b, i)) {
					negativeSum += fabs(g.getEdge(b, i));
				}
			} else {
				// nodes i and b are in different clusters
				// 2. calculates the change in the sum of external
				//    positive edges (within different clusters)
				if(g.isPositiveEdge(i, b)) {
					positiveSum += g.getEdge(i, b);
				}
				if(g.isPositiveEdge(b, i)) {
					positiveSum += g.getEdge(b, i);
				}
			}
		}
	}

	return Imbalance(positiveSum, negativeSum);
}

} /* namespace clusteringgraph */
