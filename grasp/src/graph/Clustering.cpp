/*
 * Clustering.cpp
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#include "include/Clustering.h"
#include "../problem/include/ClusteringProblemFactory.h"
#include "include/Graph.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>

#include <boost/log/trivial.hpp>
#include <boost/graph/adjacency_matrix.hpp>


using namespace std;
using namespace boost;
using namespace problem;

namespace clusteringgraph {

Clustering::Clustering() : clusterList(), clusterSize(),
		imbalance(0.0, 0.0), problemType(0), positiveSum(),
		negativeSum()
{

}

Clustering::Clustering(const Clustering& clustering) :
		clusterList(clustering.clusterList.begin(), clustering.clusterList.end()),
		clusterSize(clustering.clusterSize.begin(), clustering.clusterSize.end()),
		imbalance(clustering.imbalance), problemType(clustering.problemType),
		positiveSum(clustering.positiveSum),
		negativeSum(clustering.negativeSum) {

}

Clustering::~Clustering() {
	// cout << "Freeing memory of Clustering object..." << endl;
}

unsigned long Clustering::getNumberOfClusters() const {
	return this->clusterList.size();
}

BoolArray Clustering::addCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i) {
	// 1. Create a new cluster in the list
	BoolArray array(g.getN());

	// Add i to the newly created cluster
	////BOOST_LOG_TRIVIAL(trace) <<  "Adding vertex " << i << " to a new cluster.";
	array[i] = true;
	this->clusterList.push_back(array);
	this->clusterSize.push_back(1);
	if(p.getType() == ClusteringProblem::CC_PROBLEM) {
		this->imbalance += p.calculateDeltaPlusObjectiveFunction(g, *this, clusterList.size()-1, i);
	} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		this->imbalance = p.calculateDeltaPlusObjectiveFunction(g, *this, clusterList.size()-1, i);
	}
	return array;
}

void Clustering::addCluster(BoolArray b) {
	clusterList.push_back(b);
	clusterSize.push_back(b.count());
}

BoolArray& Clustering::getCluster(unsigned long clusterNumber) {
	return clusterList.at(clusterNumber);
}

void Clustering::setCluster(unsigned long clusterNumber, BoolArray b) {
	clusterList[clusterNumber] = b;
}

int Clustering::getBiggestClusterIndex() {
	ClusterList::iterator pos = std::min_element(clusterList.begin(), clusterList.end(),
			ClusterSizeComparison(false));
	return std::distance(clusterList.begin(), pos);
}

void Clustering::addNodeToCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, const unsigned long& k) {
	//BOOST_LOG_TRIVIAL(trace) << "Adding vertex " << i << " to cluster " << k;
	this->getCluster(k)[i] = true;
	this->clusterSize[k]++;
	if(p.getType() == ClusteringProblem::CC_PROBLEM) {
		this->imbalance += p.calculateDeltaPlusObjectiveFunction(g, *this, k, i);
	} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		this->imbalance = p.calculateDeltaPlusObjectiveFunction(g, *this, k, i);
	}
}

void Clustering::removeCluster(SignedGraph& g, unsigned long k) {
	// Swaps the k-th and the last element, to avoid linear-time removal
	//swap(clusterList[k], clusterList[clusterList.size() - 1]);
	//clusterList.erase(clusterList.end() - 1);
	clusterList.erase(clusterList.begin()+k);
	clusterSize.erase(clusterSize.begin()+k);
}

unsigned long Clustering::getClusterSize(unsigned long k) {
	return this->clusterSize[k];
}

void Clustering::setClusterSize(unsigned long k, unsigned long size) {
	this->clusterSize[k] = size;
}

void Clustering::removeNodeFromCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, const unsigned long& k) {
	// verifica se o cluster eh unitario
	// TODO possivel otimizacao: verificar se pelo menos 2 bits estao setados
	//BOOST_LOG_TRIVIAL(trace) << "Removing vertex " << i << " from cluster " << k;
	if(p.getType() == ClusteringProblem::CC_PROBLEM) {
		this->imbalance -= p.calculateDeltaMinusObjectiveFunction(g, *this, k, i);
	} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		this->imbalance = p.calculateDeltaMinusObjectiveFunction(g, *this, k, i);
	}
	if(clusterSize[k] <= 1) {
		//BOOST_LOG_TRIVIAL(trace) << "Deleting cluster " << k;
		this->removeCluster(g, k);
	} else {
		this->getCluster(k)[i] = false;
		this->clusterSize[k]--;
	}
}

void Clustering::printClustering() {
	BOOST_LOG_TRIVIAL(trace) << "Clustering configuration: I(P) = " << fixed << setprecision(2)
			<< this->imbalance.getValue();
	BOOST_LOG_TRIVIAL(trace) << toString();
}

void Clustering::printClustering(ostream& os) {
	os << "Clustering configuration: " << std::endl;
	print(os, clusterList);
}

void Clustering::print(std::ostream& os, ClusterList& l)
{
	unsigned long numberOfClusters = l.size();
	for(unsigned long k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << " (" << getClusterSize(k) <<  "): [ ";
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

} /* namespace clusteringgraph */
