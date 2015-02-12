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

const long Clustering::NEW_CLUSTER;
const long Clustering::NO_CLUSTER;

Clustering::Clustering() : clusterArray(), clusterSize(),
		imbalance(0.0, 0.0), problemType(0), positiveSum(),
		negativeSum()
{

}

Clustering::Clustering(SignedGraph& g) : clusterArray(g.getN(), NO_CLUSTER), clusterSize(),
		imbalance(0.0, 0.0), problemType(0), positiveSum(), negativeSum() {

}

Clustering::Clustering(const Clustering& clustering) :
		clusterArray(clustering.clusterArray.begin(), clustering.clusterArray.end()),
		clusterSize(clustering.clusterSize.begin(), clustering.clusterSize.end()),
		imbalance(clustering.imbalance), problemType(clustering.problemType),
		positiveSum(clustering.positiveSum),
		negativeSum(clustering.negativeSum) {

}

Clustering::Clustering(ClusterArray &cArray, SignedGraph& g, ClusteringProblem &p) : clusterArray(cArray),
		clusterSize(), imbalance(0.0, 0.0), problemType(p.getType()), positiveSum(), negativeSum(){
	ClusterArray::iterator pos = std::max_element(cArray.begin(), cArray.end());
	long numberOfClusters = (*pos) + 1;
	// cout << "number of clusters = " << numberOfClusters << endl;
	std::vector< std::vector<long> > clusters(numberOfClusters, std::vector<long>());

	long n = g.getN();
	for(unsigned long i = 0; i < n; i++) {
		assert(clusterArray[i] < numberOfClusters);
		clusters[clusterArray[i]].push_back(i);
	}
	// compute clusters' size
	for(long k = 0; k < numberOfClusters; k++) {
		clusterSize.push_back(clusters[k].size());
		// cout << "Size of cluster " << k << " is " << clusterSize[k] << endl;
	}
	// compute the imbalance
	this->setImbalance(p.objectiveFunction(g, *this));
}

Clustering::~Clustering() {
	// cout << "Freeing memory of Clustering object..." << endl;
}

unsigned long Clustering::getNumberOfClusters() const {
	return this->clusterSize.size();
}

void Clustering::addCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, bool updateImbalance) {
	assert(clusterArray.size() > 0);
	// 1. Increase the number of clusters
	this->clusterSize.push_back(1);

	// Adds i to the newly created cluster
	// BOOST_LOG_TRIVIAL(trace) <<  "Adding vertex " << i << " to a new cluster.";
	this->clusterArray[i] = this->getNumberOfClusters() - 1;

	if(updateImbalance) {
		if(p.getType() == ClusteringProblem::CC_PROBLEM) {
			this->imbalance += p.calculateDeltaPlusObjectiveFunction(g, *this, clusterArray.size()-1, i);
		} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
			this->imbalance = p.calculateDeltaPlusObjectiveFunction(g, *this, clusterArray.size()-1, i);
		}
	}
}

int Clustering::getBiggestClusterIndex() {
	std::vector<unsigned long>::iterator pos = std::max_element(clusterSize.begin(), clusterSize.end());
	return std::distance(clusterSize.begin(), pos);
}

void Clustering::addNodeToCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, const unsigned long& k,
		bool updateImbalance) {
	assert(clusterArray.size() > 0);
	//BOOST_LOG_TRIVIAL(trace) << "Adding vertex " << i << " to cluster " << k;
	this->clusterArray[i] = k;
	this->clusterSize[k]++;

	if(updateImbalance) {
		if(p.getType() == ClusteringProblem::CC_PROBLEM) {
			this->imbalance += p.calculateDeltaPlusObjectiveFunction(g, *this, k, i);
		} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
			this->imbalance = p.calculateDeltaPlusObjectiveFunction(g, *this, k, i);
		}
	}
}

void Clustering::removeCluster(SignedGraph& g, unsigned long k) {
	assert(clusterArray.size() > 0);
	// clusterArray.erase(clusterArray.begin()+k);
	// TODO complete code
	// only if cluster has been removed
	clusterSize.erase(clusterSize.begin()+k);
	// if cluster k was removed, all cluster numbers above k
	// must be subtracted by one
	long n = g.getN();
	for(long i = 0; i < n; i++) {
		if(clusterArray[i] > k) {
			clusterArray[i]--;
		}
	}
}

unsigned long Clustering::getClusterSize(unsigned long k) {
	return this->clusterSize[k];
}

void Clustering::setClusterSize(unsigned long k, unsigned long size) {
	this->clusterSize[k] = size;
}

void Clustering::removeNodeFromCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, 
		const unsigned long& k, bool updateImbalance) {
	assert(clusterArray.size() > 0);
	// verifica se o cluster eh unitario
	//BOOST_LOG_TRIVIAL(trace) << "Removing vertex " << i << " from cluster " << k;
	if(updateImbalance) {
		if(p.getType() == ClusteringProblem::CC_PROBLEM) {
			this->imbalance -= p.calculateDeltaMinusObjectiveFunction(g, *this, k, i);
		} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
			this->imbalance = p.calculateDeltaMinusObjectiveFunction(g, *this, k, i);
		}
	}
	this->clusterArray[i] = NO_CLUSTER;
	if(clusterSize[k] <= 1) {
		//BOOST_LOG_TRIVIAL(trace) << "Deleting cluster " << k;
		this->removeCluster(g, k);
	} else {
		this->clusterSize[k]--;
	}
}

void Clustering::printClustering(unsigned long n) {
	BOOST_LOG_TRIVIAL(trace) << "Clustering configuration: I(P) = " << fixed << setprecision(2)
			<< this->imbalance.getValue();
	BOOST_LOG_TRIVIAL(trace) << toString(n);
}

void Clustering::printClustering(ostream& os, unsigned long n) {
	os << "Clustering configuration: " << std::endl;
	print(os, clusterArray, n);
}

void Clustering::print(std::ostream& os, ClusterArray& a, unsigned long n)
{
	unsigned long numberOfClusters = *std::max_element(a.begin(), a.end()) + 1;
	std::vector< std::vector<long> > clusters(numberOfClusters, std::vector<long>());

	for(unsigned long i = 0; i < n; i++) {
		clusters[clusterArray[i]].push_back(i);
	}

	for(unsigned long k = 0; k < numberOfClusters; k++) {
    	os << " Partition " << k << " (" << getClusterSize(k) <<  "): [ ";
    	for(unsigned long i = 0; i < clusters[k].size(); i++) {
    		os << clusters[k][i] << " ";
    	}
    	os << "] \n";
    }
}

string Clustering::toString(unsigned long n) {
	stringstream ss;
	printClustering(ss, n);
	return ss.str();
}

// TODO verificar se essa igualdade funciona
bool Clustering::equals(Clustering& c) {
	if(std::equal(&this->clusterArray, &this->clusterArray +
			this->clusterArray.size(), &c.clusterArray))
		return true;
	else
		return false;
}

} /* namespace clusteringgraph */
