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

#include <boost/log/trivial.hpp>
#include <boost/graph/adjacency_matrix.hpp>


using namespace std;
using namespace boost;
using namespace problem;

namespace clusteringgraph {

Clustering::Clustering(int n) : vertexInClusterId(n, 0), clusterIdList(), clusterIdCount(0),
		clusterSize(), imbalance(0.0, 0.0), problemType(0), positiveSum(),
		negativeSum()
{

}

Clustering::Clustering(const Clustering& clustering) :
		vertexInClusterId(clustering.vertexInClusterId.begin(), clustering.vertexInClusterId.end()),
		clusterIdList(clustering.clusterIdList.begin(), clustering.clusterIdList.end()),
		clusterIdCount(clustering.clusterIdCount),
		clusterSize(clustering.clusterSize.begin(), clustering.clusterSize.end()),
		imbalance(clustering.imbalance),
		problemType(clustering.problemType),
		positiveSum(clustering.positiveSum),
		negativeSum(clustering.negativeSum) {

}

Clustering::~Clustering() {
	// cout << "Freeing memory of Clustering object..." << endl;
}

unsigned long Clustering::getNumberOfClusters() const {
	return clusterIdList.size();
}

void Clustering::addCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i) {
	// 1. Adds a new cluster to clusterIdList
	clusterIdList.push_back(clusterIdCount++);
	clusterSize.push_back(1);

	// Adds i to the newly created cluster
	//BOOST_LOG_TRIVIAL(trace) <<  "Adding vertex " << i << " to a new cluster.";
	this->vertexInClusterId.at(i) = clusterIdCount - 1;
	if(p.getType() == ClusteringProblem::CC_PROBLEM) {
		this->imbalance += p.calculateDeltaPlusObjectiveFunction(g, *this, clusterIdCount-1, i);
	} else if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		this->imbalance = p.calculateDeltaPlusObjectiveFunction(g, *this, clusterIdCount-1, i);
	}
}

void Clustering::addNodeToCluster(SignedGraph& g, ClusteringProblem& p, const unsigned long& i, const unsigned long& k) {
	//BOOST_LOG_TRIVIAL(trace) << "Adding vertex " << i << " to cluster " << k;
	this->vertexInClusterId.at(i) = k;
	clusterSize.at(k)++;
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
	this->clusterIdList.erase(this->clusterIdList.begin()+k);
}

unsigned long Clustering::getClusterSize(unsigned long k) {
	return this->clusterSize.at(k);
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
	clusterSize.at(k)--;
	if(clusterSize.at(k) == 0) {
		//BOOST_LOG_TRIVIAL(trace) << "Deleting cluster " << k;
		this->removeCluster(g, k);
	}
}

void Clustering::printClustering() {
	//BOOST_LOG_TRIVIAL(trace) << "Clustering configuration: I(P) = " << fixed << setprecision(2)
			//<< this->imbalance.getValue();
	//BOOST_LOG_TRIVIAL(trace) << toString();
}

void Clustering::printClustering(ostream& os) {
	os << "Clustering configuration: " << std::endl;
	// TODO recreate clustering printing method based on new data structure!
	// print(os, clusterList);
}

void Clustering::print(std::ostream& os, ClusterList& l)
{
	/*
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
    */
}

string Clustering::toString() {
	stringstream ss;
	printClustering(ss);
	return ss.str();
}

// TODO verificar se essa igualdade funciona
// TODO recreate clustering equals method based on new data structure!
bool Clustering::equals(Clustering& c) {
	/*
	if(std::equal(&this->clusterList, &this->clusterList +
			this->clusterList.size(), &c.clusterList))
		return true;
	else
	*/
		return false;
}

} /* namespace clusteringgraph */
