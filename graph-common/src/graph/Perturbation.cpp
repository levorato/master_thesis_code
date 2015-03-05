/*
 * Perturbation.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: mlevorato
 */

#include "include/Perturbation.h"
#include "../util/include/RandomUtil.h"

#include <boost/log/trivial.hpp>
#include <vector>
#include <algorithm>
#include <iostream>

namespace clusteringgraph {

using namespace std;
using namespace util;
using namespace boost;

Perturbation::~Perturbation() {
	// TODO Auto-generated destructor stub
}

Clustering Perturbation::randomMove(SignedGraph* g, Clustering clustering, ClusteringProblem& p,
		unsigned long numberOfMoves) {

	// BOOST_LOG_TRIVIAL(debug)<< "Generating perturbation of level " << numberOfMoves;
	Clustering c = clustering;
	RandomUtil randomUtil;
	randomUtil.setSeed(this->_randomSeed);

	int n = g->getN();
	int nc = c.getNumberOfClusters();
	std::vector<int> nodeList(n, 0);
	std::vector<int> k2List(nc + numberOfMoves, 0);
	for(int i = 0; i < n; i++) {
		nodeList[i] = i;
	}
	std::random_shuffle(nodeList.begin(), nodeList.end());
	for(int i = 0; i < numberOfMoves; i++) {
		nc = c.getNumberOfClusters();
		int k2 = randomUtil.next(0, nc - 1);
		int idx = randomUtil.next(0, n - 1);
		// cout << "k2 = " << k2 << endl;
		c = move1optCCProblem(g, c, p, idx, k2);
	}
	return c;
}

Clustering Perturbation::move1optCCProblem(SignedGraph* g, Clustering clustering, ClusteringProblem& p, int node, int k2) {
	Clustering cTemp = clustering;
	int nc = cTemp.getNumberOfClusters();
	int n = g->getN();
	// BOOST_LOG_TRIVIAL(trace)<< "Random move 1-opt .";
	// cTemp.printClustering();
	// (A) ==== Random node i from cluster k1 ====
	// (B) ==== Random cluster k2 ====
	ClusterArray cluster = clustering.getClusterArray();
	int k1 = cluster[node];
	if(k2 == k1) {
		k2 = Clustering::NEW_CLUSTER;
	}
	// removes node from cluster k1 and inserts in cluster k2
	cTemp.removeNodeFromCluster(*g, p, node, k1);
	if(k2 >= 0) {
		// recalculates the number of clusters, as one of them may have been removed
		if((cTemp.getNumberOfClusters() < nc) && (k2 >= k1)) {
			// cluster k1 has been removed
			cTemp.addNodeToCluster(*g, p, node, k2 - 1);
		} else {
			cTemp.addNodeToCluster(*g, p, node, k2);
		}
	} else {  // inserts node into a new cluster
		cTemp.addCluster(*g, p, node);
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Generated cluster:";
	// cTemp.printClustering();
	return cTemp;
}

Clustering Perturbation::randomMove1opt(SignedGraph* g, Clustering clustering, ClusteringProblem& p) {
	Clustering cTemp = clustering;
	// shuffles 2 clusters (origin and destination) and a vertex 'node'
	int k1 = 0, k2 = 0, node = 0;
	int nc = cTemp.getNumberOfClusters();
	int n = g->getN();
	// BOOST_LOG_TRIVIAL(trace)<< "Random move 1-opt .";
	// cTemp.printClustering();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// (A) ==== Random cluster 1 ====
	k1 = randomUtil.next(0, nc - 1);
	// Avoid removing a cluster if RCCProblem
	while((p.getType() == ClusteringProblem::RCC_PROBLEM) and (cTemp.getClusterSize(k1) == 1)) {
		// if RCCProblem, gets a new random cluster whose size is bigger than one (to avoid removing a cluster)
		k1 = randomUtil.next(0, nc - 1);
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Random k1 = " << k1;
	// (B) ==== Random node i from cluster 1 ====
	int idx_i = randomUtil.next(0, cTemp.getClusterSize(k1) - 1);
	ClusterArray myCluster = clustering.getClusterArray();
	for(int j = 0, count = 0; j < n; j++) {
		if(myCluster[j] == k1) {
			if(count == idx_i) {
				node = j;
				break;
			}
			count++;
		}
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Random node = " << node;
	// (C) ==== Random cluster 2 ====
	int startc = -1;  // cluster == -1 means new cluster
	if(p.getType() == ClusteringProblem::RCC_PROBLEM) {  // avoid creating a new cluster if RCCProblem
		startc = 0;
	}
	k2 = randomUtil.next(startc, nc - 1);
	// k2 must be different from k1
	while(k2 == k1) {
		k2 = randomUtil.next(startc, nc - 1);
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Random k2 = " << k2;
	// removes node from cluster1 and inserts in cluster2
	cTemp.removeNodeFromCluster(*g, p, node, k1);
	if(k2 >= 0) {
		// recalculates the number of clusters, as one of them may have been removed
		if((cTemp.getNumberOfClusters() < nc) && (k2 >= k1)) {
			// cluster k1 has been removed
			cTemp.addNodeToCluster(*g, p, node, k2 - 1);
		} else {
			cTemp.addNodeToCluster(*g, p, node, k2);
		}
	} else {  // inserts node into a new cluster
		cTemp.addCluster(*g, p, node);
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Generated cluster:";
	// cTemp.printClustering();
	return cTemp;
}

std::pair<int, int> Perturbation::randomMove1optCCProblem(int n, int nc) {
	// BOOST_LOG_TRIVIAL(trace)<< "Random move 1-opt .";
	// cTemp.printClustering();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// (A) ==== Random node i from any cluster ====
	int node = randomUtil.next(0, n - 1);
	// BOOST_LOG_TRIVIAL(trace)<< "Random node = " << node;
	// (B) ==== Random cluster k2 ====
	// k2 is a cluster offset to which cluster node i is being moved to
	// if k2 == cluster[node] then node is moved to a new cluster
	int k2 = randomUtil.next(0, nc - 1);
	// BOOST_LOG_TRIVIAL(trace)<< "Random k2 = " << k2;
	// removes node from cluster1 and inserts in cluster2
	return std::make_pair(node, k2);
}

std::pair< std::vector<int>, std::vector<int> > Perturbation::generateRandomMoveListCCProblem(int n,
		int nc, unsigned long numberOfMoves) {
	// BOOST_LOG_TRIVIAL(debug)<< "Generating perturbation of level " << numberOfMoves;
	std::pair< std::vector<int>, std::vector<int> > movementList;
	std::vector<int> nodeList;
	std::vector<int> k2List;
	for(int i = 0; i < numberOfMoves; i++) {
		std::pair<int, int> movement = randomMove1optCCProblem(n, nc);
		nodeList.push_back(movement.first);
		k2List.push_back(movement.second);
	}
	return std::make_pair(nodeList, k2List);
}

} /* namespace clusteringgraph */
