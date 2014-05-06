/*
 * Perturbation.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: mlevorato
 */

#include "include/Perturbation.h"
#include "../util/include/RandomUtil.h"

#include <boost/log/trivial.hpp>

namespace clusteringgraph {

using namespace util;
using namespace boost;

Perturbation::~Perturbation() {
	// TODO Auto-generated destructor stub
}

Clustering Perturbation::randomMove(SignedGraph* g, Clustering clustering, ClusteringProblem& p,
		unsigned long numberOfMoves) {

	BOOST_LOG_TRIVIAL(debug)<< "Generating perturbation of level " << numberOfMoves;
	Clustering c = clustering;
	for(int i = 0; i < numberOfMoves; i++) {
		// TODO avoid cyclic moves
		c = randomMove1opt(g, c, p);
	}
	return c;
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
	BoolArray cluster1 = cTemp.getCluster(k1);
	// (B) ==== Random node i from cluster 1 ====
	int idx_i = randomUtil.next(0, cTemp.getClusterSize(k1) - 1);
	for(int j = 0, count = 0; j < n; j++) {
		if(cluster1[j]) {
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

} /* namespace clusteringgraph */
