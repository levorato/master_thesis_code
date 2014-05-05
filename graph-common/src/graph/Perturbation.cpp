/*
 * Perturbation.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: mlevorato
 */

#include "include/Perturbation.h"

namespace clusteringgraph {


Perturbation::~Perturbation() {
	// TODO Auto-generated destructor stub
}

Clustering Perturbation::randomMove(SignedGraph* g, Clustering clustering, ClusteringProblem& p,
		unsigned long numberOfMoves) {

	Clustering c = clustering;
	for(int i = 0; i < numberOfMoves; i++) {
		// TODO avoid cyclic moves
		c = move1opt(g, c, p);
	}
	return c;
}

Clustering Perturbation::move1opt(SignedGraph* g, Clustering clustering, ClusteringProblem& p) {
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}

	// removes node i from cluster1 and inserts in cluster2
	Clustering cTemp = clustering;
	// shuffle 2 clusters (origin and destination) and a vertex i
	// ...
	int k1 = 0, k2 = 0, i = 0;
	BoolArray cluster2 = cTemp.getCluster(k2);

	//BOOST_LOG_TRIVIAL(trace) << "Option 1: Taking node " << i << " from cluster " << k1 << " to cluster " << k2;
	int nc = cTemp.getNumberOfClusters();
	cTemp.removeNodeFromCluster(*g, p, i, k1);
	// recalculates the number of clusters, as one of them may have been removed
	if((cTemp.getNumberOfClusters() < nc) && (k2 >= k1)) {
		// cluster k1 has been removed
		cTemp.addNodeToCluster(*g, p, i, k2 - 1);
	} else {
		cTemp.addNodeToCluster(*g, p, i, k2);
	}
	// numberOfTestedCombinations++;
	return cTemp;
}

} /* namespace clusteringgraph */
