/*
 * CCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/CCProblem.h"
#include "../graph/include/Clustering.h"

namespace problem {

using namespace std;
using namespace clusteringgraph;

CCProblem::CCProblem() {
	// TODO Auto-generated constructor stub

}

CCProblem::~CCProblem() {
	// TODO Auto-generated destructor stub
}

int CCProblem::getType() const {
	return ClusteringProblem::CC_PROBLEM;
}

/**
 * The Imbalance of a partition P (I(P)) is defined as the total
 * weight of negative uncut arcs and positive cut arcs.
 */
float CCProblem::objectiveFunction(SignedGraph* g, Clustering* c) const {
	float positiveSum = 0, negativeSum = 0;
	int nc = c->getNumberOfClusters();
	int n = g->getN();

	// cout << "[CCProblem] Disparando calculo da funcao objetivo." << endl;
	// c->printClustering();
	/*
	float total = 0;
	for(int k = 0; k < nc; k++) {
		float soma = 0, soma2 = 0;
		vector<int> cluster;
    	for(int i = 0; i < n; i++) {
    		if(c->getCluster(k)[i]) {
    			cluster.push_back(i);
    		}
    	}
    	for(int x = 0; x < cluster.size(); x++) {
    		int i = cluster.at(x);
    		for(int y = 0; y < cluster.size(); y++) {
    			int j = cluster.at(y);
    			if(g->getEdge(i, j) < 0) {
    				soma += abs(g->getEdge(i, j));
    			}
    		}
    		for(int a = 0; a < n; a++) {
    			if((!c->getCluster(k)[a]) && (g->getEdge(i, a) > 0)) {
    				soma2 += g->getEdge(i, a);
    			}
    			if((!c->getCluster(k)[a]) && (g->getEdge(a, i) > 0)) {
					soma2 += g->getEdge(a, i);
				}
    		}
    	}
    	// cout << "Internal negative sum of cluster " << k << ": " << soma << endl;
    	// cout << "External positive sum of cluster " << k << ": " << soma2 << endl;
    	total += soma + soma2;
    }
    // cout << "Calculated I(P) = " << total << endl;
    //c->printClustering();
	return total;
	*/
	// For each vertex i
	for(int i = 0; i < n; i++) {
		// For each vertex j
		for(int j = 0; j < n; j++) {
			if(i != j) {
				bool sameCluster = false;
				for(int k = 0; k < nc; k++) {
					BoolArray cluster = c->getCluster(k);
					if(cluster[i] && cluster[j]) {
						sameCluster = true;
					}
				}
				if((g->getEdge(i, j) > 0) && (not sameCluster)) {
					positiveSum += g->getEdge(i, j);
				} else if((g->getEdge(i, j) < 0) && sameCluster) {
					negativeSum += abs(g->getEdge(i, j));
				}
			}
		}
	}
	// cout << "Valor calculado: " << (positiveSum + negativeSum) << endl;
	return (positiveSum + negativeSum);
}

} /* namespace problem */
