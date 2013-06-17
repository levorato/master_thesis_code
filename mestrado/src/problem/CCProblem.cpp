/*
 * CCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/CCProblem.h"
#include "../graph/include/Clustering.h"
#include <cmath>

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
Imbalance CCProblem::objectiveFunction(SignedGraph* g, Clustering* c) const {
	double positiveSum = 0, negativeSum = 0;
	int nc = c->getNumberOfClusters();
	int n = g->getN();

	// cout << "[CCProblem] Disparando calculo da funcao objetivo." << endl;
	// c->printClustering();
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
						break;
					}
				}
				// cout <<  g->getEdge(i, j) << " " << sameCluster << " ";
				if((g->isPositiveEdge(i, j)) && (not sameCluster)) {
					// cout <<  g->getEdge(i, j) << " is pos edge.\n";
					positiveSum += g->getEdge(i, j);
				} else if((g->isNegativeEdge(i, j)) && sameCluster) {
					// cout <<  g->getEdge(i, j) << " is neg edge.\n";
					negativeSum += fabs(g->getEdge(i, j));
				}
			}
		}
	}

	// cout << "Valor calculado: " << (positiveSum + negativeSum) << endl;
	return Imbalance(positiveSum, negativeSum);
}

} /* namespace problem */
