/*
 * CCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/CCProblem.h"

namespace problem {

CCProblem::CCProblem() {
	// TODO Auto-generated constructor stub

}

CCProblem::~CCProblem() {
	// TODO Auto-generated destructor stub
}

/**
 * The Imbalance of a partition P (I(P)) is defined as the total
 * weight of negative uncut arcs and positive cut arcs.
 */
int CCProblem::objectiveFunction(SignedGraph* g, Clustering* c) {
	int positiveSum = 0, negativeSum = 0;
	int nc = c->getNumberOfClusters();
	int n = c->getNumberOfNodes();
	// For each cluster i
	for(int i = 0; i < nc; i++) {
		BoolArray isInClusterI = *(c->getCluster(i));

		// calculates the sum of internal negative edges (within the same cluster)
		for(int a = 0; a < n; a++) {
			if(isInClusterI[a]) {
				for(int b = 0; b < n; b++) {
					if(isInClusterI[b]) {
						// nodes a and b are in the same cluster
						if(g->getEdge(a, b) < 0)
							negativeSum += g->getEdge(a, b);
					}
				}
			}
		}
		// calculates the sum of external positive edges (within different clusters)
		// For each cluster j
		for(int j = 0; j < nc; j++) {
			if(i == j)  continue;
			BoolArray isInClusterJ = *(c->getCluster(j));

			for(int a = 0; a < n; a++) {
				if(isInClusterI[a]) {
					for(int b = 0; b < n; b++) {
						if(isInClusterJ[b]) {
							// nodes a and b are in different clusters
							if(g->getEdge(a, b) > 0)
								positiveSum += g->getEdge(a, b);
						}
					}
				}
			}
		}
	}
	return (positiveSum + negativeSum);
}

} /* namespace problem */
