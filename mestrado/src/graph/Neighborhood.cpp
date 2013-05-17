/*
 * Neighborhood.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/Neighborhood.h"
#include <limits>


namespace clusteringgraph {

NeighborhoodList::NeighborhoodList(Clustering *c, int n) :
		numberOfNodes(n), clusteringPtr(c) {

}

// TODO: Implementar de acordo com o especificado pelo Yuri
// para 1-opt e 2-opt
Clustering* NeighborhoodList::generateNeighborhood(int l, SignedGraph* g, ClusteringProblem* problem) {
	int nc = this->clusteringPtr->getNumberOfClusters();
	int n = this->clusteringPtr->getNumberOfNodes();
	float originalQuality = problem->objectiveFunction(g, this->clusteringPtr.get());

	if(l == 1) {  // 1-opt
		for(int k1 = 0; k1 < nc; k1++) {
			// cluster(k1)
			const BoolArray cluster1 = this->clusteringPtr->getCluster(k1);
			// For each node i in cluster(k1)
			for(int i = 0; i < n; i++) {
				if(cluster1[i]) {
					for(int k2 = 0; k2 < nc; k2++) {
						if(k1 != k2) {
							// cluster(k2)
							const BoolArray cluster2 = this->clusteringPtr->getCluster(k2);
							// removes node i from cluster1 and inserts in cluster2
							Clustering *cTemp = new Clustering(this->clusteringPtr.get(), n);
							// TODO check if the removal of node i destroys cluster1
							cTemp->removeNodeFromCluster(i, k1);
							// TODO program an alternative: node i becomes a new cluster, all by itself
							cTemp->addNodeToCluster(i, k2);
							float quality = problem->objectiveFunction(g, cTemp);
							if(quality > originalQuality) {
								return cTemp;
							}
						}
					}
				}
			}
		}
	} else {  // 2-opt
		for(int k1 = 0; k1 < nc; k1++) {
			// cluster(k1)
			const BoolArray cluster1 = this->clusteringPtr->getCluster(k1);
			for(int k2 = k1 + 1; k2 < nc; k2++) {
				// cluster(k2)
				const BoolArray cluster2 = this->clusteringPtr->getCluster(k2);
				// For each node i in cluster(k1)
				for(int i = 0; i < n; i++) {
					if(cluster1[i]) {
						// For each node j in cluster(k2)
						for(int j = 0; j < n; j++) {
							if(cluster2[j]) {
								for(int k3 = 0; k3 < nc; k3++) {
									if(k1 != k3) {
										// cluster(k3)
										const BoolArray cluster3 = this->clusteringPtr->getCluster(k3);
										for(int k4 = k3 + 1; k4 < nc; k4++) {
											if(k2 != k4) {
												// cluster(k4)
												const BoolArray cluster4 = this->clusteringPtr->getCluster(k4);
												Clustering *cTemp = new Clustering(this->clusteringPtr.get(), n);
												// removes node i from cluster1 and inserts in cluster3
												// TODO check if the removal of node i destroys cluster1
												cTemp->removeNodeFromCluster(i, k1);
												// TODO program an alternative: node i becomes a new cluster, all by itself
												cTemp->addNodeToCluster(i, k3);
												// removes node j from cluster2 and inserts in cluster4
												cTemp->removeNodeFromCluster(j, k2);
												cTemp->addNodeToCluster(j, k4);

												float quality = problem->objectiveFunction(g, cTemp);
												if(quality > originalQuality) {
													return cTemp;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return this->clusteringPtr.get();
}

} /* namespace clusteringgraph */
