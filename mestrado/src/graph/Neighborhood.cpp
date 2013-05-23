/*
 * Neighborhood.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/Neighborhood.h"
#include <limits>
#include <boost/make_shared.hpp>


namespace clusteringgraph {

NeighborhoodList::NeighborhoodList(int n) :
		numberOfNodes(n) {

}

ClusteringPtr NeighborhoodList::process2optCombination(Clustering *clustering, int k1, int k2, int k3,
		int k4, int n, int i, int j) {

	ClusteringPtr cTemp = make_shared<Clustering>(*clustering, n);
	// removes node i from cluster1 and inserts in cluster3
	// TODO check if the removal of node i destroys cluster1
	cTemp->removeNodeFromCluster(i, k1);
	if(k3 >= 0) {
		// inserts i in existing cluster k3
		cTemp->addNodeToCluster(i, k3);
	} else {
		// inserts i in a new cluster (alone)
		int nodeArray[1] = {i};
		cTemp->addCluster(nodeArray, 1);
	}
	// removes node j from cluster2 and inserts in cluster4
	cTemp->removeNodeFromCluster(j, k2);
	if(k4 >= 0) {
		// inserts j in existing cluster k4
		cTemp->addNodeToCluster(j, k4);
	} else {
		// inserts j in a new cluster (alone)
		int nodeArray[1] = {j};
		cTemp->addCluster(nodeArray, 1);
	}

	return cTemp;
}

ClusteringPtr NeighborhoodList::generateNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem) {
	int nc = clustering->getNumberOfClusters();
	int n = clustering->getNumberOfNodes();
	float bestValue = problem.objectiveFunction(g, clustering);

	if(l == 1) {  // 1-opt
		for(int k1 = 0; k1 < nc; k1++) {
			// cluster(k1)
			BoolArray* cluster1 = clustering->getCluster(k1);
			// For each node i in cluster(k1)
			for(int i = 0; i < n; i++) {
				if((*cluster1)[i]) {
					// Option 1: node i is moved to another existing cluster k2
					for(int k2 = 0; k2 < nc; k2++) {
						if(k1 != k2) {
							// cluster(k2)
							// removes node i from cluster1 and inserts in cluster2
							cout << "New clustering combination generated." << endl;
							ClusteringPtr cTemp = make_shared<Clustering>(*clustering, n);
							// TODO check if the removal of node i destroys cluster1
							cTemp->removeNodeFromCluster(i, k1);
							cTemp->addNodeToCluster(i, k2);
							// cTemp->printClustering();
							float objective = problem.objectiveFunction(g, cTemp.get());
							if(objective < bestValue) {
								cout << "Better solution found in 1-neighborhood." << endl;
								return cTemp;
							}
						}
					}
					// Option 2: node i is moved to a new cluster, alone
					// removes node i from cluster1 and inserts in newCluster
					cout << "New clustering combination generated." << endl;
					ClusteringPtr cTemp = make_shared<Clustering>(*clustering, n);
					cTemp->removeNodeFromCluster(i, k1);
					int nodeArray[1] = {i};
					cTemp->addCluster(nodeArray, 1);
					// cTemp->printClustering();
					float objective = problem.objectiveFunction(g, cTemp.get());
					if(objective < bestValue) {
						cout << "Better solution found in 1-neighborhood." << endl;
						return cTemp;
					}
				}
			}
		}
	} else {  // 2-opt
		for(int k1 = 0; k1 < nc; k1++) {
			// cluster(k1)
			BoolArray* cluster1 = clustering->getCluster(k1);
			for(int k2 = k1 + 1; k2 < nc; k2++) {
				// cluster(k2)
				BoolArray* cluster2 = clustering->getCluster(k2);
				// For each node i in cluster(k1)
				for(int i = 0; i < n; i++) {
					if((*cluster1)[i]) {
						// For each node j in cluster(k2)
						for(int j = 0; j < n; j++) {
							if((*cluster2)[j]) {
								// Option 1: node i is moved to another existing cluster k3
								for(int k3 = 0; k3 < nc; k3++) {
									if(k1 != k3) {
										// cluster(k3)
										// Option 1: node i is moved to another existing cluster k3 and
										//           node j is moved to another existing cluster k4
										for(int k4 = k3 + 1; k4 < nc; k4++) {
											if(k2 != k4) {
												// cluster(k4)
												ClusteringPtr cTemp = process2optCombination(clustering, k1, k2, k3, k4, n, i, j);
												// cTemp->printClustering();
												float objective = problem.objectiveFunction(g, cTemp.get());
												if(objective < bestValue) {
													cout << "Better solution found in 2-neighborhood." << endl;
													return cTemp;
												}
											}
										}
										// Option 2: node j is moved to a new cluster, alone
										ClusteringPtr cTemp = process2optCombination(clustering, k1, k2, k3, Clustering::NEW_CLUSTER,
												n, i, j);
										// cTemp->printClustering();
										float objective = problem.objectiveFunction(g, cTemp.get());
										if(objective < bestValue) {
											cout << "Better solution found in 2-neighborhood." << endl;
											return cTemp;
										}
									}
								}
								// Option 3: node i is moved to a new cluster, alone, and
								//           node j is moved to another existing cluster k4
								for(int k4 = 0; k4 < nc; k4++) {
									if(k2 != k4) {
										// cluster(k4)
										ClusteringPtr cTemp = process2optCombination(clustering, k1, k2, Clustering::NEW_CLUSTER, k4,
												n, i, j);
										// cTemp->printClustering();
										float objective = problem.objectiveFunction(g, cTemp.get());
										if(objective < bestValue) {
											cout << "Better solution found in 2-neighborhood." << endl;
											return cTemp;
										}
									}
								}
								// Option 4: nodes i and j are moved to new clusters, and left alone
								ClusteringPtr cTemp = process2optCombination(clustering, k1, k2, Clustering::NEW_CLUSTER, Clustering::NEW_CLUSTER,
										n, i, j);
								// cTemp->printClustering();
								float objective = problem.objectiveFunction(g, cTemp.get());
								if(objective < bestValue) {
									cout << "Better solution found in 2-neighborhood." << endl;
									return cTemp;
								}
							}
						}
					}
				}
			}
		}
	}
	return ClusteringPtr();
}

} /* namespace clusteringgraph */
