/*
 * SequentialNeighborhoodSearch.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/Neighborhood.h"
#include <limits>
#include <boost/make_shared.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <cassert>
#include "include/SequentialNeighborhoodSearch.h"

namespace clusteringgraph {

SequentialNeighborhoodSearch::SequentialNeighborhoodSearch() {

}

ClusteringPtr SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank) {
	unsigned long nc = clustering->getNumberOfClusters();
	return SequentialNeighborhoodSearch::searchNeighborhood(l, g, clustering, problem,
			timeSpentSoFar, timeLimit, randomSeed, 0, nc - 1);
}

ClusteringPtr SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex) {

	unsigned long n = g->getN();
	assert(initialClusterIndex < clustering->getNumberOfClusters());
	assert(finalClusterIndex < clustering->getNumberOfClusters());
	unsigned long numberOfClustersInInterval = finalClusterIndex - initialClusterIndex;
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	ClusteringPtr cBest = make_shared<Clustering>(*clustering);
	// random number generators used in loop randomization
	boost::uniform_int<> distnc(initialClusterIndex, finalClusterIndex);
	boost::uniform_int<> distN(0, n-1);
	boost::minstd_rand generator(randomSeed);
	generator.seed(boost::random::random_device()());
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uninc(generator, distnc);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uniN(generator, distN);

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	if (l == 1) {  // 1-opt
		// TODO change code to best improvement in 1-opt
		for (unsigned long k1 = uninc(), cont1 = 0; cont1 < numberOfClustersInInterval; cont1++) {
			// cluster(k1)
			BoolArray cluster1 = clustering->getCluster(k1);
			// For each node i in cluster(k1)
			for (unsigned long i = uniN(), cont2 = 0; cont2 < n; cont2++, i = (i + 1) % n) {
				if (cluster1[i]) {
					// Option 1: node i is moved to another existing cluster k2
					for (unsigned long k2 = uninc(), cont3 = 0; cont3 < numberOfClustersInInterval; cont3++) {
						if (k1 != k2) {
							// cluster(k2)
							// removes node i from cluster1 and inserts in cluster2
							ClusteringPtr cTemp = make_shared < Clustering
									> (*clustering);
							BoolArray cluster2 = cTemp->getCluster(k2);

							// cout << "Taking node " << i << " from " << k1 << " to " << k2 << endl;
							cTemp->addNodeToCluster(*g, i, k2);
							cTemp->removeNodeFromCluster(*g, i, k1);

							// cTemp->printClustering();
							Imbalance newImbalance = cTemp->getImbalance();
							Imbalance bestImbalance = cBest->getImbalance();
							if (newImbalance < bestImbalance) {
								// cout << "Better solution found in 1-neighborhood: " << setprecision(2) << objective << "\n";
								// First improvement for 1-opt neighborhood
								cBest.reset();
								cBest = cTemp;
								return cBest;
							}
						}
						// return if time limit is exceeded
						boost::timer::cpu_times end_time = timer.elapsed();
						double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
						if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
						// increment rule
						k2++;
						if(k2 >= totalNumberOfClusters) {
							k2 = 0;
						}
					}
					// Option 2: node i is moved to a new cluster, alone
					// removes node i from cluster1 and inserts in newCluster
					// cout << "New clustering combination generated." << endl;
					ClusteringPtr cTemp = make_shared < Clustering
							> (*clustering);
					// cout << "Taking node " << i << " from " << k1 << " to new cluster." << endl;
					cTemp->removeNodeFromCluster(*g, i, k1);
					BoolArray cluster2 = cTemp->addCluster(*g, i);
					// cTemp->printClustering();
					Imbalance newImbalance = cTemp->getImbalance();
					Imbalance bestImbalance = cBest->getImbalance();
					if (newImbalance < bestImbalance) {
						// cout << "Better solution found in 1-neighborhood: " << setprecision(2) << objective << "\n";
						// First improvement for 1-opt neighborhood
						cBest.reset();
						cBest = cTemp;
						return cBest;
					}
				}
			}
			// increment rule
			k1++;
			if(k1 > finalClusterIndex) {
				k1 = initialClusterIndex;
			}
		}
		// returns the best combination found in 1-opt
		return cBest;
	} else {  // 2-opt
		// TODO insert logic to deal with first improvement in parallel VNS
		// cout << "Generating 2-opt neighborhood..." << endl;
		Imbalance bestImbalance = cBest->getImbalance();
		for (unsigned long k1 = uninc(), contk1 = 0; contk1 < numberOfClustersInInterval; contk1++) {
			// cluster(k1)
			BoolArray cluster1 = clustering->getCluster(k1);
			for (unsigned long k2 = uninc(), contk2 = 0; contk2 < numberOfClustersInInterval; contk2++) {
				// cluster(k2)
				BoolArray cluster2 = clustering->getCluster(k2);
				// For each node i in cluster(k1)
				for (unsigned long i = uniN(), conti = 0; conti < n; conti++, i = (i + 1) % n) {
					if (cluster1[i]) {
						// For each node j in cluster(k2)
						for (unsigned long j = uniN(), contj = 0; contj < n; contj++, j = (j + 1) % n) {
							if(j == i)  continue;
							if (cluster2[j]) {
								// Option 1: node i is moved to another existing cluster k3
								for (unsigned long k3 = uninc(), contk3 = 0; contk3 < numberOfClustersInInterval; contk3++) {
									if (k1 != k3) {
										// cluster(k3)
										// Option 1: node i is moved to another existing cluster k3 and
										//           node j is moved to another existing cluster k4
										// cout << "Option 1" << endl;
										// TODO colocar laço aleatorio aqui! Verificar com o Yuri
										for (unsigned long k4 = k3 + 1; k4 < totalNumberOfClusters; k4++) {
											if (k2 != k4) {
												// cluster(k4)
												ClusteringPtr cTemp =
														process2optCombination(*g,
																clustering, k1,
																k2, k3, k4, n,
																i, j);
												// cTemp->printClustering();
												Imbalance newImbalance = cTemp->getImbalance();
												if (newImbalance < bestImbalance) {
													// cout << "Better solution found in 2-neighborhood." << endl;
													return cTemp;
												}
												// return if time limit is exceeded
												boost::timer::cpu_times end_time = timer.elapsed();
												double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
												if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
											}
										}
										// Option 2: node j is moved to a new cluster, alone
										// cout << "Option 2" << endl;
										ClusteringPtr cTemp =
												process2optCombination(*g,
														clustering, k1, k2, k3,
														Clustering::NEW_CLUSTER,
														n, i, j);
										// cTemp->printClustering();
										Imbalance newImbalance = cTemp->getImbalance();
										if (newImbalance < bestImbalance) {
											// cout << "Better solution found in 2-neighborhood." << endl;
											return cTemp;
										}
									}
									// return if time limit is exceeded
									boost::timer::cpu_times end_time = timer.elapsed();
									double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
									if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
									// increment rule
									k3++;
									if(k3 >= totalNumberOfClusters) {
										k3 = 0;
									}
								}
								// Option 3: node i is moved to a new cluster, alone, and
								//           node j is moved to another existing cluster k4
								// cout << "Option 3" << endl;
								for (unsigned long k4 = uninc(), contk4 = 0; contk4 < numberOfClustersInInterval; contk4++) {
									if (k2 != k4) {
										// cluster(k4)
										ClusteringPtr cTemp =
												process2optCombination(*g,
														clustering, k1, k2,
														Clustering::NEW_CLUSTER,
														k4, n, i, j);
										// cTemp->printClustering();
										Imbalance newImbalance = cTemp->getImbalance();
										if (newImbalance < bestImbalance) {
											// cout << "Better solution found in 2-neighborhood." << endl;
											return cTemp;
										}
									}
									// return if time limit is exceeded
									boost::timer::cpu_times end_time = timer.elapsed();
									double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
									if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
									// increment rule
									k4++;
									if(k4 >= totalNumberOfClusters) {
										k4 = 0;
									}
								}
								// Option 4: nodes i and j are moved to new clusters, and left alone
								// cout << "Option 4" << endl;
								ClusteringPtr cTemp = process2optCombination(*g,
										clustering, k1, k2,
										Clustering::NEW_CLUSTER,
										Clustering::NEW_CLUSTER, n, i, j);
								// cTemp->printClustering();
								Imbalance newImbalance = cTemp->getImbalance();
								if (newImbalance < bestImbalance) {
									// cout << "Better solution found in 2-neighborhood." << endl;
									return cTemp;
								}
								// return if time limit is exceeded
								boost::timer::cpu_times end_time = timer.elapsed();
								double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
								if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
							}
						}
					}
				}
				// increment rule
				k2++;
				if(k2 >= totalNumberOfClusters) {
					k2 = 0;
				}
			}
			// increment rule
			k1++;
			if(k1 > finalClusterIndex) {
				k1 = initialClusterIndex;
			}
		}
	}
	return cBest;
}

} /* namespace clusteringgraph */

namespace clusteringgraph {

SequentialNeighborhoodSearch::~SequentialNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

} /* namespace clusteringgraph */
