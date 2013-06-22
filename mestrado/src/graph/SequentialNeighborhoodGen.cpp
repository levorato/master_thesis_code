/*
 * Neighborhood.cpp
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
#include "include/SequentialNeighborhoodGen.h"

namespace clusteringgraph {

SequentialNeighborhoodGenerator::SequentialNeighborhoodGenerator(int n) :
		NeighborhoodListGenerator(n) {

}

// TODO This method can be parallelized with MPI: neighborhood generation across several processors.
ClusteringPtr SequentialNeighborhoodGenerator::generateNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, const ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed) {
	int nc = clustering->getNumberOfClusters();
	int n = g->getN();
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	ClusteringPtr cBest = make_shared<Clustering>(*clustering, g->getN());
	// random number generators used in loop randomization
	boost::uniform_int<> distnc(0,nc-1);
	boost::uniform_int<> distN(0,n-1);
	boost::minstd_rand generator(randomSeed);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uninc(generator, distnc);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uniN(generator, distN);

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	if (l == 1) {  // 1-opt
		for (int k1 = uninc(), cont1 = 0; cont1 < nc; cont1++, k1 = (k1 + 1) % nc) {
			// cluster(k1)
			BoolArray cluster1 = clustering->getCluster(k1);
			// For each node i in cluster(k1)
			for (int i = uniN(), cont2 = 0; cont2 < n; cont2++, i = (i + 1) % n) {
				if (cluster1[i]) {
					// Option 1: node i is moved to another existing cluster k2
					for (int k2 = uninc(), cont3 = 0; cont3 < nc; cont3++, k2 = (k2 + 1) % nc) {
						if (k1 != k2) {
							// cluster(k2)
							// removes node i from cluster1 and inserts in cluster2
							ClusteringPtr cTemp = make_shared < Clustering
									> (*clustering, n);
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
					}
					// Option 2: node i is moved to a new cluster, alone
					// removes node i from cluster1 and inserts in newCluster
					// cout << "New clustering combination generated." << endl;
					ClusteringPtr cTemp = make_shared < Clustering
							> (*clustering, n);
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
		}
		// returns the best combination found in 1-opt
		return cBest;
	} else {  // 2-opt
		// cout << "Generating 2-opt neighborhood..." << endl;
		Imbalance bestImbalance = cBest->getImbalance();
		for (int k1 = uninc(), contk1 = 0; contk1 < nc; contk1++, k1 = (k1 + 1) % nc) {
			// cluster(k1)
			BoolArray cluster1 = clustering->getCluster(k1);
			for (int k2 = uninc(), contk2 = 0; contk2 < nc; contk2++, k2 = (k2 + 1) % nc) {
				// cluster(k2)
				BoolArray cluster2 = clustering->getCluster(k2);
				// For each node i in cluster(k1)
				for (int i = uniN(), conti = 0; conti < n; conti++, i = (i + 1) % n) {
					if (cluster1[i]) {
						// For each node j in cluster(k2)
						for (int j = uniN(), contj = 0; contj < n; contj++, j = (j + 1) % n) {
							if(j == i)  continue;
							if (cluster2[j]) {
								// Option 1: node i is moved to another existing cluster k3
								for (int k3 = uninc(), contk3 = 0; contk3 < nc; contk3++, k3 = (k3 + 1) % nc) {
									if (k1 != k3) {
										// cluster(k3)
										// Option 1: node i is moved to another existing cluster k3 and
										//           node j is moved to another existing cluster k4
										// cout << "Option 1" << endl;
										// TODO colocar laço aleatorio aqui! Verificar com o Yuri
										for (int k4 = k3 + 1; k4 < nc; k4++) {
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
								}
								// Option 3: node i is moved to a new cluster, alone, and
								//           node j is moved to another existing cluster k4
								// cout << "Option 3" << endl;
								for (int k4 = uninc(), contk4 = 0; contk4 < nc; contk4++, k4 = (k4 + 1) % nc) {
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
			}
		}
	}
	return cBest;
}

} /* namespace clusteringgraph */

namespace clusteringgraph {

SequentialNeighborhoodGenerator::~SequentialNeighborhoodGenerator() {
	// TODO Auto-generated destructor stub
}

} /* namespace clusteringgraph */
