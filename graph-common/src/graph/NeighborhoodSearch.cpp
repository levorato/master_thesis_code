#include <limits>
#include <boost/make_shared.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/timer/timer.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/optional.hpp>
#include <iomanip>
#include <iostream>
#include <cassert>
#include "include/NeighborhoodSearch.h"
#include "../util/include/MPIMessage.h"

using namespace boost::mpi;
using namespace std;
using namespace util;

// probes for mpi VND interrupt message (first improvement from other node) and returns if it exists
#define MPI_IPROBE_RETURN(ret_val) \
optional< mpi::status > stat = world.iprobe(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG); \
if (stat) { \
        InputMessageParallelVND imsg; \
        mpi::status stat = world.recv(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG, imsg); \
    return (ret_val); }


namespace clusteringgraph {

Clustering NeighborhoodSearch::search1opt(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialClusterIndex,
        		unsigned long finalClusterIndex, bool firstImprovement, unsigned long k) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned long n = g->getN();
	unsigned long numberOfClustersInInterval = finalClusterIndex - initialClusterIndex + 1;
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	boost::uniform_int<> distnc(0, totalNumberOfClusters - 1);
	boost::uniform_int<> distN(0, n - 1);
	boost::minstd_rand generator(randomSeed);
	generator.seed(boost::random::random_device()());
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uninc(generator, distnc);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uniN(generator, distN);
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;

	for (unsigned long k1 = uninc(), cont1 = 0; cont1 < numberOfClustersInInterval; cont1++) {
		// cluster(k1)
		// RCC Problem: must have exactly k clusters -> Cannot remove node from standalone cluster
		unsigned long s = clustering->getClusterSize(k1);
		if((problem.getType() == ClusteringProblem::RCC_PROBLEM) && (s == 1)) {
			continue;
		}
		BoolArray cluster1 = clustering->getCluster(k1);
		// For each node i in cluster(k1)
		for (unsigned long i = uniN(), cont2 = 0; cont2 < n; cont2++, i = (i + 1) % n) {
			if (cluster1[i]) {
				// Option 1: node i is moved to another existing cluster k2
				for (unsigned long k2 = uninc(), cont3 = 0; cont3 < totalNumberOfClusters; cont3++) {
					if (k1 != k2) {
						// cluster(k2)
						// removes node i from cluster1 and inserts in cluster2
						Clustering cTemp = *clustering;
						BoolArray cluster2 = cTemp.getCluster(k2);

						//BOOST_LOG_TRIVIAL(trace) << "Option 1: Taking node " << i << " from cluster " << k1 << " to cluster " << k2;
						int nc = cTemp.getNumberOfClusters();
						cTemp.removeNodeFromCluster(*g, problem, i, k1);
						// recalculates the number of clusters, as one of them may have been removed
						if((cTemp.getNumberOfClusters() < nc) && (k2 >= k1)) {
							// cluster k1 has been removed
							cTemp.addNodeToCluster(*g, problem, i, k2 - 1);
						} else {
							cTemp.addNodeToCluster(*g, problem, i, k2);
						}
						numberOfTestedCombinations++;
						// cTemp->printClustering();
						Imbalance newImbalance = cTemp.getImbalance();
						Imbalance bestImbalance = cBest.getImbalance();
						if (newImbalance < bestImbalance) {
							//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 1-neighborhood: " << setprecision(2) << newImbalance.getValue();
							// First improvement for 1-opt neighborhood
							cBest = cTemp;
							if(firstImprovement) {
								return cBest;
							}
						}
					}
					// return if time limit is exceeded
					boost::timer::cpu_times end_time = timer.elapsed();
					double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
					// std::cout << timeSpentSoFar + localTimeSpent << endl;
					if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
					// loop increment rule
					k2++;
					if(k2 >= totalNumberOfClusters) {
						k2 = 0;
					}
				}
				// Option 2: node i is moved to a new cluster, alone
				// removes node i from cluster1 and inserts in newCluster
				// cout << "New clustering combination generated." << endl;
				if((problem.getType() == ClusteringProblem::CC_PROBLEM) && (s > 1)) {
					// this code is not executed for RCC Problem, as it must have exactly k clusters
					// this code is not executed if node i is to be moved from a standalone cluster to another standalone cluster (s == 1)
					Clustering cTemp = *clustering;
					//BOOST_LOG_TRIVIAL(trace) << "Option 2: Taking node " << i << " from " << k1 << " to new cluster.";
					cTemp.removeNodeFromCluster(*g, problem, i, k1);
					BoolArray cluster2 = cTemp.addCluster(*g, problem, i);
					numberOfTestedCombinations++;
					// cTemp->printClustering();
					Imbalance newImbalance = cTemp.getImbalance();
					Imbalance bestImbalance = cBest.getImbalance();
					if (newImbalance < bestImbalance) {
						//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 1-neighborhood: " << setprecision(2) << newImbalance.getValue();
						// First improvement for 1-opt neighborhood
						cBest = cTemp;
						if(firstImprovement) {
							return cBest;
						}
					}
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
}

Clustering NeighborhoodSearch::search2opt(SignedGraph* g,
        Clustering* clustering, ClusteringProblem* problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex, bool firstImprovement, unsigned long k) {

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	//BOOST_LOG_TRIVIAL(trace) << "Generating 2-opt neighborhood...";
	mpi::communicator world;
	unsigned long n = g->getN();
	unsigned long numberOfClustersInInterval = finalClusterIndex - initialClusterIndex + 1;
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	boost::uniform_int<> distn1(initialClusterIndex, finalClusterIndex);
	boost::uniform_int<> distnc(0, totalNumberOfClusters - 1);
	boost::uniform_int<> distN(0, n - 1);
	boost::minstd_rand generator(randomSeed);
	generator.seed(boost::random::random_device()());
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > unin1(generator, distn1);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uninc(generator, distnc);
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uniN(generator, distN);
	// stores the best clustering combination generated (minimum imbalance)
	Clustering cBest = *clustering;
	Clustering cTemp = cBest;

	Imbalance bestImbalance = cBest.getImbalance();
	for (unsigned long k1 = uninc(), contk1 = 0; contk1 < numberOfClustersInInterval; contk1++) {
		// cluster(k1)
		unsigned long s1 = clustering->getClusterSize(k1);
		if((problem->getType() == ClusteringProblem::RCC_PROBLEM) && (s1 == 1)) {
			// no clusters can be removed if RCC Problem
			continue;
		}
		BoolArray cluster1 = clustering->getCluster(k1);
		// if(problem->getType() == ClusteringProblem::RCC_PROBLEM) assert(cluster1.size() > 1);
		for (unsigned long k2 = uninc(), contk2 = 0; contk2 < totalNumberOfClusters; contk2++) {
			// cluster(k2)
			unsigned long s2 = clustering->getClusterSize(k2);
			if((problem->getType() == ClusteringProblem::RCC_PROBLEM) && (s2 == 1)) {
				// no clusters can be removed if RCC Problem
				continue;
			}
			if(k1 == k2) {  // vertices i and j are being removed from the same cluster k1
				if((problem->getType() == ClusteringProblem::RCC_PROBLEM) && (s2 <= 2)) {
					// no clusters can be removed if RCC Problem
					continue;
				} else if((problem->getType() == ClusteringProblem::CC_PROBLEM) && (s2 < 2)) {
					// cluster k1 must have at least 2 elements
					continue;
				}
			}
			BoolArray cluster2 = clustering->getCluster(k2);
			// if(problem->getType() == ClusteringProblem::RCC_PROBLEM) assert(cluster2.size() > 1);
			// For each node i in cluster(k1)
			for (unsigned long i = uniN(), conti = 0; conti < n; conti++, i = (i + 1) % n) {
				if (cluster1[i]) {
					// For each node j in cluster(k2)
					for (unsigned long j = uniN(), contj = 0; contj < n; contj++, j = (j + 1) % n) {
						if(j == i)  continue;
						if (cluster2[j]) {
							// Option 1: node i is moved to another existing cluster k3
							for (unsigned long k3 = uninc(), contk3 = 0; contk3 < totalNumberOfClusters; contk3++) {
								MPI_IPROBE_RETURN(cBest)

								if (k1 != k3) {
									// cluster(k3)
									// Option 1: node i is moved to another existing cluster k3 and
									//           node j is moved to another existing cluster k4
									for (unsigned long k4 = k3 + 1; k4 < totalNumberOfClusters; k4++) {
										if (k2 != k4) {
											// cluster(k4)
											// BOOST_LOG_TRIVIAL(trace) << "Option 1: node " << i << " into cluster " << k3 << " and node " << j << " into cluster " << k4;
											cTemp = *clustering;
											bool success = NeighborhoodSearch::process2optCombination(*g,
															cTemp, problem, k1,
															k2, k3, k4, n,
															i, j);
											// cTemp->printClustering();
											if(not success) {
												continue;
											}
											if((problem->getType() == ClusteringProblem::RCC_PROBLEM) && (cTemp.getNumberOfClusters() > clustering->getNumberOfClusters())) {
												BOOST_LOG_TRIVIAL(error) << "Cluster with more than k clusters generated!";
											}

											Imbalance newImbalance = cTemp.getImbalance();
											if (newImbalance < bestImbalance) {
												//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 2-neighborhood.";
												cBest = cTemp;
												bestImbalance = cBest.getImbalance();
												if(firstImprovement) {
													return cBest;
												}
											}
											// return if time limit is exceeded
											boost::timer::cpu_times end_time = timer.elapsed();
											double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
											if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
										}
									}
									// Option 2: node j is moved to a new cluster, alone
									if((problem->getType() == ClusteringProblem::CC_PROBLEM) && (s2 > 1)) {
										// this code is not executed in RCC Problem, as it must have exactly k clusters
										// this code is not executed if node j is to be moved from a standalone cluster to another standalone cluster
										//BOOST_LOG_TRIVIAL(trace) << "Option 2: node " << i << " into cluster " << k3 << " and node " << j << " into new cluster.";
										cTemp = *clustering;
										NeighborhoodSearch::process2optCombination(*g,
														cTemp, problem, k1, k2, k3,
														Clustering::NEW_CLUSTER,
														n, i, j);
										// cTemp->printClustering();
										Imbalance newImbalance = cTemp.getImbalance();
										if (newImbalance < bestImbalance) {
											//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 2-neighborhood.";
											cBest = cTemp;
											bestImbalance = cBest.getImbalance();
											if(firstImprovement) {
												return cBest;
											}
										}
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
							if((problem->getType() == ClusteringProblem::CC_PROBLEM) && (s1 > 1)) {
								// this code is not executed in RCC Problem, as it must have exactly k clusters
								// this code is not executed if node i is to be moved from a standalone cluster to another standalone cluster
								for (unsigned long k4 = uninc(), contk4 = 0; contk4 < numberOfClustersInInterval; contk4++) {
									if (k2 != k4) {
										// cluster(k4)
										//BOOST_LOG_TRIVIAL(trace) << "Option 3: " << i << " into new cluster and " << j << " into cluster " << k4;
										cTemp = *clustering;
										NeighborhoodSearch::process2optCombination(*g,
														cTemp, problem, k1, k2,
														Clustering::NEW_CLUSTER,
														k4, n, i, j);
										//cTemp->printClustering();
										Imbalance newImbalance = cTemp.getImbalance();
										if (newImbalance < bestImbalance) {
											//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 2-neighborhood.";
											cBest = cTemp;
											bestImbalance = cBest.getImbalance();
											if(firstImprovement) {
												return cBest;
											}
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
							}
							// Option 4: nodes i and j are moved to 2 new clusters, and left alone
							if((problem->getType() == ClusteringProblem::CC_PROBLEM) && (s1 > 1) && (s2 > 1)) {
								// this code is not executed in RCC Problem, as it must have exactly k clusters
								// this code is not executed if node i and/or j is to be moved from a standalone cluster to another standalone cluster
								//BOOST_LOG_TRIVIAL(trace) << "Option 4: nodes " << i << " and " << j << " into new clusters.";
								cTemp = *clustering;
								NeighborhoodSearch::process2optCombination(*g,
										cTemp, problem, k1, k2,
										Clustering::NEW_CLUSTER,
										Clustering::NEW_CLUSTER, n, i, j);
								// cTemp->printClustering();
								Imbalance newImbalance = cTemp.getImbalance();
								if (newImbalance < bestImbalance) {
									//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 2-neighborhood.";
									cBest = cTemp;
									bestImbalance = cBest.getImbalance();
									if(firstImprovement) {
										return cBest;
									}
								}
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
	return cBest;
}

}


