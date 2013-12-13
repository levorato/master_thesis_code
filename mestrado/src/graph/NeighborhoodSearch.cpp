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
#include "../resolution/grasp/include/ParallelGrasp.h"
#include "../util/include/MPIMessage.h"

using namespace boost::mpi;
using namespace resolution::grasp;
using namespace std;
using namespace util;

// probes for mpi VNS interrupt message (first improvement from other node) and returns if it exists
#define MPI_IPROBE_RETURN(ret_val) \
optional< mpi::status > stat = world.iprobe(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VNS_TAG); \
if (stat) { \
	InputMessageParallelVNS imsg; \
	mpi::status stat = world.recv(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VNS_TAG, imsg); \
	BOOST_LOG_TRIVIAL(trace) << "VNS interrupt message received."; \
    return ret_val; }


namespace clusteringgraph {

ClusteringPtr NeighborhoodSearch::search1opt(SignedGraph* g,
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
	ClusteringPtr cBest = make_shared<Clustering>(*clustering);

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
						cTemp->addNodeToCluster(*g, problem, i, k2);
						cTemp->removeNodeFromCluster(*g, problem, i, k1);
						// Full recalculation of objective value if RCC Problem
						if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
							Imbalance imb = problem.objectiveFunction(*g, *cTemp);
							if(imb.getValue() != cTemp->getImbalance().getValue()) {
								BOOST_LOG_TRIVIAL(error) << "RCC obj function and delta do not match!";
							} else {
								BOOST_LOG_TRIVIAL(error) << "RCC obj function and delta MATCH!";
							}
							cTemp->setImbalance(imb);
						}
						numberOfTestedCombinations++;

						// cTemp->printClustering();
						Imbalance newImbalance = cTemp->getImbalance();
						Imbalance bestImbalance = cBest->getImbalance();
						if (newImbalance < bestImbalance) {
							// cout << "Better solution found in 1-neighborhood: " << setprecision(2) << objective << "\n";
							// First improvement for 1-opt neighborhood
							cBest.reset();
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
				if((problem.getType() == ClusteringProblem::CC_PROBLEM) ||
						( (problem.getType() == ClusteringProblem::RCC_PROBLEM) && (clustering->getNumberOfClusters() < k) )) {
					ClusteringPtr cTemp = make_shared < Clustering
							> (*clustering);
					// cout << "Taking node " << i << " from " << k1 << " to new cluster." << endl;
					cTemp->removeNodeFromCluster(*g, problem, i, k1);
					BoolArray cluster2 = cTemp->addCluster(*g, problem, i);
					// Full recalculation of objective value if RCC Problem
					if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
						Imbalance imb = problem.objectiveFunction(*g, *cTemp);
						if(imb.getValue() != cTemp->getImbalance().getValue()) {
							BOOST_LOG_TRIVIAL(error) << "RCC obj function and delta do not match!";
						} else {
							BOOST_LOG_TRIVIAL(error) << "RCC obj function and delta MATCH!";
						}
						cTemp->setImbalance(imb);
					}
					numberOfTestedCombinations++;
					// cTemp->printClustering();
					Imbalance newImbalance = cTemp->getImbalance();
					Imbalance bestImbalance = cBest->getImbalance();
					if (newImbalance < bestImbalance) {
						// cout << "Better solution found in 1-neighborhood: " << setprecision(2) << objective << "\n";
						// First improvement for 1-opt neighborhood
						cBest.reset();
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

ClusteringPtr NeighborhoodSearch::search2opt(SignedGraph* g,
        Clustering* clustering, ClusteringProblem& problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex, bool firstImprovement, unsigned long k) {

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// cout << "Generating 2-opt neighborhood..." << endl;
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
	ClusteringPtr cBest = make_shared<Clustering>(*clustering);

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
								MPI_IPROBE_RETURN(cBest)

								if (k1 != k3) {
									// cluster(k3)
									// Option 1: node i is moved to another existing cluster k3 and
									//           node j is moved to another existing cluster k4
									// cout << "Option 1" << endl;
									for (unsigned long k4 = k3 + 1; k4 < totalNumberOfClusters; k4++) {
										if (k2 != k4) {
											// cluster(k4)
											ClusteringPtr cTemp =
													NeighborhoodSearch::process2optCombination(*g,
															clustering, problem, k1,
															k2, k3, k4, n,
															i, j);
											// cTemp->printClustering();
											Imbalance newImbalance = cTemp->getImbalance();
											if (newImbalance < bestImbalance) {
												// cout << "Better solution found in 2-neighborhood." << endl;
												cBest.reset();
												cBest = cTemp;
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
									if((problem.getType() == ClusteringProblem::CC_PROBLEM) ||
															( (problem.getType() == ClusteringProblem::RCC_PROBLEM)
																	&& (clustering->getNumberOfClusters() < k) )) {
										// cout << "Option 2" << endl;
										ClusteringPtr cTemp =
												NeighborhoodSearch::process2optCombination(*g,
														clustering, problem, k1, k2, k3,
														Clustering::NEW_CLUSTER,
														n, i, j);
										// cTemp->printClustering();
										Imbalance newImbalance = cTemp->getImbalance();
										if (newImbalance < bestImbalance) {
											// cout << "Better solution found in 2-neighborhood." << endl;
											cBest.reset();
											cBest = cTemp;
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
							if((problem.getType() == ClusteringProblem::CC_PROBLEM) ||
									( (problem.getType() == ClusteringProblem::RCC_PROBLEM)
											&& (clustering->getNumberOfClusters() < k) )) {
								// cout << "Option 3" << endl;
								for (unsigned long k4 = uninc(), contk4 = 0; contk4 < numberOfClustersInInterval; contk4++) {
									if (k2 != k4) {
										// cluster(k4)
										ClusteringPtr cTemp =
												NeighborhoodSearch::process2optCombination(*g,
														clustering, problem, k1, k2,
														Clustering::NEW_CLUSTER,
														k4, n, i, j);
										// cTemp->printClustering();
										Imbalance newImbalance = cTemp->getImbalance();
										if (newImbalance < bestImbalance) {
											// cout << "Better solution found in 2-neighborhood." << endl;
											cBest.reset();
											cBest = cTemp;
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
							if((problem.getType() == ClusteringProblem::CC_PROBLEM) ||
									( (problem.getType() == ClusteringProblem::RCC_PROBLEM)
											&& (clustering->getNumberOfClusters() + 2 <= k) )) {
								// cout << "Option 4" << endl;
								ClusteringPtr cTemp = NeighborhoodSearch::process2optCombination(*g,
										clustering, problem, k1, k2,
										Clustering::NEW_CLUSTER,
										Clustering::NEW_CLUSTER, n, i, j);
								// cTemp->printClustering();
								Imbalance newImbalance = cTemp->getImbalance();
								if (newImbalance < bestImbalance) {
									// cout << "Better solution found in 2-neighborhood." << endl;
									cBest.reset();
									cBest = cTemp;
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

ClusteringPtr NeighborhoodSearch::process2optCombination(SignedGraph& g, Clustering* clustering,
		ClusteringProblem& problem, int k1, int k2, int k3, int k4, int n, int i, int j) {

         // cout << "2-opt-comb: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << ", " << i << ", " << j << endl;
         // clustering->printClustering();
         ClusteringPtr cTemp = make_shared < Clustering > (*clustering);
         int nc = cTemp->getNumberOfClusters();
	 // increments number of tested combinations
	 numberOfTestedCombinations++;
         // the offset caused by cluster deletions
         // removes node i from cluster1 and inserts in cluster3
         // TODO check if the removal of node i destroys cluster1
         // cout << "k3" << endl;
         cTemp->removeNodeFromCluster(g, problem, i, k1);
         // recalculates the number of clusters, as one of them may have been removed
         int newnc1 = cTemp->getNumberOfClusters();
         if(newnc1 < nc) {
                 // cluster k1 has been removed
                 if(k2 >= k1) { k2--; assert(k2 >= 0); }
                 if(k3 >= k1) { k3--; assert(k3 >= 0); }
                 if(k4 >= k1) { k4--; /* assert(k4 >= 0); */ }
         }
         if (k3 > k1) {
                 // inserts i in existing cluster k3
                 cTemp->addNodeToCluster(g, problem, i, k3);
         } else {
                 // inserts i in a new cluster (alone)
                 cTemp->addCluster(g, problem, i);
         }
         // cout << "k4" << endl;
         // removes node j from cluster2 and inserts in cluster4
         cTemp->removeNodeFromCluster(g, problem, j, k2);
         int newnc2 = cTemp->getNumberOfClusters();
         if(newnc2 < newnc1) {
                 // cout << "cluster k2 has been removed" << endl;
                 if(k4 >= k2) { k4--; assert(k4 >= 0); }
         }
         // cout << "Node removed" << endl;
         if (k4 > k2) {
                 // inserts j in existing cluster k4
                 cTemp->addNodeToCluster(g, problem, j, k4);
         } else {
                 // inserts j in a new cluster (alone)
                 cTemp->addCluster(g, problem, j);
         }
         // cout << "Return" << endl;
         // Full recalculation of objective value if RCC Problem
         if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			Imbalance imb = problem.objectiveFunction(g, *cTemp);
			if(imb.getValue() != cTemp->getImbalance().getValue()) {
				BOOST_LOG_TRIVIAL(error) << "RCC obj function and delta do not match!";
			}
			cTemp->setImbalance(imb);
		}
         return cTemp;
 }

}
