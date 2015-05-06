#include <limits>
#include <boost/make_shared.hpp>
#include <boost/timer/timer.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/optional.hpp>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <vector>

#include "include/NeighborhoodSearch.h"
#include "../util/include/MPIMessage.h"
#include "../util/include/RandomUtil.h"
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

using namespace boost::mpi;
using namespace std;
using namespace util;

#define EPS 10e-6

// probes for mpi VND interrupt message (first improvement from other node) and returns if it exists
// or probes for vnd output message, if this search is being done by a master process
#define MPI_IPROBE_RETURN(ret_val) \
optional< mpi::status > stat = world.iprobe(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG); \
if (stat) { \
        InputMessageParallelVND imsg; \
        mpi::status stat = world.recv(mpi::any_source, MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG, imsg); \
	return (ret_val); } \
optional< mpi::status > stat2 = world.iprobe(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_VND_TAG); \
if (stat2) { \
	return (ret_val); \
} 


namespace clusteringgraph {

Clustering NeighborhoodSearch::search1opt(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	mpi::communicator world;

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	unsigned long nc = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	ClusterArray myCluster = clustering->getClusterArray();

	if(problem.getType() == ClusteringProblem::CC_PROBLEM) {  // runs optimized local search for CC Problem
		// rebuilds the auxiliary matrices in case it is needed (between metaheuristic iterations)
		if(h_vertexClusterPosSum.empty() or h_vertexClusterNegSum.empty() or h_isNeighborCluster.empty()) {
			// pre-calculates, in a list, for each vertex, which clusters are neighbors of it (i.e. has edges)
			h_vertexClusterPosSum.resize(n * (nc + 1), double(0));
			h_vertexClusterNegSum.resize(n * (nc + 1), double(0));
			h_isNeighborCluster.resize(n * (nc + 1), 0);
			// BOOST_LOG_TRIVIAL(trace) << "Generating local search auxiliary matrices.";
			updateVertexClusterSumArrays(g, h_vertexClusterPosSum, h_vertexClusterNegSum, h_isNeighborCluster, clustering);
		}
		// BOOST_LOG_TRIVIAL(trace) << "Invoking optimized 1-opt CC local search...";
		return search1optCCProblem(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
		                myRank, initialSearchIndex, finalSearchIndex, firstImprovement, k,
		        		myCluster, h_isNeighborCluster, h_vertexClusterPosSum, h_vertexClusterNegSum);
	}
	numberOfTestedCombinations = 0;
	// for each vertex i, tries to move i to another cluster in myNeighborClusterList[i]
	// For each node i in cluster(k1)
	for (unsigned long i = randomUtil.next(initialSearchIndex, finalSearchIndex), cont2 = 0; cont2 < numberOfVerticesInInterval; cont2++) {
		// vertex i is in cluster(k1)
		unsigned long k1 = myCluster[i];
		// RCC Problem: must have exactly k clusters -> Cannot remove node from standalone cluster
		unsigned long s = clustering->getClusterSize(k1);
		if((problem.getType() == ClusteringProblem::RCC_PROBLEM) && (s == 1)) {
			continue;
		}
		// Option 1: node i is moved from k1 to another existing cluster k2 != k1
		for (unsigned long k2 = 0; k2 < nc; k2++) {
			if( (k1 != k2) && (h_isNeighborCluster[i+k2*n] > 0) ) {
				MPI_IPROBE_RETURN(cBest)
				// removes node i from cluster1 and inserts in cluster2
				Clustering cTemp = *clustering;

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
				// return if time limit is exceeded
				boost::timer::cpu_times end_time = timer.elapsed();
				double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
				// std::cout << timeSpentSoFar + localTimeSpent << endl;
				if(timeSpentSoFar + localTimeSpent >= timeLimit)  return (cBest);
			}
		}
		// Option 2: node i is moved to a new cluster, alone
		// removes node i from cluster k1 and inserts in newCluster
		// cout << "New clustering combination generated." << endl;
		if((problem.getType() == ClusteringProblem::CC_PROBLEM) && (s > 1)) {
			// this code is not executed for RCC Problem, as it must have exactly k clusters
			// this code is not executed if node i is to be moved from a standalone cluster to another standalone cluster (s == 1)
			Clustering cTemp = *clustering;
			//BOOST_LOG_TRIVIAL(trace) << "Option 2: Taking node " << i << " from " << k1 << " to new cluster.";
			cTemp.removeNodeFromCluster(*g, problem, i, k1);
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
		// increment rule
		i++;
		if(i > finalSearchIndex) {
			i = initialSearchIndex;
		}
	}
	// returns the best combination found in 1-opt
	return cBest;
}

Clustering NeighborhoodSearch::search2opt(SignedGraph* g,
        Clustering* clustering, ClusteringProblem* problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, unsigned long initialSearchIndex,
		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k) {

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	//BOOST_LOG_TRIVIAL(trace) << "Generating 2-opt neighborhood...";
	mpi::communicator world;
	unsigned long n = g->getN();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	unsigned long nc = clustering->getNumberOfClusters();
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance)
	Clustering cBest = *clustering;
	Clustering cTemp = cBest;
	Imbalance bestImbalance = cBest.getImbalance();

	ClusterArray myCluster = clustering->getClusterArray();
	std::vector<double> h_VertexClusterPosSum(n * (nc + 1), double(0));
	std::vector<double> h_VertexClusterNegSum(n * (nc + 1), double(0));
	// pre-calculates, in a list, for each vertex, which clusters are neighbors of it (i.e. has edges)
	std::vector< unordered_set<unsigned long> > myNeighborClusterList(n, unordered_set<unsigned long>());
	for(int i = 0; i < n; i++) {  // for each vertex i
		DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {
			int j = target(*f, g->graph);
			double weight = ((Edge*)f->get_property())->weight;
			if(myCluster[i] != myCluster[j]) {  // different cluster
				myNeighborClusterList[i].insert(myCluster[j]);
			}
			if(weight > 0) {
				h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			}
		}
		// iterates over in-edges of vertex i
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			if(myCluster[i] != myCluster[j]) {  // different cluster
				myNeighborClusterList[i].insert(myCluster[j]);
			}
			if(weight > 0) {
					h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			}
		}
	}
	if(problem->getType() == ClusteringProblem::CC_PROBLEM) {  // runs optimized local search for CC Problem
		return search2optCCProblem(g, clustering, *problem, timeSpentSoFar, timeLimit, randomSeed,
						myRank, initialSearchIndex, finalSearchIndex, firstImprovement, k,
						myCluster, myNeighborClusterList, h_VertexClusterPosSum, h_VertexClusterNegSum);
	}
	// for each vertex i, tries to move i to another cluster in myNeighborClusterList[i]
	// For each node i in cluster(k1)
	for (unsigned long i = randomUtil.next(initialSearchIndex, finalSearchIndex), conti = 0; conti < numberOfVerticesInInterval; conti++) {
		// cluster(k1)
		unsigned long k1 = myCluster[i];
		unsigned long s1 = clustering->getClusterSize(k1);
		if((problem->getType() == ClusteringProblem::RCC_PROBLEM) && (s1 == 1)) {
			// no clusters can be removed if RCC Problem
			continue;
		}
		// For each node j in cluster(k2)
		for (unsigned long j = randomUtil.next(0, n - 1), contj = 0; contj < n; contj++, j = (j + 1) % n) {
			if(j == i)  continue;
			// cluster(k2)
			unsigned long k2 = myCluster[j];
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
			// Option 1: node i is moved to another existing cluster k3 != k1
			for (unordered_set<unsigned long>::iterator itr = myNeighborClusterList[i].begin(); itr != myNeighborClusterList[i].end(); ++itr) {
				MPI_IPROBE_RETURN(cBest)
				// cluster(k3)
				unsigned long k3 = *itr;
				// Option 1: node i is moved to another existing cluster k3 and
				//           node j is moved to another existing cluster k4
				for (unordered_set<unsigned long>::iterator itr2 = myNeighborClusterList[j].begin(); itr2 != myNeighborClusterList[j].end(); ++itr2) {
					// cluster(k4)
					unsigned long k4 = *itr2;
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
			timer.stop();
			boost::timer::cpu_times end_time = timer.elapsed();
			double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
			timer.resume();
			// Option 3: node i is moved to a new cluster, alone, and
			//           node j is moved to another existing cluster k4
			if((problem->getType() == ClusteringProblem::CC_PROBLEM) && (s1 > 1)) {
				// this code is not executed in RCC Problem, as it must have exactly k clusters
				// this code is not executed if node i is to be moved from a standalone cluster to another standalone cluster
				for (unordered_set<unsigned long>::iterator itr2 = myNeighborClusterList[j].begin(); itr2 != myNeighborClusterList[j].end(); ++itr2) {
					// cluster(k4)
					unsigned long k4 = *itr2;
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
					// return if time limit is exceeded
					timer.stop();
					boost::timer::cpu_times end_time = timer.elapsed();
					double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
					if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
					timer.resume();
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
		}
		// return if time limit is exceeded
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
		if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
		timer.resume();
		// increment rule
		i++;
		if(i > finalSearchIndex) {
			i = initialSearchIndex;
		}
	}
	return cBest;
}

Clustering NeighborhoodSearch::search1optCCProblem(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k,
        		ClusterArray& myCluster,
				std::vector<long>& isNeighborCluster,
        		std::vector<double>& vertexClusterPosSum,
        		std::vector<double>& vertexClusterNegSum) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	mpi::communicator world;

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long nc = clustering->getNumberOfClusters();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	// random number generators used in loop randomization
	RandomUtil randomUtil;

	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	numberOfTestedCombinations = 0;
	int bestDestCluster = -1;
	int bestSrcVertex = -1;
	Imbalance bestImbalance2 = clustering->getImbalance();
	bool foundBetterSolution = false;
	// for each vertex i, tries to move i to another cluster in myNeighborClusterList[i]
	// For each node i in cluster(k1)
	for (unsigned long i = randomUtil.next(initialSearchIndex, finalSearchIndex), cont2 = 0; cont2 < numberOfVerticesInInterval; cont2++) {
		// vertex i is in cluster(k1)
		ulong k1 = myCluster[i];
		// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
		// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
		// REMOVAL of vertex i from cluster k1 -> avoids recalculating
		//   the same thing in the k2 (destination cluster) loop
		double _negativeSumK1 = -vertexClusterNegSum[i + k1 * n];
		bool timeLimitExceeded = false;
		// Node i is moved from k1 to another existing cluster k2 != k1, include a new cluster (k2 = nc)
		for (unsigned long k2 = 0; k2 <= nc; k2++) {
			if( (k1 != k2) && (h_isNeighborCluster[i+k2*n] > 0) ) {
				MPI_IPROBE_RETURN(cBest)
				// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
				double negativeSum = _negativeSumK1 + vertexClusterNegSum[i + k2 * n];
				double positiveSum = -(- vertexClusterPosSum[i + k1 * n]) + (-vertexClusterPosSum[i + k2 * n]);
				numberOfTestedCombinations++;
				if(clustering->getImbalance().getValue() + positiveSum + negativeSum < bestImbalance2.getValue()) {  // improvement in imbalance
					bestImbalance2 = Imbalance(positiveSum + clustering->getImbalance().getPositiveValue(),
							negativeSum + clustering->getImbalance().getNegativeValue());
					bestDestCluster = k2;
					bestSrcVertex = i;
					foundBetterSolution = true;
					if(firstImprovement) {
						break;
					}
				}
				// return if time limit is exceeded
				timer.stop();
				boost::timer::cpu_times end_time = timer.elapsed();
				double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
				if(timeSpentSoFar + localTimeSpent >= timeLimit){
					timeLimitExceeded = true;
					break;
				}
				timer.resume();
			}
		}
		if((firstImprovement and foundBetterSolution) or timeLimitExceeded) {
			break;
		}
		// increment rule
		i++;
		if(i > finalSearchIndex) {
			i = initialSearchIndex;
		}
	}
	// Reproduce the best clustering found using host data structures
	int new_nc = nc;
	Clustering newClustering(*clustering);
	if(bestSrcVertex >= 0) {
		// BOOST_LOG_TRIVIAL(debug) << "[New local search] Processing complete. Best result: vertex " << bestSrcVertex << " (from cluster " << myCluster[bestSrcVertex]
		//			<< ") goes to cluster " << bestDestCluster << " with I(P) = " << bestImbalance2.getValue() << " " << bestImbalance2.getPositiveValue() << " " << bestImbalance2.getNegativeValue();
		int k1 = myCluster[bestSrcVertex];
		int k2 = bestDestCluster;
		bool newClusterK2 = (k2 == nc);
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex, k1);
		if(not newClusterK2) {  // existing cluster k2
			if((newClustering.getNumberOfClusters() < nc) && (k2 >= k1)) {
				// cluster k1 has been removed
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2 - 1);
			} else {
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2);
			}
		} else {  // new cluster k2
			newClustering.addCluster(*g, problem, bestSrcVertex);
		}
		cBest = newClustering;

		// BOOST_LOG_TRIVIAL(trace) << "Updating auxiliary matrices with delta...";
		long nc = bestClustering.getNumberOfClusters();
		long new_nc = cBest.getNumberOfClusters();
		// CASO ESPECIAL 1: novo cluster
		if(new_nc > nc) {  // acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
			h_vertexClusterPosSum.resize(n * (new_nc + 1), double(0));
			h_vertexClusterNegSum.resize(n * (new_nc + 1), double(0));
			h_isNeighborCluster.resize(n * (new_nc + 1), 0);
		}
		this->updateVertexClusterSumArraysDelta(g, h_vertexClusterPosSum, h_vertexClusterNegSum, h_isNeighborCluster,
				cBest, nc, new_nc, bestSrcVertex, k1, k2);
		bestClustering = cBest;
		// CASO ESPECIAL 2: cluster removido
		if(new_nc < nc) {
			h_vertexClusterPosSum.resize(n * (new_nc + 1), double(0));
			h_vertexClusterNegSum.resize(n * (new_nc + 1), double(0));
			h_isNeighborCluster.resize(n * (new_nc + 1), 0);
		}
		// BOOST_LOG_TRIVIAL(trace) << "Updated.";
	} else {
		// BOOST_LOG_TRIVIAL(debug) << "[New local search] Validation. No improvement.";
	}

	// returns the best combination found in 1-opt
	return cBest;
}

Clustering NeighborhoodSearch::search2optCCProblem(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k,
        		ClusterArray& myCluster,
				std::vector< unordered_set<unsigned long> >& myNeighborClusterList,
				std::vector<double>& vertexClusterPosSum,
				std::vector<double>& vertexClusterNegSum) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	RandomUtil randomUtil;
	mpi::communicator world;
	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	unsigned long nc = clustering->getNumberOfClusters();
	for(int i = 0; i < n; i++) {  // inserts nc (aka new cluster) into the list of cluster neighbors of each vertex
		myNeighborClusterList[i].insert(nc);
	}
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	numberOfTestedCombinations = 0;
	// best destination clusters
	int bestDestCluster1 = -1;
	int bestDestCluster2 = -1;
	double originalPosImbalance = clustering->getImbalance().getPositiveValue();
	double originalNegImbalance = clustering->getImbalance().getNegativeValue();
	int bestSrcVertex1 = -1;
	int bestSrcVertex2 = -1;
	Imbalance bestImbalance = clustering->getImbalance();
	bool foundBetterSolution = false;
	// for each vertex i, tries to move i to another cluster in myNeighborClusterList[i]
	// For each node i in cluster(k1)
	for (unsigned long i = randomUtil.next(initialSearchIndex, finalSearchIndex), conti = 0; conti < numberOfVerticesInInterval; conti++) {
		// For each node j in cluster(k2)
		for (unsigned long j = randomUtil.next(0, n - 1), contj = 0; contj < n; contj++, j = (j + 1) % n) {
			if(j == i)  continue;
			// vertex i is in cluster(k1)
			ulong k1 = myCluster[i];
			// vertex j is in cluster(k2)
			ulong k2 = myCluster[j];
			unsigned long s2 = clustering->getClusterSize(k2);
			if(k1 == k2 && (s2 < 2)) {
				// source cluster k1 must have at least 2 elements
				continue;
			}
			// Option 1: vertex i is moved from k1 to another existing cluster k3 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k3 = nc)
			// REMOVAL of vertex i from cluster k1 -> avoids recalculating
			//   the same thing in the k2 (destination cluster) loop
			double negativeSumK1 = -vertexClusterNegSum[i*(nc+1) + k1];
			// Option 1: node i is moved to another existing cluster k3 != k1
			for (unordered_set<unsigned long>::iterator itr = myNeighborClusterList[i].begin(); itr != myNeighborClusterList[i].end(); ++itr) {
				MPI_IPROBE_RETURN(cBest)
				// cluster(k3)
				unsigned long k3 = *itr;
				if(k3 != k1) {
					// calculates the cost of removing vertex i from cluster k1 and inserting into cluster k3
					double negativeSum = negativeSumK1 + vertexClusterNegSum[i*(nc+1) + k3];
					double positiveSum = -(- vertexClusterPosSum[i*(nc+1) + k1]) + (-vertexClusterPosSum[i*(nc+1) + k3]);
					for (unordered_set<unsigned long>::iterator itr2 = myNeighborClusterList[j].begin(); itr2 != myNeighborClusterList[j].end(); ++itr2) {
						// cluster(k4)
						unsigned long k4 = *itr2;
						if(k4 != k2) {
							// calculates the cost of removing vertex j from cluster k2 and inserting into cluster k4
							double negativeSum2 = negativeSum - vertexClusterNegSum[j*(nc+1) + k2] + vertexClusterNegSum[j*(nc+1) + k4];
							double positiveSum2 = positiveSum + -(- vertexClusterPosSum[j*(nc+1) + k2]) + (-vertexClusterPosSum[j*(nc+1) + k4]);
							numberOfTestedCombinations++;
							if(originalPosImbalance + originalNegImbalance + positiveSum2 + negativeSum2 < bestImbalance.getValue()) {  // improvement in imbalance
								bestImbalance = Imbalance(positiveSum2 + originalPosImbalance, negativeSum2 + originalNegImbalance);
								bestDestCluster1 = k3;
								bestDestCluster2 = k4;
								bestSrcVertex1 = i;
								bestSrcVertex2 = j;
								foundBetterSolution = true;
								if(firstImprovement) {
									break;
								}
							}
							// return if time limit is exceeded
							timer.stop();
							boost::timer::cpu_times end_time = timer.elapsed();
							double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
							if(timeSpentSoFar + localTimeSpent >= timeLimit){
									return cBest;
							}
							timer.resume();
						}
					}
					if(firstImprovement and foundBetterSolution) {
						break;
					}
				}
			}
			if(firstImprovement and foundBetterSolution) {
				break;
			}
		}
		if(firstImprovement and foundBetterSolution) {
			break;
		}
		// increment rule
		i++;
		if(i > finalSearchIndex) {
			i = initialSearchIndex;
		}
	}
	// Reproduce the best clustering found using host data structures
	if(bestSrcVertex1 >= 0) {
		// BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Processing complete. Best result: vertex " << bestSrcVertex1 << " (from cluster " << myCluster[bestSrcVertex1]
		// 			<< ") goes to cluster " << bestDestCluster1 << " with I(P) = " << bestImbalance.getValue() << " " << bestImbalance.getPositiveValue() << " " << bestImbalance.getNegativeValue();
		Clustering newClustering(*clustering);
		int k1 = myCluster[bestSrcVertex1];
		int k3 = bestDestCluster1;
		int k2 = myCluster[bestSrcVertex2];
		int k4 = bestDestCluster2;
		bool newClusterK3 = (k3 == nc);
		bool newClusterK4 = (k4 == nc);
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex1, k1);
		int newnc1 = newClustering.getNumberOfClusters();
		if(not newClusterK3) {  // existing cluster k3
			if((newnc1 < nc) && (k3 >= k1)) {
				// cluster k1 has been removed
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex1, k3 - 1);
			} else {
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex1, k3);
			}
		} else {  // new cluster k3
			newClustering.addCluster(*g, problem, bestSrcVertex1);
		}
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex2, k2);
		int newnc2 = newClustering.getNumberOfClusters();
		if(not newClusterK4) {  // existing cluster k4
			if((newnc2 < newnc1) && (k4 >= k1)) {
				// cluster k2 has been removed
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex2, k4 - 1);
			} else {
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex2, k4);
			}
		} else {  // new cluster k4
			newClustering.addCluster(*g, problem, bestSrcVertex2);
		}
		cBest = newClustering;
		/*
		if(newClustering.getImbalance().getValue() != bestImbalance.getValue()) {
			BOOST_LOG_TRIVIAL(error) << "New-calculation and old-calculation objective function values DO NOT MATCH!";
		}*/
		// BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Validation. Best result: I(P) = " << newClustering.getImbalance().getValue() << " "
		// 		<< newClustering.getImbalance().getPositiveValue() << " " << newClustering.getImbalance().getNegativeValue();
	} else {
		// BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Validation. No improvement.";
	}
	// returns the best combination found in 2-opt
	return cBest;
}

void NeighborhoodSearch::updateVertexClusterSumArrays(SignedGraph* g, std::vector<double>& vertexClusterPosSumArray,
		std::vector<double>& vertexClusterNegSumArray, std::vector<long>& isNeighborClusterArray,	Clustering* clustering) {

	long n = g->getN();
	ClusterArray clusterArray = clustering->getClusterArray();
	long nc = clustering->getNumberOfClusters();

    for(int i = 0; i < n; i++) {
        // For each vertex i, stores the sum of edge weights between vertex i and all clusters
        for(int k = 0; k <= nc; k++) {
        	vertexClusterPosSumArray[k * n + i] = double(0);
        	vertexClusterNegSumArray[k * n + i] = double(0);
        	isNeighborClusterArray[k * n + i] = 0;
        }
        // in/out-edges of vertex i
        int ki = clusterArray[i];
        DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {
			int j = target(*f, g->graph);
			double weight = ((Edge*)f->get_property())->weight;
			int kj = clusterArray[j];

			if(kj >= 0) {
				if(ki != kj) {  // different cluster
					isNeighborClusterArray[kj * n + i]++;  // vertex i now has external connection to cluster kj
					isNeighborClusterArray[ki * n + j]++;  // vertex j now has external connection to cluster ki
				}
				if(weight > 0) {
					vertexClusterPosSumArray[kj * n + i] += fabs(weight);
				} else {
					vertexClusterNegSumArray[kj * n + i] += fabs(weight);
				}
			}
		}
		// iterates over in-edges of vertex i
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			int kj = clusterArray[j];

			if(kj >= 0) {
				if(ki != kj) {  // different cluster
					isNeighborClusterArray[kj * n + i]++;  // vertex i now has external connection to cluster kj
					isNeighborClusterArray[ki * n + j]++;  // vertex j now has external connection to cluster ki
				}
				if(weight > 0) {
					vertexClusterPosSumArray[kj * n + i] += fabs(weight);
				} else {
					vertexClusterNegSumArray[kj * n + i] += fabs(weight);
				}
			}
		}
		// preenche a possibilidade de se mover o vertice i para um cluster novo (k = nc)
		isNeighborClusterArray[nc * n + i] = 1;
    }
}

void NeighborhoodSearch::updateVertexClusterSumArraysDelta(SignedGraph* g, std::vector<double>& vertexClusterPosSumArray,
		std::vector<double>& vertexClusterNegSumArray, std::vector<long>& isNeighborClusterArray,
		Clustering& clustering, uint nc, uint new_nc, int i, int k1, int k2) {  // vertex i is being moved from cluster k1 to k2

	int n = g->getN();
	ClusterArray clusterArray = clustering.getClusterArray();
	if(new_nc > nc) {  // vertex i is being moved to a new cluster
		// move a fileira correspondente ao cluster k = nc na matriz de soma, shiftando os dados para a direita (nc + 1)
		for(int v = 0; v < n; v++) {
			vertexClusterPosSumArray[(new_nc) * n + v] = vertexClusterPosSumArray[(nc) * n + v];
			vertexClusterNegSumArray[(new_nc) * n + v] = vertexClusterNegSumArray[(nc) * n + v];
		}
		// zera a fileira movida anteriormente
		for(int v = 0; v < n; v++) {
			vertexClusterPosSumArray[(nc) * n + v] = 0.0;
			vertexClusterNegSumArray[(nc) * n + v] = 0.0;
		}
		// preenche a possibilidade de se mover todos os vertices para um cluster novo (k = new_nc)
		for(int v = 0; v < n; v++) {
			isNeighborClusterArray[v + new_nc * n] = 1;
			isNeighborClusterArray[v + nc * n] = 0;
		}
	}
	// in/out-edges of vertex i
	// isNeighborClusterArray[i+k1*n] = 1;  // vertex i now has external connection to cluster k1
	/* for(ulong k = 0; k < nc; k++) {
	 	isNeighborClusterArray[i + k * n] = 0;
	} */
	DirectedGraph::out_edge_iterator f, l;
	// For each out edge of i
	for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {
		int j = target(*f, g->graph);
		double weight = ((Edge*)f->get_property())->weight;
		int kj = clusterArray[j];

		isNeighborClusterArray[j + k1 * n] -= 2;
		if(kj != k2) {
			isNeighborClusterArray[i + kj * n]++;  // vertex i now has external connection to cluster kj
			isNeighborClusterArray[j + k2 * n]++;  // vertex j now has external connection to cluster k2
		} else {  // cluster[i] == cluster[j] == k2
			isNeighborClusterArray[i + kj * n] = 0;  // vertex i now has NO external connections to cluster kj
			isNeighborClusterArray[j + k2 * n] = 0;  // vertex j now has NO external connections to cluster k2
		}
		if(weight > 0) {
			vertexClusterPosSumArray[j + k1 * n] -= fabs(weight);
			vertexClusterPosSumArray[j + k2 * n] += fabs(weight);
		} else {
			vertexClusterNegSumArray[j + k1 * n] -= fabs(weight);
			vertexClusterNegSumArray[j + k2 * n] += fabs(weight);
		}
	}
	// iterates over in-edges of vertex i
	DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
	for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
		double weight = ((Edge*)f2->get_property())->weight;
		int j = source(*f2, g->graph);
		int kj = clusterArray[j];

		isNeighborClusterArray[j + k1 * n] -= 2;
		if(kj != k2) {
			isNeighborClusterArray[i + kj * n]++;  // vertex i now has external connection to cluster kj
			isNeighborClusterArray[j + k2 * n]++;  // vertex j now has external connection to cluster k2
		} else {  // cluster[i] == cluster[j] == k2
			isNeighborClusterArray[i + kj * n] = 0;  // vertex i now has NO external connections to cluster kj
			isNeighborClusterArray[j + k2 * n] = 0;  // vertex j now has NO external connections to cluster k2
		}
		if(weight > 0) {
			vertexClusterPosSumArray[j + k1 * n] -= fabs(weight);
			vertexClusterPosSumArray[j + k2 * n] += fabs(weight);
		} else {
			vertexClusterNegSumArray[j + k1 * n] -= fabs(weight);
			vertexClusterNegSumArray[j + k2 * n] += fabs(weight);
		}
	}

	if(new_nc < nc) {
		// remove a fileira correspondente ao cluster k1 na matriz de soma, shiftando os dados para a esquerda
		for(int k = k1 + 1; k <= nc; k++) {
			for(int v = 0; v < n; v++) {
				vertexClusterPosSumArray[(k - 1) * n + v] = vertexClusterPosSumArray[(k) * n + v];
				vertexClusterNegSumArray[(k - 1) * n + v] = vertexClusterNegSumArray[(k) * n + v];
				isNeighborClusterArray[(k - 1) * n + v] = isNeighborClusterArray[(k) * n + v];
			}
		}
	}
}

}

