/*
 * SequentialNeighborhoodSearch.cpp
 *
 *  Created on: May 25, 2014
 *      Author: mlevorato
 */

#include "include/SequentialNeighborhoodSearch.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include "include/CUDASearch.h"
#include "util/include/RandomUtil.h"

#include <boost/log/trivial.hpp>
#include <boost/timer/timer.hpp>

#include <limits>
#include <thrust/host_vector.h>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

namespace clusteringgraph {

using namespace thrust;
using namespace util;


SequentialNeighborhoodSearch::SequentialNeighborhoodSearch() {
	// TODO Auto-generated constructor stub

}

SequentialNeighborhoodSearch::~SequentialNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

Clustering SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, bool firstImprovementOnOneNeig) {
	unsigned long nc = clustering->getNumberOfClusters();
	numberOfTestedCombinations = 0;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
	}
	return SequentialNeighborhoodSearch::searchNeighborhood(l, g, clustering, problem,
			timeSpentSoFar, timeLimit, randomSeed, myRank, 0, nc, firstImprovementOnOneNeig, k);
}

Clustering SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialClusterIndex,
		unsigned long finalClusterIndex, bool firstImprovementOnOneNeig, unsigned long k) {

	if (l == 1) {  // 1-opt
		// Sequential search does first improvement in 1-opt, depending on parameter value
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, firstImprovementOnOneNeig, k);
	} else {  // 2-opt is always first improvement
		return this->search2opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialClusterIndex, finalClusterIndex, true, k);
	}
}

Clustering SequentialNeighborhoodSearch::search1opt(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long nc = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	// pre-calculates, in an array, to which cluster each vertex belongs
	std::vector<unsigned long> myCluster(n);
	for(unsigned long k = 0; k < nc; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace) << "[New local search] Begin...";
	// Array that stores the sum of edge weights between vertex i and all clusters
	std::vector<double> h_VertexClusterPosSum(n * (nc + 1));
	std::vector<double> h_VertexClusterNegSum(n * (nc + 1));
	for(int i = 0; i < n * (nc + 1); i++) {
		h_VertexClusterPosSum[i] = double(0);
		h_VertexClusterNegSum[i] = double(0);
	}
	// For each vertex, creates a list of in and out edges
	matrix<double> matrix (n, nc + 1);
	for(int i = 0; i < n; i++) {
		for(int k = 0; k < nc + 1; k++) {
			matrix(i, k) = double(0);
		}
	}
	int i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			count++; edge++;
			if(weight > 0) {
				h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
				matrix(i, myCluster[j]) += fabs(weight);
			}
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			count++; edge++;
			if(weight > 0) {
					h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
					matrix(i, myCluster[j]) += fabs(weight);
			}
		}
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k2) of local search
	std::vector<uint> h_randomIndex(n);
	for(uint i = 0; i < n; i++) {
		h_randomIndex[i] = RandomUtil::next(initialSearchIndex, finalSearchIndex);
	}
	BOOST_LOG_TRIVIAL(trace) << "[New local search] Ready to calculate...";

	numberOfTestedCombinations = 0;
	int bestDestCluster = -1;
	int bestSrcVertex = -1;
	Imbalance bestImbalance = clustering->getImbalance();
	bool foundBetterSolution = false;
	for(uint i = 0; i < n; i++) {  // for each vertex i
		// vertex i is in cluster(k1)
		ulong k1 = myCluster[i];
		// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
		// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
		// REMOVAL of vertex i from cluster k1 -> avoids recalculating
		//   the same thing in the k2 (destination cluster) loop
		double _negativeSumK1 = -h_VertexClusterNegSum[i*(nc+1) + k1];
		double _positiveSumK1 = 0.0;
		for(uint k = 0; k < nc; k++) {
			if(k != k1) {
					_positiveSumK1 -= h_VertexClusterPosSum[i*(nc+1) + k];
			}
		}
		// Random initial vertex
		uint k2 = h_randomIndex[i];
		uint range = finalSearchIndex - initialSearchIndex + 1;
		for(uint countK2 = 0; countK2 < range; countK2++) {  // cluster(k2)
			if(k2 != k1) {
				// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
				double negativeSum = _negativeSumK1 + h_VertexClusterNegSum[i*(nc+1) + k2];
				double positiveSum = _positiveSumK1;
				for(uint k = 0; k < nc; k++) {
					if(k != k2) {
						positiveSum += h_VertexClusterPosSum[i*(nc+1) + k];
					}
				}
				numberOfTestedCombinations++;
				if(clustering->getImbalance().getValue() + positiveSum + negativeSum < bestImbalance.getValue()) {  // improvement in imbalance
					bestImbalance = Imbalance(positiveSum + clustering->getImbalance().getPositiveValue(), negativeSum + clustering->getImbalance().getNegativeValue());
					bestDestCluster = k2;
					bestSrcVertex = i;
					foundBetterSolution = true;
					if(firstImprovement) {
						break;
					}
				}
			}
			// loop increment rule
			k2++;
			if(k2 > finalSearchIndex) {
				k2 = initialSearchIndex;
			}
		}
		if(firstImprovement and foundBetterSolution) {
			break;
		}
	}
	// Reproduce the best clustering found using host data structures
	if(bestSrcVertex >= 0) {
		std::vector<unsigned long> myCluster2(myCluster);
		int previousCluster = myCluster2[bestSrcVertex];
		myCluster2[bestSrcVertex] = bestDestCluster;
		int h_nc = nc;
		// validate the current number of clusters in the solution
		if(bestDestCluster == h_nc) {  // Option 1 (easier): destcluster == nc, i.e., bestSrcVertex is going to a standalone cluster
			h_nc++;  // just increment the number of clusters
		}
		// Option 2 (harder): bestSrcVertex's previous cluster no longer exists (removed its only element)
		// Check if there are other vertices in the previous cluster of bestSrcVertex
		bool found = false;
		for(int i = 0; i < n; i++) {
			if(myCluster2[i] == previousCluster) {
				found = true;
				break;
			}
		}
		if(not found) {  // previousCluster is empty, has to be removed
			for(int i = 0; i < n; i++) {
				if(myCluster2[i] > previousCluster) {
					myCluster2[i]--;
				}
			}
			h_nc--;
		}

		/*
		cout << "[New local search] Processing complete. Best result: vertex " << bestSrcVertex << " (from cluster " << myCluster[bestSrcVertex]
					<< ") goes to cluster " << bestDestCluster << " with I(P) = " << bestImbalance.getValue() << " " << bestImbalance.getPositiveValue() << " " << bestImbalance.getNegativeValue() << endl;
					*/
		Clustering newClustering(*clustering);
		int k1 = myCluster[bestSrcVertex];
		int k2 = bestDestCluster;
		bool newClusterK2 = (k2 == nc);
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex, k1);
		if(not newClusterK2) {  // existing cluster k2
			if((newClustering.getNumberOfClusters() < nc) && (k2 >= k1)) {
				// cluster k1 has been removed
				// cout << "Cluster " << k1 << " removed." << endl;
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2 - 1);
			} else {
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2);
			}
		} else {  // new cluster k2
			newClustering.addCluster(*g, problem, bestSrcVertex);
		}
		cBest = newClustering;

		std::vector<unsigned long> myCluster3(n);
		for(unsigned long k = 0; k < newClustering.getNumberOfClusters(); k++) {  // for each cluster k
			BoolArray clusterK = newClustering.getCluster(k);
			for(unsigned long i = 0; i < n; i++) {  // for each vertex i
				if(clusterK[i]) {  // vertex i is in cluster k
					myCluster3[i] = k;
				}
			}
		}
		// assert that myCluster3 equals myCluster2
		bool equals = true;
		for(int i = 0; i < n; i++) {
			if(myCluster2[i] != myCluster3[i]) {
				equals = false;
				cerr << "Failed on vertex " << i << " " << myCluster2[i] << " " << myCluster3[i] << endl;
				break;
			}
		}
		assert(equals);
		/*
		if(newClustering.getImbalance().getValue() != bestImbalance.getValue()) {
			BOOST_LOG_TRIVIAL(error) << "New and old objective function values DO NOT MATCH! Correct = " << newClustering.getImbalance().getValue()
					<< " vs obtained = " << bestImbalance.getValue();
		}*/
		BOOST_LOG_TRIVIAL(debug) << "[New local search] Validation. Best result: I(P) = " << newClustering.getImbalance().getValue() << " "
				<< newClustering.getImbalance().getPositiveValue() << " " << newClustering.getImbalance().getNegativeValue();
	} else {
		BOOST_LOG_TRIVIAL(debug) << "[New local search] Validation. No improvement.";
	}
	// returns the best combination found in 1-opt
	return cBest;
}

Clustering SequentialNeighborhoodSearch::search2opt(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	unsigned long nc = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	// pre-calculates, in an array, to which cluster each vertex belongs
	std::vector<unsigned long> myCluster(n);
	for(unsigned long k = 0; k < nc; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace) << "[New local search 2-opt] Begin...";
	std::vector<double> h_VertexClusterPosSum(n * (nc + 1));
	std::vector<double> h_VertexClusterNegSum(n * (nc + 1));
	for(int i = 0; i < n * (nc + 1); i++) {
		h_VertexClusterPosSum[i] = double(0);
		h_VertexClusterNegSum[i] = double(0);
	}
	// For each vertex, creates a list of in and out edges
	matrix<double> matrix (n, nc + 1);
	for(int i = 0; i < n; i++) {
		for(int k = 0; k < nc + 1; k++) {
			matrix(i, k) = double(0);
		}
	}
	int i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			count++; edge++;
			if(weight > 0) {
				h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
				matrix(i, myCluster[j]) += fabs(weight);
			}
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			count++; edge++;
			if(weight > 0) {
					h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
					matrix(i, myCluster[j]) += fabs(weight);
			}
		}
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k3) of local search
	std::vector<uint> h_randomIndex(n);
	for(uint idx = 0; idx < n; idx++) {
		h_randomIndex[idx] = RandomUtil::next(0, nc);
	}

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
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			// vertex i is in cluster(k1)
			ulong k1 = myCluster[i];
			// vertex j is in cluster(k2)
			ulong k2 = myCluster[j];
			// Option 1: vertex i is moved from k1 to another existing cluster k3 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k3 = nc)
			// REMOVAL of vertex i from cluster k1 -> avoids recalculating
			//   the same thing in the k2 (destination cluster) loop
			double negativeSumK1 = -h_VertexClusterNegSum[i*(nc+1) + k1];
			double positiveSumK1 = 0.0;
			for(uint k = 0; k < nc; k++) {
					if(k != k1) {
							positiveSumK1 -= h_VertexClusterPosSum[i*(nc+1) + k];
					}
			}
			// Random initial vertex
			uint k3 = h_randomIndex[i];
			for(uint countK3 = 0; countK3 <= nc; countK3++) {  // cluster(k3)
				if(k3 != k1) {
					// calculates the cost of removing vertex i from cluster k1 and inserting into cluster k3
					double negativeSum = negativeSumK1 + h_VertexClusterNegSum[i*(nc+1) + k3];
					double positiveSum = positiveSumK1;
					for(uint k = 0; k < nc; k++) {
						if(k != k3) {
							positiveSum += h_VertexClusterPosSum[i*(nc+1) + k];
						}
					}
					for(uint k4 = k3 + 1; k4 <= nc; k4++) {  // cluster(k4)
						if(k4 != k2) {
							// calculates the cost of removing vertex j from cluster k2 and inserting into cluster k4
							double negativeSum2 = negativeSum + h_VertexClusterNegSum[j*(nc+1) + k4];
							double positiveSum2 = positiveSum;
							for(uint k = 0; k < nc; k++) {
								if(k != k4) {
									positiveSum2 += h_VertexClusterPosSum[j*(nc+1) + k];
								}
							}
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
						}
					}
					if(firstImprovement and foundBetterSolution) {
						break;
					}
				}
				// loop increment rule
				k3++;
				if(k3 > nc) {
					k3 = 0;
				}
			}
			if(firstImprovement and foundBetterSolution) {
				break;
			}
		}
		if(firstImprovement and foundBetterSolution) {
			break;
		}
	}
	// Reproduce the best clustering found using host data structures
	if(bestSrcVertex1 >= 0) {
		BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Processing complete. Best result: vertex " << bestSrcVertex1 << " (from cluster " << myCluster[bestSrcVertex1]
					<< ") goes to cluster " << bestDestCluster1 << " with I(P) = " << bestImbalance.getValue() << " " << bestImbalance.getPositiveValue() << " " << bestImbalance.getNegativeValue();
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
		BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Validation. Best result: I(P) = " << newClustering.getImbalance().getValue() << " "
				<< newClustering.getImbalance().getPositiveValue() << " " << newClustering.getImbalance().getNegativeValue();
	} else {
		BOOST_LOG_TRIVIAL(debug) << "[New local search 2-opt] Validation. No improvement.";
	}
	// returns the best combination found in 1-opt
	return cBest;
}

} /* namespace clusteringgraph */
