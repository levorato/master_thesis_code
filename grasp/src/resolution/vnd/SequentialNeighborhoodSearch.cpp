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
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	// pre-calculates, in an array, to which cluster each vertex belongs
	std::vector<unsigned long> myCluster(n);
	for(unsigned long k = 0; k < totalNumberOfClusters; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace) << "[New local search] Begin...";
	// graph structure (adapted adjacency list)
	std::vector<double> h_weights(2 * m);  // in/out edge weights
	std::vector<int> h_dest(2 * m);  // edge destination (vertex j)
	std::vector<int> h_numedges(n);  // number of edges of each vertex i
	std::vector<int> h_offset(n);  // initial edge number for vertex i
	// Array that stores the sum of edge weights between vertex i and all clusters
	std::vector<double> h_VertexClusterPosSum(n * totalNumberOfClusters);
	std::vector<double> h_VertexClusterNegSum(n * totalNumberOfClusters);
	for(int i = 0; i < n * totalNumberOfClusters; i++) {
		h_VertexClusterPosSum[i] = 0.0;
		h_VertexClusterNegSum[i] = 0.0;
	}
	// For each vertex, creates a list of in and out edges
	std::vector< std::vector<string> > arestasI(n);
	int i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		h_offset[i] = offset;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
			if(weight > 0) {
				h_VertexClusterPosSum[i * totalNumberOfClusters + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * totalNumberOfClusters + myCluster[j]] += fabs(weight);
				stringstream ss;
				ss << j << ", " << weight;
				arestasI[j].push_back(ss.str());
			}
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
			if(weight > 0) {
					h_VertexClusterPosSum[i * totalNumberOfClusters + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * totalNumberOfClusters + myCluster[j]] += fabs(weight);
					stringstream ss;
					ss << j << ", " << weight;
					arestasI[j].push_back(ss.str());
			}
		}
		h_numedges[i] = count;
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k2) of local search
	ulong nc = clustering->getNumberOfClusters();
	std::vector<uint> h_randomIndex(n);
	for(uint i = 0; i < n; i++) {
		h_randomIndex[i] = RandomUtil::next(initialSearchIndex, finalSearchIndex);
	}
	BOOST_LOG_TRIVIAL(trace) << "[New local search] Ready to calculate...";

	numberOfTestedCombinations = 0;
	int bestDestCluster = -1;
	int bestSrcVertex = -1;
	Imbalance bestImbalance = clustering->getImbalance();
	for(uint i = 0; i < n; i++) {  // for each vertex i
		// vertex i is in cluster(k1)
		ulong k1 = myCluster[i];
		uint numedges = h_numedges[i];
		uint offset = h_offset[i];
		// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
		// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
		// REMOVAL of vertex i from cluster k1 -> avoids recalculating
		//   the same thing in the k2 (destination cluster) loop
		double negativeSumK1 = 0.0, positiveSumK1 = 0.0;
		/*
		double negativeSumK1 = -h_VertexClusterNegSum[i*nc + k1];
		double positiveSumK1 = 0.0;
		for(uint k = 0; k < nc; k++) {
				if(k != k1) {
						positiveSumK1 -= h_VertexClusterPosSum[i*nc + k];
				}
		}*/
		// in/out-edges of vertex i

		ulong count = offset + numedges;
		for (ulong edgenum = offset; edgenum < count; edgenum++) {
			int targ = h_dest[edgenum];
			double weight = h_weights[edgenum];
			// REMOVAL from cluster k1: subtracts imbalance
			if(myCluster[targ] == k1) {  // same cluster
				if(weight < 0) {
					negativeSumK1 -= fabs(weight);
				}
			} else {  // diff cluster
				if(weight > 0) {
					positiveSumK1 -= weight;
				}
			}
		}
		// Random initial vertex
		uint k2 = h_randomIndex[i];
		uint range = finalSearchIndex - initialSearchIndex + 1;
		for(uint countK2 = 0; countK2 < range; countK2++) {  // cluster(k2)
			if(k2 != k1) {
				// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
				double negativeSum = negativeSumK1, positiveSum = positiveSumK1;

				double negativeSum2 = negativeSumK1 + h_VertexClusterNegSum[i*nc + k2];
				double positiveSum2 = positiveSumK1;
				for(uint k = 0; k < nc; k++) {
					if(k != k2) {
						positiveSum2 += h_VertexClusterPosSum[i*nc + k];
					}
				}
				std::vector<string> arestas;
				ulong count = offset + numedges;
				// in/out-edges of vertex i
				for (ulong edgenum = offset; edgenum < count; edgenum++) {
					int targ = h_dest[edgenum];
					double weight = h_weights[edgenum];
					// ADDITION to cluster k2 != k1: adds imbalance
					if(myCluster[targ] == k2) {  // same cluster
						if(weight < 0) {
							negativeSum += fabs(weight);
							stringstream ss;
							ss << targ + ", ";
							arestas.push_back(ss.str());
						}
					} else {  // diff cluster
						if(weight > 0) {
							positiveSum += weight;
						}
					}
				}
				if(negativeSum != negativeSum2) {
					cout << "calculo original\n";
					for(int z = 0; z < arestas.size(); z++) {
						cout << arestas[z] << "\n";
					}
					cout << "calculo novo\n";
					for(int z = 0; z < arestasI[i].size(); z++) {
						cout << arestasI[i][z] << "\n";
					}
				}
				numberOfTestedCombinations++;
				if(clustering->getImbalance().getValue() + positiveSum + negativeSum < bestImbalance.getValue()) {  // improvement in imbalance
					bestImbalance = Imbalance(positiveSum + clustering->getImbalance().getPositiveValue(), negativeSum + clustering->getImbalance().getNegativeValue());
					bestDestCluster = k2;
					bestSrcVertex = i;
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
	}
	// Reproduce the best clustering found using host data structures
	if(bestSrcVertex >= 0) {
		BOOST_LOG_TRIVIAL(debug) << "[New local search] Processing complete. Best result: vertex " << bestSrcVertex << " (from cluster " << myCluster[bestSrcVertex]
					<< ") goes to cluster " << bestDestCluster << " with I(P) = " << bestImbalance.getValue() << " " << bestImbalance.getPositiveValue() << " " << bestImbalance.getNegativeValue();
		Clustering newClustering(*clustering);
		int k1 = myCluster[bestSrcVertex];
		int k2 = bestDestCluster;
		bool newClusterK2 = (k2 == totalNumberOfClusters);
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex, k1);
		if(not newClusterK2) {  // existing cluster k2
			if((newClustering.getNumberOfClusters() < totalNumberOfClusters) && (k2 >= k1)) {
				// cluster k1 has been removed
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2 - 1);
			} else {
				newClustering.addNodeToCluster(*g, problem, bestSrcVertex, k2);
			}
		} else {  // new cluster k2
			newClustering.addCluster(*g, problem, bestSrcVertex);
		}
		cBest = newClustering;
		if(newClustering.getImbalance().getValue() != bestImbalance.getValue()) {
			BOOST_LOG_TRIVIAL(error) << "New and old objective function values DO NOT MATCH! Correct = " << newClustering.getImbalance().getValue()
					<< " vs obtained = " << bestImbalance.getValue();
		}
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
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	// pre-calculates, in an array, to which cluster each vertex belongs
	boost::shared_ptr<unsigned long[]> myCluster(new unsigned long[n]);
	for(unsigned long k = 0; k < totalNumberOfClusters; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}

	BOOST_LOG_TRIVIAL(trace) << "[CUDA 2-opt] Begin transfer to device...";
	// A -> Transfer to device
	// transfers the myClusters array to CUDA device
	unsigned long numberOfChunks = n * n;  // the search space for each pair of vertices (dest cluster)
	unsigned long* cluster = myCluster.get();
	thrust::host_vector<unsigned long> h_mycluster(cluster, cluster + n);
	// objective function value
	thrust::host_vector<float> h_functionValue(2);
	h_functionValue[0] = clustering->getImbalance().getPositiveValue();
	h_functionValue[1] = clustering->getImbalance().getNegativeValue();
	// destination (result) host vectors
	thrust::host_vector<unsigned long> h_destcluster1(numberOfChunks);  // destination cluster (k3)
	thrust::host_vector<unsigned long> h_destcluster2(numberOfChunks);  // destination cluster (k4)
	thrust::host_vector<float> h_destPosFunctionValue(numberOfChunks);  // positive imbalance value
	thrust::host_vector<float> h_destNegFunctionValue(numberOfChunks);  // negative imbalance value
	thrust::host_vector<unsigned long> h_destNumComb(numberOfChunks);  // number of combinations
	// graph structure (adapted adjacency list)
	thrust::host_vector<float> h_weights(2 * m);  // in/out edge weights
	thrust::host_vector<int> h_dest(2 * m);  // edge destination (vertex j)
	thrust::host_vector<int> h_numedges(n);  // number of edges of each vertex i
	thrust::host_vector<int> h_offset(n);  // initial edge number for vertex i
	// For each vertex, creates a list of in and out edges
	int i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		h_offset[i] = offset;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		h_numedges[i] = count;
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k3) of local search
	thrust::host_vector<uint> h_randomIndex(numberOfChunks);
	ulong nparts = n;
	ulong chunkSize = n;
	for(uint idx = 0; idx < numberOfChunks; idx++) {
		h_randomIndex[idx] = RandomUtil::next(0, totalNumberOfClusters);
	}
	// TODO transform into class constant
	// number of threads per block
	unsigned short threadsCount = 512;

	// Pass raw array and its size to kernel
	run2optSearchKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, n, m,
			h_destcluster1, h_destcluster2, h_destPosFunctionValue, h_destNegFunctionValue, threadsCount, totalNumberOfClusters,
			numberOfChunks, firstImprovement, h_destNumComb, h_randomIndex);

	// Returns the best value found, iterating through all threads results
	int bestSrcVertex1 = -1;
	int bestSrcVertex2 = -1;
	int bestDestCluster1 = -1;
	int bestDestCluster2 = -1;
	double bestImbalance = clustering->getImbalance().getValue();
	// To simulate sequential first improvement, chooses a random initial index for the following loop
	for(int i = RandomUtil::next(0, numberOfChunks - 1), cont = 0; cont < numberOfChunks; cont++, i = (i + 1) % numberOfChunks) {
		if(h_destPosFunctionValue[i] + h_destNegFunctionValue[i] < bestImbalance) {
			bestSrcVertex1 = i % n;
			bestSrcVertex2 = i / n;
			bestDestCluster1 = h_destcluster1[i];
			bestDestCluster2 = h_destcluster2[i];
			bestImbalance = h_destPosFunctionValue[i] + h_destNegFunctionValue[i];
			if(firstImprovement)  break;
		}
	}
	// Sums the number of combinations visited
	numberOfTestedCombinations = 0;
	for(int i = 0; i < numberOfChunks; i++) {
		numberOfTestedCombinations += h_destNumComb[i];
	}

	// Reproduce the best clustering found using host data structures
	if(bestSrcVertex1 >= 0) {
		BOOST_LOG_TRIVIAL(debug) << "[CUDA] Processing complete. Best result: vertex " << bestSrcVertex1 << " (from cluster " << myCluster[bestSrcVertex1]
					<< ") goes to cluster " << bestDestCluster1 << " with I(P) = " << bestImbalance << " " << h_destPosFunctionValue[bestSrcVertex1] << " " << h_destNegFunctionValue[bestSrcVertex1];
		Clustering newClustering(*clustering);
		int k1 = myCluster[bestSrcVertex1];
		int k3 = bestDestCluster1;
		int k2 = myCluster[bestSrcVertex2];
		int k4 = bestDestCluster2;
		bool newClusterK3 = (k3 == totalNumberOfClusters);
		bool newClusterK4 = (k4 == totalNumberOfClusters);
		newClustering.removeNodeFromCluster(*g, problem, bestSrcVertex1, k1);
		int newnc1 = newClustering.getNumberOfClusters();
		if(not newClusterK3) {  // existing cluster k3
			if((newnc1 < totalNumberOfClusters) && (k3 >= k1)) {
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
		if(newClustering.getImbalance().getValue() != bestImbalance) {
			BOOST_LOG_TRIVIAL(error) << "CUDA and CPU objective function values DO NOT MATCH!";
		}
		BOOST_LOG_TRIVIAL(debug) << "[CUDA] Validation. Best result: I(P) = " << newClustering.getImbalance().getValue() << " "
				<< newClustering.getImbalance().getPositiveValue() << " " << newClustering.getImbalance().getNegativeValue();
	} else {
		BOOST_LOG_TRIVIAL(debug) << "[CUDA] Validation. No improvement.";
	}
	// returns the best combination found in 1-opt
	return cBest;
}

} /* namespace clusteringgraph */
