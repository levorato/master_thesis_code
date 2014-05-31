/*
 * CUDANeighborhoodSearch.cpp
 *
 *  Created on: May 25, 2014
 *      Author: mlevorato
 */

#include "include/CUDANeighborhoodSearch.h"
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


CUDANeighborhoodSearch::CUDANeighborhoodSearch() {
	// TODO Auto-generated constructor stub

}

CUDANeighborhoodSearch::~CUDANeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

Clustering CUDANeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
        Clustering* clustering, ClusteringProblem& problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, bool firstImprovementOnOneNeig) {

	// Resets the number of combinations tested on neighborhood search
	numberOfTestedCombinations = 0;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
	}

	unsigned long threadsCount = 120;

	// Splits the processing in (n / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = 0;
	unsigned long remainingVertices = 0;
	if(g->getN() > threadsCount) {
		sizeOfChunk = (unsigned long)((double)g->getN() / threadsCount);
		remainingVertices = g->getN() % threadsCount;
	} else {
		sizeOfChunk = 1;
		threadsCount = g->getN();
		remainingVertices = 0;
	}
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	// invoke the 1-opt CUDA method
	BOOST_LOG_TRIVIAL(debug) << "Waiting for CUDA...\n";
	Clustering bestClustering;
	double bestValue = numeric_limits<double>::infinity();
	bestClustering = searchNeighborhood(l, g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed, myRank, 0, threadsCount*sizeOfChunk - 1,
			firstImprovementOnOneNeig, k);

	// the leader (me) does its part of the work too
	BOOST_LOG_TRIVIAL(trace) << "Total number of vertices is " << g->getN() << endl;
	BOOST_LOG_TRIVIAL(trace) << "CUDA Parallelization status: " << threadsCount <<
			" search slaves (threads) will process " << sizeOfChunk << " vertice(s) each one." << endl;

	BOOST_LOG_TRIVIAL(debug) << "[Parallel Search CUDA] Best solution found: Obj = " << bestClustering.getImbalance().getValue() << endl;
	bestClustering.printClustering();
	return bestClustering;
}

Clustering CUDANeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialSearchIndex,
		unsigned long finalSearchIndex, bool firstImprovementOnOneNeig, unsigned long k) {
	// Unused method
	assert(initialSearchIndex < g->getN());
	assert(finalSearchIndex < g->getN());

	if (l == 1) {  // 1-opt
		// Parallel search always does best improvement in 1-opt
		// Therefore, parameter firstImprovementOnOneNeig is ignored
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, false, k);
	} else {  // 2-opt
		// TODO implement global first improvement in parallel 2-opt
		// first improvement found in one process must break all other processes' loop
		// IMPORTANT: Parallel VND does first improvement on 2-opt if CC Problem
		// Or best improvement on 2-opt if RCC Problem
		bool firstImprovementOn2Opt = true;
		if(problem.getType() == ClusteringProblem::CC_PROBLEM) {
			firstImprovementOn2Opt = true;
		} else if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			// TODO check if this is still necessary with GRASP RCC
			// firstImprovementOn2Opt = false;
		}
		return this->search2opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, firstImprovementOn2Opt, k);
	}
}

Clustering CUDANeighborhoodSearch::search1opt(SignedGraph* g,
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
	boost::shared_ptr<unsigned long[]> myCluster(new unsigned long[n]);
	for(unsigned long k = 0; k < nc; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}

	BOOST_LOG_TRIVIAL(trace) << "[CUDA] Begin transfer to device...";
	// A -> Transfer to device
	// transfers the myClusters array to CUDA device
	unsigned long numberOfChunks = n * 4;  // the search space for each vertex (dest cluster) will be split into 4 chunks
	unsigned long* cluster = myCluster.get();
	thrust::host_vector<unsigned long> h_mycluster(cluster, cluster + n);
	// objective function value
	thrust::host_vector<float> h_functionValue(2);
	h_functionValue[0] = clustering->getImbalance().getPositiveValue();
	h_functionValue[1] = clustering->getImbalance().getNegativeValue();
	// destination (result) host vectors
	thrust::host_vector<unsigned long> h_destcluster(numberOfChunks);  // destination cluster (k2)
	thrust::host_vector<float> h_destPosFunctionValue(numberOfChunks);  // positive imbalance value
	thrust::host_vector<float> h_destNegFunctionValue(numberOfChunks);  // negative imbalance value
	thrust::host_vector<unsigned long> h_destNumComb(numberOfChunks);  // number of combinations
	// graph structure (adapted adjacency list)
	thrust::host_vector<float> h_weights(2 * m);  // in/out edge weights
	thrust::host_vector<int> h_dest(2 * m);  // edge destination (vertex j)
	thrust::host_vector<int> h_numedges(n);  // number of edges of each vertex i
	thrust::host_vector<int> h_offset(n);  // initial edge number for vertex i
	// Array that stores the sum of edge weights between vertex i and all clusters
	thrust::host_vector<float> h_VertexClusterPosSum(n * (nc+1));
	thrust::host_vector<float> h_VertexClusterNegSum(n * (nc+1));
	for(int i = 0; i < n * (nc+1); i++) {
		h_VertexClusterPosSum[i] = 0.0;
		h_VertexClusterNegSum[i] = 0.0;
	}
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
			if(weight > 0) {
				h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
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
					h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			}
		}
		h_numedges[i] = count;
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k2) of local search
	thrust::host_vector<uint> h_randomIndex(numberOfChunks);
	ulong nparts = numberOfChunks / n;
	ulong chunkSize = ulong(ceil((float)(nc + 1.0) / nparts));
	for(uint idx = 0; idx < numberOfChunks; idx++) {
		uint part = idx / n;
		uint initialK2 = part * chunkSize;
		uint finalK2 = (part + 1) * chunkSize - 1;
		if(initialK2 < nc + 1) {
			if(finalK2 >= nc + 1) {
				finalK2 = nc;
			}
			h_randomIndex[idx] = RandomUtil::next(initialK2, finalK2);
		} else {
			h_randomIndex[idx] = 0;
		}
	}
	// TODO transform into class constant
	// number of threads per block
	unsigned short threadsCount = 256;

	// Pass raw array and its size to kernel
	run1optSearchKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, n, m,
			h_destcluster, h_destPosFunctionValue, h_destNegFunctionValue, threadsCount, nc,
			numberOfChunks, firstImprovement, h_destNumComb, h_randomIndex, h_VertexClusterPosSum, h_VertexClusterNegSum);

	// Finds the best value found, iterating through all threads results
	int bestSrcVertex = -1;
	int bestDestCluster = -1;
	double bestImbalance = clustering->getImbalance().getValue();
	// To simulate sequential first improvement, chooses a random initial index for the following loop
	for(int i = RandomUtil::next(0, numberOfChunks - 1), cont = 0; cont < numberOfChunks; cont++, i = (i + 1) % numberOfChunks) {
		if(h_destPosFunctionValue[i] + h_destNegFunctionValue[i] < bestImbalance) {
			bestSrcVertex = i % n;
			bestDestCluster = h_destcluster[i];
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
	if(bestSrcVertex >= 0) {
		BOOST_LOG_TRIVIAL(debug) << "[CUDA] Processing complete. Best result: vertex " << bestSrcVertex << " (from cluster " << myCluster[bestSrcVertex]
					<< ") goes to cluster " << bestDestCluster << " with I(P) = " << bestImbalance << " " << h_destPosFunctionValue[bestSrcVertex] << " " << h_destNegFunctionValue[bestSrcVertex];
		Clustering newClustering(*clustering);
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

Clustering CUDANeighborhoodSearch::search2opt(SignedGraph* g,
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
	boost::shared_ptr<unsigned long[]> myCluster(new unsigned long[n]);
	for(unsigned long k = 0; k < nc; k++) {  // for each cluster k
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
	// Array that stores the sum of edge weights between vertex i and all clusters
	thrust::host_vector<float> h_VertexClusterPosSum(n * (nc+1));
	thrust::host_vector<float> h_VertexClusterNegSum(n * (nc+1));
	for(int i = 0; i < n * (nc+1); i++) {
		h_VertexClusterPosSum[i] = 0.0;
		h_VertexClusterNegSum[i] = 0.0;
	}
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
			if(weight > 0) {
				h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
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
					h_VertexClusterPosSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum[i * (nc+1) + myCluster[j]] += fabs(weight);
			}
		}
		h_numedges[i] = count;
		offset += count;
	}
	// random number vector used for initial cluster index (destination -> k3) of local search
	thrust::host_vector<uint> h_randomIndex(numberOfChunks);
	ulong nparts = n;
	ulong chunkSize = n;
	for(uint idx = 0; idx < numberOfChunks; idx++) {
		h_randomIndex[idx] = RandomUtil::next(0, nc);
	}
	// TODO transform into class constant
	// number of threads per block
	unsigned short threadsCount = 512;

	// Pass raw array and its size to kernel
	run2optSearchKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, n, m,
			h_destcluster1, h_destcluster2, h_destPosFunctionValue, h_destNegFunctionValue, threadsCount, nc,
			numberOfChunks, firstImprovement, h_destNumComb, h_randomIndex, h_VertexClusterPosSum, h_VertexClusterNegSum);

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
