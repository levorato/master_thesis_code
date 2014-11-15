/*
 * CUDAVariableNeighborhoodDescent.cpp
 *
 *  Created on: 5/6/2014
 *      Author: Mario Levorato
 */

#include "include/CUDAVariableNeighborhoodDescent.h"
#include "problem/include/ClusteringProblem.h"
#include "graph/include/Imbalance.h"
#include "util/include/RandomUtil.h"
// #include "../include/LocalSearch.h"
#include "graph/include/Graph.h"
#include "graph/include/SequentialNeighborhoodSearch.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/round.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/log/trivial.hpp>
#include "include/CUDASearch.h"
#include "util/include/RandomUtil.h"
#include <limits>
#include <cstdio>
#include <thrust/host_vector.h>
#include <vector>

using namespace problem;
using namespace boost;
using namespace clusteringgraph;
using namespace thrust;
using namespace util;

namespace resolution {
namespace vnd {

CUDAVariableNeighborhoodDescent::CUDAVariableNeighborhoodDescent(NeighborhoodSearch &neighborhoodSearch,
		unsigned long seed, const int &lsize, const bool& firstImprovement1Opt, const long &tlimit) :
		VariableNeighborhoodDescent(neighborhoodSearch, seed, lsize, firstImprovement1Opt, tlimit),
				timeResults(), timeSum(0.0) {
	// TODO Auto-generated constructor stub

}

CUDAVariableNeighborhoodDescent::~CUDAVariableNeighborhoodDescent() {
	// TODO Auto-generated destructor stub
}

Clustering CUDAVariableNeighborhoodDescent::localSearch(SignedGraph *g, Clustering& Cc, const int& graspIteration,
		ClusteringProblem& problem, const long& timeSpentSoFar, const int& myRank) {
	this->timeSpentOnLocalSearch = 0.0;
	this->numberOfTestedCombinations = 0;
	BOOST_LOG_TRIVIAL(debug)<< "CUDA local search VND...";

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long nc = Cc.getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best Cc combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = Cc;
	ClusterArray myCluster = Cc.getClusterArray();

	BOOST_LOG_TRIVIAL(trace) << "[CUDA] Begin transfer to device...";
	// A -> Transfer to device
	// transfers the myClusters array to CUDA device
	unsigned long numberOfChunks = n * (nc + 1);  // the search space for each vertex (dest cluster) will be split into 4 chunks
	thrust::host_vector<unsigned long> h_mycluster(myCluster);
	// objective function value
	thrust::host_vector<float> h_functionValue(1);
	h_functionValue[0] = Cc.getImbalance().getValue();
	// destination (result) host vectors
	thrust::host_vector<unsigned long> h_destcluster(numberOfChunks);  // destination cluster (k2)
	thrust::host_vector<float> h_destPosFunctionValue(numberOfChunks);  // positive imbalance value
	thrust::host_vector<float> h_destNegFunctionValue(numberOfChunks);  // negative imbalance value
	thrust::host_vector<unsigned long> h_destNumComb(numberOfChunks);  // number of combinations
	// Array that stores the sum of edge weights between vertex i and all clusters
	thrust::host_vector<float> h_VertexClusterPosSum(n * (nc+1));
	thrust::host_vector<float> h_VertexClusterNegSum(n * (nc+1));
	for(int i = 0; i < n * (nc+1); i++) {
		h_VertexClusterPosSum[i] = 0.0;
		h_VertexClusterNegSum[i] = 0.0;
	}
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
			if(weight > 0) {
				h_VertexClusterPosSum[myCluster[j] * n + i] += fabs(weight);
			} else {
				h_VertexClusterNegSum[myCluster[j] * n + i] += fabs(weight);
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
					h_VertexClusterPosSum[myCluster[j] * n + i] += fabs(weight);
			} else {
					h_VertexClusterNegSum[myCluster[j] * n + i] += fabs(weight);
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
	unsigned short threadsCount = 256;  // limited by shared memory size

	// Pass raw array and its size to kernel
	std::vector<uint> sourceVertexList;
	std::vector<uint> destinationClusterList;
	float bestImbalance = Cc.getImbalance().getValue();
	int l = 1;
	runVNDKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, n, m, threadsCount, nc,
			numberOfChunks, firstImprovementOnOneNeig, h_randomIndex, h_VertexClusterPosSum, h_VertexClusterNegSum,
			sourceVertexList, destinationClusterList, bestImbalance, timeSpentSoFar, l);

	// validate arrays calculation
/*
	thrust::host_vector<float> h_VertexClusterPosSum2(n * (nc+1));
	thrust::host_vector<float> h_VertexClusterNegSum2(n * (nc+1));
	for(int i = 0; i < n * (nc+1); i++) {
		h_VertexClusterPosSum2[i] = 0.0;
		h_VertexClusterNegSum2[i] = 0.0;
	}
	i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			count++; edge++;
			if(weight > 0) {
				h_VertexClusterPosSum2[i * (nc+1) + h_mycluster[j]] += fabs(weight);
			} else {
				h_VertexClusterNegSum2[i * (nc+1) + h_mycluster[j]] += fabs(weight);
			}
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			count++; edge++;
			if(weight > 0) {
					h_VertexClusterPosSum2[i * (nc+1) + h_mycluster[j]] += fabs(weight);
			} else {
					h_VertexClusterNegSum2[i * (nc+1) + h_mycluster[j]] += fabs(weight);
			}
		}
		offset += count;
	}
	bool equal = true;
	for(int i = 0; i < n; i++) {
		// printf("i = %d\n", i);
		for(int k = 0; k <= nc; k++) {
			// printf("k = %d\n", k);
			// printf("Pos: correct: %.2f CUDA: %.2f ; ", h_VertexClusterPosSum2[i * (nc + 1) + k], h_VertexClusterPosSum[i * (nc + 1) + k]);
			// printf("Neg: %.2f %.2f\n", h_VertexClusterNegSum2[i * (nc + 1) + k], h_VertexClusterNegSum[i * (nc + 1) + k]);
			if(h_VertexClusterPosSum2[i * (nc + 1) + k] != h_VertexClusterPosSum[i * (nc + 1) + k]) {
				equal = false;
				// printf("Failed on pos i = %d and k = %d\n", i, k);

				// break;
			}
			if(h_VertexClusterNegSum2[i * (nc + 1) + k] != h_VertexClusterNegSum[i * (nc + 1) + k]) {
				equal = false;
				// printf("Failed on neg i = %d and k = %d\n", i, k);

				// break;
			}
		}
	}
	assert(equal);
*/
	// To simulate sequential first improvement, chooses a random initial index for the following loop
	/*
	for(int i = RandomUtil::next(0, numberOfChunks - 1), cont = 0; cont < numberOfChunks; cont++, i = (i + 1) % numberOfChunks) {
		if(h_destPosFunctionValue[i] + h_destNegFunctionValue[i] < bestImbalance) {
			bestSrcVertex = i % n;
			bestDestCluster = h_destcluster[i];
			bestImbalance = h_destPosFunctionValue[i] + h_destNegFunctionValue[i];
			if(firstImprovement)  break;
		}
	}*/

	// Reproduce the best Cc found using host data structures
	if(bestImbalance < Cc.getImbalance().getValue()) {
		// BOOST_LOG_TRIVIAL(info) << "Number of clusters is " << nc;
		// cout << "Number of clusters is " << nc << endl;
		/*
		for(int i = 0; i < n; i++) {
			cout << h_mycluster[i] << " ";
		}
		cout << endl;
		*/
		Clustering newClustering(Cc);
		for(int x = 0; x < sourceVertexList.size(); x++) {
			unsigned long bestSrcVertex = sourceVertexList[x];
			unsigned long bestDestCluster = destinationClusterList[x];
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
			nc = newClustering.getNumberOfClusters();
		}
                if(newClustering.getImbalance().getValue() != bestImbalance) {
                        BOOST_LOG_TRIVIAL(error) << "CUDA and CPU objective function values DO NOT MATCH! " << bestImbalance;
                }
                BOOST_LOG_TRIVIAL(debug) << "[CUDA] Validation. Best result: I(P) = " << newClustering.getImbalance().getValue() << " "
                                << newClustering.getImbalance().getPositiveValue() << " " << newClustering.getImbalance().getNegativeValue();

		cBest = newClustering;
	} else {
		BOOST_LOG_TRIVIAL(debug) << "[CUDA] Validation. No improvement.";
	}
	// returns the result of VND
	return cBest;
}

void CUDAVariableNeighborhoodDescent::measureTimeResults(Clustering &CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	Imbalance imbalance = CStar.getImbalance();
	timeResults << fixed << setprecision(4) << (timeSpentInSearch + timeSpentOnLocalSearch) << "," << imbalance.getValue()
			<< "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << "," << CStar.getNumberOfClusters()
			<< "," << (iteration+1) << "\n";
}

void CUDAVariableNeighborhoodDescent::notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	measureTimeResults(CStar, timeSpentOnLocalSearch, iteration);
}

} /* namespace vnd */
} /* namespace resolution */
