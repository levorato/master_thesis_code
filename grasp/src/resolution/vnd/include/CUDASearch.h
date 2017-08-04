#pragma once

#include <thrust/host_vector.h>
#include <vector>
#include "problem/include/ClusteringProblem.h"
#include "../../construction/include/ConstructClustering.h"
#include "graph/include/Graph.h"

namespace clusteringgraph {

using namespace thrust;
using namespace problem;
using namespace resolution::construction;
using namespace clusteringgraph;




	bool run1optSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				unsigned long n, unsigned long m, unsigned short threadsCount, unsigned long& nc, 
				unsigned long numberOfChunks, bool firstImprovement,
				thrust::host_vector<unsigned int>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum, unsigned int& bestSrcVertex, unsigned int& destcluster,
				float& destFunctionValue, const long& timeSpentSoFar, const unsigned int& l);

	bool runVNDKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				unsigned long n, unsigned long m, unsigned short threadsCount, unsigned long& nc, 
				unsigned long numberOfChunks, bool firstImprovement,
				thrust::host_vector<unsigned int>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum, thrust::host_vector<unsigned int> &h_neighbor_cluster,
				std::vector<unsigned int>& sourceVertexList, std::vector<unsigned int>& destinationClusterList,
				float& destPositiveImbalance, float& destNegativeImbalance, const long& timeSpentSoFar, const unsigned int& l);

	extern "C" bool run2optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				unsigned long n, unsigned long m, thrust::host_vector<unsigned long>& h_destcluster1,
				thrust::host_vector<unsigned long>& h_destcluster2,
				thrust::host_vector<float>& h_destFunctionValue, unsigned short threadsCount, 
				unsigned long nc, unsigned long numberOfChunks,
				bool firstImprovement,
				thrust::host_vector<unsigned int>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum);

	bool runGRASPKernel(ClusteringProblem& problem, ConstructClustering &construct,
				SignedGraph *g, int processRank, unsigned long timeLimit, int iter,
				thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				unsigned long n, unsigned long m, unsigned short threadsCount, bool firstImprovement,
				Clustering& result, int &totalIterations, double& timeSpentConstruct, double& timeSpentGRASP,
				stringstream &constructivePhaseResults, stringstream &iterationResults);

	bool runILSKernel(ClusteringProblem& problem, ConstructClustering &construct,
				SignedGraph *g, int processRank, unsigned long timeLimit,
				const int& iterMax, const int& iterMaxILS, const int& perturbationLevelMax,
				thrust::host_vector<float>& h_weights, thrust::host_vector<long>& h_dest,
				thrust::host_vector<long>& h_numedges, thrust::host_vector<long>& h_offset,
				unsigned long n, unsigned long m, unsigned short threadsCount, bool firstImprovement,
				Clustering& result, long &totalIterations, double& timeSpentConstruct, double& timeSpentILS,
				stringstream &constructivePhaseResults, stringstream &iterationResults);

	bool runConstructKernel(unsigned long randomSeed, thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				unsigned long n, unsigned long m, unsigned long nc, unsigned short threadsCount, 
				thrust::host_vector<unsigned long>& h_newcluster,
				double& imbalance);

}
