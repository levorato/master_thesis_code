#pragma once

#include <thrust/host_vector.h>
#include <vector>

namespace clusteringgraph {

using namespace thrust;

	typedef unsigned char ubyte;
	typedef unsigned short ushort;
	typedef unsigned int uint;
	typedef unsigned long ulong;


	bool run1optSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				ulong n, ulong m, ushort threadsCount, ulong& nc, ulong numberOfChunks, bool firstImprovement,
				thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum, uint& bestSrcVertex, uint& destcluster,
				float& destFunctionValue, const long& timeSpentSoFar, const unsigned int& l);

	bool runVNDKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
				thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
				thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				ulong n, ulong m, ushort threadsCount, ulong& nc, ulong numberOfChunks, bool firstImprovement,
				thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum, std::vector<uint>& sourceVertexList,
				std::vector<uint>& destinationClusterList,
				float& destFunctionValue, const long& timeSpentSoFar, const unsigned int& l);

	extern "C" bool run2optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				ulong n, ulong m, thrust::host_vector<unsigned long>& h_destcluster1,
				thrust::host_vector<unsigned long>& h_destcluster2,
				thrust::host_vector<float>& h_destFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
				bool firstImprovement,
				thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum);

}
