#pragma once

#include <thrust/host_vector.h>

namespace clusteringgraph {

using namespace thrust;

	typedef unsigned char ubyte;
	typedef unsigned short ushort;
	typedef unsigned int uint;
	typedef unsigned long ulong;


	extern "C" bool run1optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m, ushort threadsCount, ulong nc, ulong numberOfChunks, bool firstImprovement,
			thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
			thrust::host_vector<float>& h_VertexClusterNegSum, uint& bestSrcVertex, uint& destcluster, float& destFunctionValue);

	extern "C" bool run2optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
				ulong n, ulong m, thrust::host_vector<unsigned long>& h_destcluster1,
				thrust::host_vector<unsigned long>& h_destcluster2, thrust::host_vector<float>& h_destPosFunctionValue,
				thrust::host_vector<float>& h_destNegFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
				bool firstImprovement, thrust::host_vector<unsigned long>& h_destNumComb,
				thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
				thrust::host_vector<float>& h_VertexClusterNegSum);

}
