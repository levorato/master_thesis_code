#pragma once

#include <thrust/host_vector.h>

namespace clusteringgraph {

using namespace thrust;

	typedef unsigned char ubyte;
	typedef unsigned short ushort;
	typedef unsigned int uint;
	typedef unsigned long ulong;


	extern "C" bool runSimpleSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m,
			thrust::host_vector<unsigned long>& h_destcluster, thrust::host_vector<float>& h_destPosFunctionValue,
			thrust::host_vector<float>& h_destNegFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
			bool firstImprovement, thrust::host_vector<unsigned long>& h_destNumComb);

}
