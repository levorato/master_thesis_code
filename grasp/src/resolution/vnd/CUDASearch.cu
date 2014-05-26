#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>

// CUDA facts:
//
// On devices of compute capability 2.x and beyond, 32-bit integer multiplication is natively supported,
// but 24-bit integer multiplication is not. __[u]mul24 is therefore implemented using multiple instructions
// and should not be used.
//
// Integer division and modulo operation are costly: below 20 instructions on devices of compute capability 2.x and
// higher. They can be replaced with bitwise operations in some cases: If n is a power of 2, (i/n) is equivalent to
// (i>>log2(n)) and (i%n) is equivalent to (i&(n-1)); the compiler will perform these conversions if n is literal.

namespace clusteringgraph {
	/// CUDA kernel for simple byte-per-cell world evaluation.
	///
	/// @param lifeData  Linearized 2D array of life data with byte-per-cell density.
	/// @param worldWidth  Width of life world in cells (bytes).
	/// @param worldHeight  Height of life world in cells (bytes).
	/// @param resultLifeData  Result buffer in the same format as input.
	__global__ void simpleSearchKernel(const double* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, const double* funcArray, uint n, uint m,
		ulong* destClusterArray, double* destFuncArray) {
		

		// compute new objective function value
		
/*
		uint worldSize = worldWidth * worldHeight;
		
		for (uint cellId = blockIdx.x * blockDim.x + threadIdx.x;
				cellId < worldSize;
				cellId += blockDim.x * gridDim.x) {

			uint x = cellId % worldWidth;
			uint yAbs = cellId - x;
			
			// Count alive cells.
			uint aliveCells = lifeData[xLeft + yAbsUp] + lifeData[x + yAbsUp] + lifeData[xRight + yAbsUp]
				+ lifeData[xLeft + yAbs] + lifeData[xRight + yAbs]
				+ lifeData[xLeft + yAbsDown] + lifeData[x + yAbsDown] + lifeData[xRight + yAbsDown];

			resultLifeData[x + yAbs] = aliveCells == 3 || (aliveCells == 2 && lifeData[x + yAbs]) ? 1 : 0;
		} */
	}
	
	/// Runs a kernel for simple byte-per-cell world evaluation.
	bool runSimpleSearchKernel(thrust::host_vector<double>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<double>& h_functionValue,
			ulong n, ulong m,
			thrust::host_vector<unsigned long>& h_destcluster, thrust::host_vector<double>& h_destFunctionValue,
			ushort threadsCount) {

		thrust::device_vector<double> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		thrust::device_vector<double> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		// destination vectors
		thrust::device_vector<double> d_destFunctionValue(2);
		thrust::device_vector<unsigned long> d_destCluster(n);
	
		double* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		double* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		unsigned long* destClusterArray = thrust::raw_pointer_cast( &d_destCluster[0] );
		double* destFuncArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );

		/*
		if ((worldWidth * worldHeight) % threadsCount != 0) {
			return false;
		}*/

		size_t reqBlocksCount = n / threadsCount;
		ushort blocksCount = (ushort)std::min((size_t)32768, reqBlocksCount);

		simpleSearchKernel<<<blocksCount, threadsCount>>>(weightArray, destArray, numArray, offsetArray, 
				clusterArray, funcArray, uint(n), uint(m), destClusterArray, destFuncArray);
		
		checkCudaErrors(cudaDeviceSynchronize());

		return true;
	}
		
}
