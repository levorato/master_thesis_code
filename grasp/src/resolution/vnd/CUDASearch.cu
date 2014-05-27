#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <stdio.h>

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
	__global__ void simpleSearchKernel(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, const float* funcArray, uint n, uint m,
		ulong* destClusterArray, float* destPosImbArray, float* destNegImbArray, ulong nc) {

		// compute new objective function value
		uint i = blockIdx.x * blockDim.x + threadIdx.x;
		
		if(i < n) {  // kernel executes local search for vertex i
			// vertex i is in cluster(k1)
			ulong k1 = clusterArray[i];
			uint bestDestCluster = k1;
			// stores the best imbalance found so far
			float bestValuePos, bestValueNeg;
			bestValuePos = funcArray[0];
			bestValueNeg = funcArray[1];
			float originalImbalance = funcArray[0] + funcArray[1];
			// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
			for (ulong k2 = 0; k2 <= nc; k2++) {  // cluster(k2)
				if(k2 != k1) {
					// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
					float negativeSum = 0.0, positiveSum = 0.0;
					ulong count = offsetArray[i] + numArray[i];
					// out-edges of vertex i
					for (ulong edgenum = offsetArray[i]; edgenum < count; edgenum++) {   
						int targ = destArray[edgenum];
						float weight = weightArray[edgenum];
						// REMOVAL from cluster k1: subtracts imbalance
						if(clusterArray[targ] == k1) {  // same cluster
							if(weight < 0) {
								negativeSum -= fabs(weight);
							}
						} else {  // diff cluster
							if(weight > 0) {
								positiveSum -= weight;
							}
						}
					}
					// in-edges of vertex i
					for(uint j = 0; j < n; j++) {
						ulong count2 = offsetArray[j] + numArray[j];
						for (ulong edgenum = offsetArray[j]; edgenum < count2; edgenum++) {
							if(destArray[edgenum] == i) {  // edge (j, i)   
								float weight = weightArray[edgenum];
								// REMOVAL from cluster k1: subtracts imbalance
								if(clusterArray[j] == k1) {  // same cluster
									if(weight < 0) {
										negativeSum -= fabs(weight);
									}
								} else {  // diff cluster
									if(weight > 0) {
										positiveSum -= weight;
									}
								}
							}
						}
					}
					// temporatily changes vertex i's cluster from k1 to k2
					
					count = offsetArray[i] + numArray[i];
					// out-edges of vertex i
					for (ulong edgenum = offsetArray[i]; edgenum < count; edgenum++) {   
						int targ = destArray[edgenum];
						float weight = weightArray[edgenum];
						// ADDITION to cluster k2 != k1: adds imbalance
						if(clusterArray[targ] == k2) {  // same cluster
							if(weight < 0) {
								negativeSum += fabs(weight);
							}
						} else {  // diff cluster
							if(weight > 0) {
								positiveSum += weight;
							}
						}
					}
					// in-edges of vertex i
					for(uint j = 0; j < n; j++) {
						ulong count2 = offsetArray[j] + numArray[j];
						for (ulong edgenum = offsetArray[j]; edgenum < count2; edgenum++) {
							if(destArray[edgenum] == i) {  // edge (j, i)   
								float weight = weightArray[edgenum];
								// ADDITION to cluster k2 != k1: adds imbalance
								if(clusterArray[j] == k2) {  // same cluster
									if(weight < 0) {
										negativeSum += fabs(weight);
									}
								} else {  // diff cluster
									if(weight > 0) {
										positiveSum += weight;
									}
								}
							}
						}
					}
					if(originalImbalance + positiveSum + negativeSum < bestValuePos + bestValueNeg) {  // improvement in imbalance
						bestValuePos = positiveSum + funcArray[0];
						bestValueNeg = negativeSum + funcArray[1];
						bestDestCluster = k2;
					}
				}
			}
			// updates thread i / vertex i imbalance result
			destPosImbArray[i] = bestValuePos;
			destNegImbArray[i] = bestValueNeg;
			destClusterArray[i] = bestDestCluster;
		}
	}
	
	/// Runs a kernel for simple byte-per-cell world evaluation.
	bool runSimpleSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m,
			thrust::host_vector<unsigned long>& h_destcluster, thrust::host_vector<float>& h_destPosFunctionValue,
			thrust::host_vector<float>& h_destNegFunctionValue, ushort threadsCount, ulong nc) {

		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		// destination vectors
		thrust::device_vector<float> d_destPosFunctionValue = h_destPosFunctionValue;
		thrust::device_vector<float> d_destNegFunctionValue = h_destNegFunctionValue;
		thrust::device_vector<unsigned long> d_destCluster = h_destcluster;
	
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		unsigned long* destClusterArray = thrust::raw_pointer_cast( &d_destCluster[0] );
		float* destPosImbArray = thrust::raw_pointer_cast( &d_destPosFunctionValue[0] );
		float* destNegImbArray = thrust::raw_pointer_cast( &d_destNegFunctionValue[0] );

		// size_t reqBlocksCount = (n * (nc - 1)) / threadsCount;
		// ushort blocksCount = (ushort)std::min((size_t)32768, reqBlocksCount);
		int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
    	// printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsCount);
		// <<<blocksCount, threadsCount>>>
		simpleSearchKernel<<<blocksPerGrid, threadsCount>>>(weightArray, destArray, numArray, offsetArray, 
				clusterArray, funcArray, uint(n), uint(m), destClusterArray, destPosImbArray, destNegImbArray, ulong(nc));
		
		checkCudaErrors(cudaDeviceSynchronize());

		h_destcluster = d_destCluster;
		h_destPosFunctionValue = d_destPosFunctionValue;
		h_destNegFunctionValue = d_destNegFunctionValue;
		
		return true;
	}
		
}
