#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <stdio.h>
#include <curand_kernel.h>

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
		ulong* destClusterArray, float* destPosImbArray, float* destNegImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, ulong* destNumCombArray, uint* randomIndexArray) {

		uint idx = blockIdx.x * blockDim.x + threadIdx.x;
		uint i = idx % n;  // kernel executes local search for vertex i
		uint part = idx / n;  // kernel executes local search for vertex i, for the 'part' partition
		uint nparts = numberOfChunks / n;
		int offset = offsetArray[i];
		int numedges = numArray[i];
		// shared memory
		extern __shared__ long s[];       // n longs
		long *s_cluster = s;              // n longs
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		
		if(idx < numberOfChunks) {
			ulong numberOfTestedCombinations = 0;
			// vertex i is in cluster(k1)
			ulong k1 = s_cluster[i];
			uint bestDestCluster = k1;
			// stores the best imbalance found so far
			float bestValuePos, bestValueNeg;
			float originalPosImbalance = funcArray[0];
			float originalNegImbalance = funcArray[1];
			bestValuePos = originalPosImbalance;
			bestValueNeg = originalNegImbalance;
			// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
			ulong chunkSize = ulong(ceil((float)(nc + 1.0) / nparts));
			uint initialK2 = part * chunkSize;
			uint finalK2 = (part + 1) * chunkSize - 1;
			if(initialK2 < nc + 1 && i < n) {
				if(finalK2 >= nc + 1) {  // last chunk
					finalK2 = nc;
				}
				// REMOVAL of vertex i from cluster k1 -> avoids recalculating 
				//   the same thing in the k2 (destination cluster) loop
				float negativeSumK1 = 0.0, positiveSumK1 = 0.0;
				// in/out-edges of vertex i
				ulong count = offset + numedges;
				for (ulong edgenum = offset; edgenum < count; edgenum++) {   
					int targ = destArray[edgenum];
					float weight = weightArray[edgenum];
					// REMOVAL from cluster k1: subtracts imbalance
					if(s_cluster[targ] == k1) {  // same cluster
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
				uint k2 = randomIndexArray[idx];
				for(uint countK2 = 0; countK2 < chunkSize; countK2++) {  // cluster(k2)
					if(k2 != k1) {
						// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
						float negativeSum = negativeSumK1, positiveSum = positiveSumK1;
						ulong count = offset + numedges;
						// in/out-edges of vertex i
						for (ulong edgenum = offset; edgenum < count; edgenum++) {   
							int targ = destArray[edgenum];
							float weight = weightArray[edgenum];
							// ADDITION to cluster k2 != k1: adds imbalance
							if(s_cluster[targ] == k2) {  // same cluster
								if(weight < 0) {
									negativeSum += fabs(weight);
								}
							} else {  // diff cluster
								if(weight > 0) {
									positiveSum += weight;
								}
							}
						}
						numberOfTestedCombinations++;
						if(originalPosImbalance + originalNegImbalance + positiveSum + negativeSum < bestValuePos + bestValueNeg) {  // improvement in imbalance
							bestValuePos = positiveSum + originalPosImbalance;
							bestValueNeg = negativeSum + originalNegImbalance;
							bestDestCluster = k2;
							if(firstImprovement) {
								destPosImbArray[idx] = bestValuePos;
								destNegImbArray[idx] = bestValueNeg;
								destClusterArray[idx] = bestDestCluster;
								destNumCombArray[idx] = numberOfTestedCombinations;
								return;
							}
						}
					}
					// loop increment rule
					k2++;
					if(k2 > finalK2) {
						k2 = initialK2;
					}
				}
				// updates thread idx / vertex i imbalance result
				destPosImbArray[idx] = bestValuePos;
				destNegImbArray[idx] = bestValueNeg;
				destClusterArray[idx] = bestDestCluster;
				destNumCombArray[idx] = numberOfTestedCombinations;
			} else {
				destPosImbArray[idx] = originalPosImbalance;
				destNegImbArray[idx] = originalNegImbalance;
				destClusterArray[idx] = k1;
				destNumCombArray[idx] = 0;
			}
		}
	}
	
	/// Runs a kernel for simple byte-per-cell world evaluation.
	bool runSimpleSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m,
			thrust::host_vector<unsigned long>& h_destcluster, thrust::host_vector<float>& h_destPosFunctionValue,
			thrust::host_vector<float>& h_destNegFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
			bool firstImprovement, thrust::host_vector<unsigned long>& h_destNumComb,
			thrust::host_vector<uint>& h_randomIndex) {

		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<uint> d_randomIndex = h_randomIndex;
		// destination vectors
		thrust::device_vector<float> d_destPosFunctionValue = h_destPosFunctionValue;
		thrust::device_vector<float> d_destNegFunctionValue = h_destNegFunctionValue;
		thrust::device_vector<unsigned long> d_destCluster = h_destcluster;
		thrust::device_vector<unsigned long> d_destNumComb = h_destNumComb;
	
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		unsigned long* destClusterArray = thrust::raw_pointer_cast( &d_destCluster[0] );
		float* destPosImbArray = thrust::raw_pointer_cast( &d_destPosFunctionValue[0] );
		float* destNegImbArray = thrust::raw_pointer_cast( &d_destNegFunctionValue[0] );
		ulong* destNumCombArray = thrust::raw_pointer_cast( &d_destNumComb[0] );

		// size_t reqBlocksCount = (n * (nc - 1)) / threadsCount;
		// ushort blocksCount = (ushort)std::min((size_t)32768, reqBlocksCount);
		int blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
    	// printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsCount);
		// <<<blocksCount, threadsCount>>>
		simpleSearchKernel<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray, offsetArray, 
				clusterArray, funcArray, uint(n), uint(m), destClusterArray, destPosImbArray, destNegImbArray, 
				ulong(nc), ulong(numberOfChunks), firstImprovement, destNumCombArray, randomIndexArray);
		
		checkCudaErrors(cudaDeviceSynchronize());

		h_destcluster = d_destCluster;
		h_destPosFunctionValue = d_destPosFunctionValue;
		h_destNegFunctionValue = d_destNegFunctionValue;
		h_destNumComb = d_destNumComb;
		
		return true;
	}
		
}
