#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
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

#define BLOCK_SIZE 512

namespace clusteringgraph {
	
	/// CUDA Kernel to update the cluster array and the objective function value.
	__global__ void updateClustering(const int bestSrcVertex, const int destcluster, const float destFunctionValue, 
			ulong* clusterArray, float* funcArray, int n, uint* nc) {
		funcArray[0] = destFunctionValue;
		int previousCluster = clusterArray[bestSrcVertex];
		clusterArray[bestSrcVertex] = destcluster;
		
		// validate the current number of clusters in the solution
		if(destcluster == nc[0]) {  // Option 1 (easier): destcluster == nc, i.e., bestSrcVertex is going to a standalone cluster
			nc[0]++;  // just increment the number of clusters
		}
		// Option 2 (harder): bestSrcVertex's previous cluster no longer exists (removed its only element)
		// Check if there are other vertices in the previous cluster of bestSrcVertex
		bool found = false;
		for(int i = 0; i < n; i++) {
			if(clusterArray[i] == previousCluster) {
				found = true;
				break;
			}
		}
		if(not found) {  // previousCluster is empty, has to be removed
			for(int i = 0; i < n; i++) {
				if(clusterArray[i] > previousCluster) {
					clusterArray[i]--;
				}
			}
			nc[0]--;
		}
		
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArrays2(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int n, int nc) {
		for(int i = 0; i < n; i++) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[i * (nc+1) + k] = 0.0;
            	vertexClusterNegSumArray[i * (nc+1) + k] = 0.0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					vertexClusterPosSumArray[i * (nc+1) + clusterArray[j]] += fabs(weight);
				} else {
					vertexClusterNegSumArray[i * (nc+1) + clusterArray[j]] += fabs(weight);
				}
			}
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArrays(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int n, int nc) {
		int i = blockDim.x*blockIdx.x + threadIdx.x;
		// shared memory
		extern __shared__ long sc[];       // n longs
		long *s_cluster = sc;              // n longs
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
	    if(i < n) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[i * (nc+1) + k] = 0.0;
            	vertexClusterNegSumArray[i * (nc+1) + k] = 0.0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					vertexClusterPosSumArray[i * (nc+1) + s_cluster[j]] += fabs(weight);
				} else {
					vertexClusterNegSumArray[i * (nc+1) + s_cluster[j]] += fabs(weight);
				}
			}
        }
	}

	/// CUDA kernel for simple 1-opt search.
	__global__ void simpleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		ulong* destClusterArray, float* destImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, uint* randomIndexArray, float* vertexClusterPosSumArray,
		float* vertexClusterNegSumArray, int threadsPerBlock) {

		uint idx = blockIdx.x * blockDim.x + threadIdx.x;
		uint i = idx % n;  // kernel executes local search for vertex i
		uint part = idx / n;  // kernel executes local search for vertex i, for the 'part' partition
		uint nparts = numberOfChunks / n;
		// shared memory
		extern __shared__ float s[];       // 2 * blockDim.x * nc floats
		float *s_clusterPosSum = s;        // blockDim.x * nc floats
		float *s_clusterNegSum = &s_clusterPosSum[threadsPerBlock * (nc+1)];   // blockDim.x * nc floats
		int base = i*(nc+1);
		int tbase = threadIdx.x * (nc+1);
		int tbase2 = threadsPerBlock * (nc+1) + threadIdx.x * (nc+1); // blockDim.x * nc + threadIdx.x * nc;
		for(int k = 0; k <= nc; k++) {
			s_clusterPosSum[tbase + k] = vertexClusterPosSumArray[base + k];
			s_clusterPosSum[tbase2 + k] = vertexClusterNegSumArray[base + k];
		}
		// ensure that all threads have loaded their values into
		// shared memory; otherwise, one thread might be computing
		// on unitialized data.
		// __syncthreads();
		
		if(idx < numberOfChunks) {
			ulong numberOfTestedCombinations = 0;
			// vertex i is in cluster(k1)
			ulong k1 = clusterArray[i];
			uint bestDestCluster = k1;
			// stores the best imbalance found so far
			float originalImbalance = funcArray[0];
			float bestImbalance = originalImbalance;
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
				float negativeSumK1 = -s_clusterNegSum[tbase+k1];
                float positiveSumK1 = 0.0;
                for(uint k = 0; k < nc; k++) {
                        if(k != k1) {
                                positiveSumK1 -= s_clusterPosSum[tbase+k];
                        }
                }
				// Random initial vertex
				uint k2 = randomIndexArray[idx];
				for(uint countK2 = 0; countK2 < chunkSize; countK2++) {  // cluster(k2)
					if(k2 != k1) {
						// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
						float negativeSum = negativeSumK1 + s_clusterNegSum[tbase+k2];
						float positiveSum = positiveSumK1;
						for(uint k = 0; k < nc; k++) {
							if(k != k2) {
								positiveSum += s_clusterPosSum[tbase+k];
							}
						}
						numberOfTestedCombinations++;
						if(originalImbalance + positiveSum + negativeSum < bestImbalance) {  // improvement in imbalance
							bestImbalance = positiveSum + negativeSum + originalImbalance;
							bestDestCluster = k2;
							if(firstImprovement) {
								destImbArray[idx] = bestImbalance;
								destClusterArray[idx] = bestDestCluster;
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
				destImbArray[idx] = bestImbalance;
				destClusterArray[idx] = bestDestCluster;
			} else {
				destImbArray[idx] = originalImbalance;
				destClusterArray[idx] = k1;
			}
		}
	}
	
		/// CUDA kernel for 2-opt search.
	__global__ void doubleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		ulong* destClusterArray1, ulong* destClusterArray2, float* destPosImbArray, float* destNegImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, ulong* destNumCombArray, uint* randomIndexArray, float* vertexClusterPosSumArray,
		float* vertexClusterNegSumArray, int threadsPerBlock) {

		uint idx = blockIdx.x * blockDim.x + threadIdx.x;
		uint i = idx % n;  // kernel executes local search for vertex i
		uint j = idx / n;  // kernel executes local search for vertex i, and vertex j
		// shared memory
		extern __shared__ float s[];       										// 4 * blockDim.x * nc+1 floats
		float *s_clusterPosSumI = s;        								 	// blockDim.x * nc+1 floats
		float *s_clusterNegSumI = &s_clusterPosSumI[threadsPerBlock * (nc+1)];  // blockDim.x * nc+1 floats
		float *s_clusterPosSumJ = &s_clusterNegSumI[threadsPerBlock * (nc+1)];  // blockDim.x * nc+1 floats
		float *s_clusterNegSumJ = &s_clusterPosSumJ[threadsPerBlock * (nc+1)];	// blockDim.x * nc+1 floats
		int baseI = i*(nc+1);
		int baseJ = j*(nc+1);
		int tbase = threadIdx.x * (nc+1);
		for(ulong k = 0; k <= nc; k++) {
			s_clusterPosSumI[tbase + k] = vertexClusterPosSumArray[baseI + k];
			s_clusterNegSumI[tbase + k] = vertexClusterNegSumArray[baseI + k];
			s_clusterPosSumJ[tbase + k] = vertexClusterPosSumArray[baseJ + k];
			s_clusterNegSumJ[tbase + k] = vertexClusterNegSumArray[baseJ + k];
		}
		// ensure that all threads have loaded their values into
		// shared memory; otherwise, one thread might be computing
		// on unitialized data.
		__syncthreads();
		
		if(idx < numberOfChunks) {
			ulong numberOfTestedCombinations = 0;
			// vertex i is in cluster(k1)
			ulong k1 = clusterArray[i];
			// vertex j is in cluster(k2)
			ulong k2 = clusterArray[j];
			// best destination clusters
			uint bestDestCluster1 = k1;
			uint bestDestCluster2 = k2;
			// stores the best imbalance found so far
			float bestValuePos, bestValueNeg;
			float originalPosImbalance = funcArray[0];
			float originalNegImbalance = funcArray[1];
			bestValuePos = originalPosImbalance;
			bestValueNeg = originalNegImbalance;
			// Option 1: vertex i is moved from k1 to another existing cluster k3 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k3 = nc)
			// REMOVAL of vertex i from cluster k1 -> avoids recalculating 
			//   the same thing in the destination cluster loop
			float negativeSumK1 = -(s_clusterNegSumI[tbase+k1] + s_clusterNegSumJ[tbase+k2]);
            float positiveSumK1 = 0.0;
            for(uint k = 0; k < nc; k++) {
                    if(k != k1) {
                            positiveSumK1 -= s_clusterPosSumI[tbase+k];
                    }
                    if(k != k2) {
                            positiveSumK1 -= s_clusterPosSumJ[tbase+k];
                    }
            }
			// Random initial vertex
			uint k3 = randomIndexArray[idx];
			for(uint countK3 = 0; countK3 <= nc; countK3++) {  // cluster(k3)
				if(k3 != k1) {
					// calculates the cost of removing vertex i from cluster k1 and inserting into cluster k3
					float negativeSum = negativeSumK1 + s_clusterNegSumI[tbase+k3];
					float positiveSum = positiveSumK1;
					for(uint k = 0; k < nc; k++) {
						if(k != k3) {
							positiveSum += s_clusterPosSumI[tbase+k];
						}
					}
					for(uint k4 = k3 + 1; k4 <= nc; k4++) {  // cluster(k4)
						if(k4 != k2) {
							// calculates the cost of removing vertex j from cluster k2 and inserting into cluster k4
							float negativeSum2 = negativeSum + s_clusterNegSumJ[tbase+k4];
							float positiveSum2 = positiveSum;
							for(uint k = 0; k < nc; k++) {
								if(k != k4) {
									positiveSum2 += s_clusterPosSumJ[tbase+k];
								}
							}
							numberOfTestedCombinations++;
							if(originalPosImbalance + originalNegImbalance + positiveSum2 + negativeSum2 < bestValuePos + bestValueNeg) {  // improvement in imbalance
								bestValuePos = positiveSum2 + originalPosImbalance;
								bestValueNeg = negativeSum2 + originalNegImbalance;
								bestDestCluster1 = k3;
								bestDestCluster2 = k4;
								if(firstImprovement) {
									destPosImbArray[idx] = bestValuePos;
									destNegImbArray[idx] = bestValueNeg;
									destClusterArray1[idx] = bestDestCluster1;
									destClusterArray2[idx] = bestDestCluster2;
									destNumCombArray[idx] = numberOfTestedCombinations;
									return;
								}
							}
						}
					}
				}
				// loop increment rule
				k3++;
				if(k3 > nc) {
					k3 = 0;
				}
			}
			// updates thread idx / vertex i imbalance result
			destPosImbArray[idx] = bestValuePos;
			destNegImbArray[idx] = bestValueNeg;
			destClusterArray1[idx] = bestDestCluster1;
			destClusterArray2[idx] = bestDestCluster2;
			destNumCombArray[idx] = numberOfTestedCombinations;
		}
	}
	
	/// Runs a kernel for 1-opt search.
	bool run1optSearchKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m, ushort threadsCount, ulong nc, ulong numberOfChunks, bool firstImprovement,
			thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
			thrust::host_vector<float>& h_VertexClusterNegSum, uint& bestSrcVertex, uint& destcluster, 
			float& destFunctionValue, const long& timeSpentSoFar, const unsigned int& l) {

		// declaration and initialization of variables
		// Graph data - read only
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		// current clustering data
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<uint> d_randomIndex = h_randomIndex;
		thrust::device_vector<float> d_VertexClusterPosSum = h_VertexClusterPosSum;
		thrust::device_vector<float> d_VertexClusterNegSum = h_VertexClusterNegSum;
		thrust::device_vector<uint> d_nc(1);
		// result / destination vectors
		thrust::device_vector<float> d_destFunctionValue(numberOfChunks);
		thrust::device_vector<unsigned long> d_destCluster(numberOfChunks);
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		unsigned long* destClusterArray = thrust::raw_pointer_cast( &d_destCluster[0] );
		float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
		uint* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		
		// number of clusters - changes every iteration of VND
		thrust::host_vector<uint> h_nc(1);
		h_nc[0] = nc;
		
		// VND loop
		int r = 1, iteration = 0;
		double timeSpentOnLocalSearch = 0.0;
		float bestImbalance = h_functionValue[0];
		
		while (r <= l and iteration < 20 /* && (timeSpentSoFar + timeSpentOnLocalSearch < timeLimit)*/) {
			// printf("*** Local search iteration %d\n", iteration);
			// 0. Triggers local processing time calculation
			// boost::timer::cpu_timer timer;
			// timer.start();
			// boost::timer::cpu_times start_time = timer.elapsed();
			
			int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
			updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
					offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, h_nc[0]);
			checkCudaErrors(cudaDeviceSynchronize());
			
			// printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
			blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
	    	// printf("CUDA kernel launch with %d blocks of %d threads, shmem = %d\n", blocksPerGrid, threadsCount, threadsCount*2*h_nc[0]*sizeof(float));
			simpleSearchKernel<<<blocksPerGrid, threadsCount, threadsCount*2*(h_nc[0]+1)*sizeof(float)>>>(clusterArray, funcArray, 
					uint(n), uint(m), destClusterArray, destImbArray, 
					ulong(h_nc[0]), ulong(numberOfChunks), firstImprovement, randomIndexArray,
					vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
			checkCudaErrors(cudaDeviceSynchronize());
			
			thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.end());
			float min_val = *iter;
			
			if(min_val < bestImbalance) {  // improvement in imbalance
				bestImbalance = min_val;
				destFunctionValue = min_val;
				uint position = iter - d_destFunctionValue.begin();
				int resultIdx = position;
				thrust::device_vector<unsigned long>::iterator it = d_destCluster.begin();
				for(int i = 0; i < position; i++, it++);
				destcluster = *it;
				bestSrcVertex = resultIdx % n;
				printf("The best src vertex is %d to cluster %d with I(P) = %.2f\n", bestSrcVertex, destcluster, destFunctionValue);
				
				// post-processing of improved solution
				int old_nc = h_nc[0];
				d_nc = h_nc;
				// printf("Before: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
				updateClustering<<< 1, 1 >>>(bestSrcVertex, destcluster, destFunctionValue, clusterArray, funcArray, n, ncArray);
				checkCudaErrors(cudaDeviceSynchronize());

				d_nc = h_nc;
				h_functionValue = d_functionValue;
				// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
				if(old_nc != h_nc[0]) {  // the number of clusters in the solution changed
					d_VertexClusterPosSum.resize(n * (h_nc[0]+1));
					d_VertexClusterNegSum.resize(n * (h_nc[0]+1));
					h_VertexClusterPosSum.resize(n * (h_nc[0]+1));
					h_VertexClusterNegSum.resize(n * (h_nc[0]+1));
					vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
					vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
				}		
				
				blocksPerGrid = (n + threadsCount - 1) / threadsCount;
				
				updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
					offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, h_nc[0]);
				checkCudaErrors(cudaDeviceSynchronize());
				r = 1;
			} else {  // no better result found in neighborhood
				r++;
			}
			iteration++;
	
			// => Finally: Stops the timer and stores the elapsed time
			// timer.stop();
			// boost::timer::cpu_times end_time = timer.elapsed();
			// timeSpentOnLocalSearch += (end_time.wall - start_time.wall) / double(1000000000);
		}
		h_mycluster = d_mycluster;
		
		return true;
	}
	
	/// Runs a kernel for 2-opt search.
	bool run2optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m, thrust::host_vector<unsigned long>& h_destcluster1, 
			thrust::host_vector<unsigned long>& h_destcluster2, thrust::host_vector<float>& h_destPosFunctionValue,
			thrust::host_vector<float>& h_destNegFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
			bool firstImprovement, thrust::host_vector<unsigned long>& h_destNumComb,
			thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
			thrust::host_vector<float>& h_VertexClusterNegSum) {

		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<uint> d_randomIndex = h_randomIndex;
		thrust::device_vector<float> d_VertexClusterPosSum = h_VertexClusterPosSum;
		thrust::device_vector<float> d_VertexClusterNegSum = h_VertexClusterNegSum;
		// destination vectors
		thrust::device_vector<float> d_destPosFunctionValue = h_destPosFunctionValue;
		thrust::device_vector<float> d_destNegFunctionValue = h_destNegFunctionValue;
		thrust::device_vector<unsigned long> d_destCluster1 = h_destcluster1;
		thrust::device_vector<unsigned long> d_destCluster2 = h_destcluster2;
		thrust::device_vector<unsigned long> d_destNumComb = h_destNumComb;
	
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		unsigned long* destClusterArray1 = thrust::raw_pointer_cast( &d_destCluster1[0] );
		unsigned long* destClusterArray2 = thrust::raw_pointer_cast( &d_destCluster2[0] );
		float* destPosImbArray = thrust::raw_pointer_cast( &d_destPosFunctionValue[0] );
		float* destNegImbArray = thrust::raw_pointer_cast( &d_destNegFunctionValue[0] );
		ulong* destNumCombArray = thrust::raw_pointer_cast( &d_destNumComb[0] );

		int blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
    	// printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsCount);
		doubleSearchKernel<<<blocksPerGrid, threadsCount, threadsCount*4*(nc+1)*sizeof(float)>>>(clusterArray, funcArray, 
				uint(n), uint(m), destClusterArray1, destClusterArray2, destPosImbArray, destNegImbArray, 
				ulong(nc), ulong(numberOfChunks), firstImprovement, destNumCombArray, randomIndexArray,
				vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
		
		checkCudaErrors(cudaDeviceSynchronize());

		h_destcluster1 = d_destCluster1;
		h_destcluster2 = d_destCluster2;
		h_destPosFunctionValue = d_destPosFunctionValue;
		h_destNegFunctionValue = d_destNegFunctionValue;
		h_destNumComb = d_destNumComb;
		
		return true;
	}
}
