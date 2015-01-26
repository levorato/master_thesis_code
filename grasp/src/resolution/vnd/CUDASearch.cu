#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <stdio.h>
#include <curand_kernel.h>
#include <vector>

// CUDA facts:
//
// On devices of compute capability 2.x and beyond, 32-bit integer multiplication is natively supported,
// but 24-bit integer multiplication is not. __[u]mul24 is therefore implemented using multiple instructions
// and should not be used.
//
// Integer division and modulo operation are costly: below 20 instructions on devices of compute capability 2.x and
// higher. They can be replaced with bitwise operations in some cases: If n is a power of 2, (i/n) is equivalent to
// (i>>log2(n)) and (i%n) is equivalent to (i&(n-1)); the compiler will perform these conversions if n is literal.

// Always check return codes of CUDA calls for errors. Do not use __syncthreads() in conditional code unless the condition 
// is guaranteed to evaluate identically for all threads of each block. Run your program under cuda-memcheck to detect stray 
// memory accesses. If your kernel dies for larger problem sizes, it might exceed the runtime limit and trigger the watchdog timer.

// Thread block size should always be a multiple of 32, because kernels issue instructions in warps (32 threads).
#define BLOCK_SIZE 256.0

namespace clusteringgraph {
	
	/// CUDA Kernel to update the cluster array and the objective function value on 1-opt search.
	__global__ void updateClustering1opt(const int bestSrcVertex, const int destcluster, const float destFunctionValue, 
			ulong* clusterArray, float* funcArray, int n, uint* nc) {
		funcArray[0] = destFunctionValue;
		int previousCluster = clusterArray[bestSrcVertex];
		clusterArray[bestSrcVertex] = destcluster;
		
		// validate the current number of clusters in the solution
		if(destcluster >= nc[0]) {  // Option 1 (easier): destcluster == nc, i.e., bestSrcVertex is going to a standalone cluster
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
	
	/// CUDA Kernel to update the cluster array and the objective function value on 2-opt search.
	__global__ void updateClustering2opt(const int bestSrcVertex1, const int bestSrcVertex2, const int destcluster1, 
			const int destcluster2, const float destFunctionValue, 
			ulong* clusterArray, float* funcArray, int n, uint* nc) {
		funcArray[0] = destFunctionValue;
		// move the first vertex from k1 to k3
		int k1 = clusterArray[bestSrcVertex1];
		clusterArray[bestSrcVertex1] = destcluster1;
		int k3 = destcluster1;
		// move the second vertex from k2 to k4
		int k2 = clusterArray[bestSrcVertex2];
		clusterArray[bestSrcVertex2] = destcluster2;
		int k4 = destcluster2;
		uint initial_nc = nc[0];
		
		// --- FIRST SWAP
		// validate the current number of clusters in the solution
		if(k3 >= initial_nc) {  // Option 1 (easier): destcluster == nc, i.e., bestSrcVertex is going to a standalone cluster
			nc[0]++;  // just increment the number of clusters
		}
		// Option 2 (harder): bestSrcVertex's previous cluster no longer exists (removed its only element)
		// Check if there are other vertices in the previous cluster of bestSrcVertex
		bool found = false;
		for(int i = 0; i < n; i++) {
			if(clusterArray[i] == k1) {
				found = true;
				break;
			}
		}
		if(not found) {  // k1 is empty, has to be removed
			for(int i = 0; i < n; i++) {
				if(clusterArray[i] > k1) {
					clusterArray[i]--;
				}
			}
			nc[0]--;
		}
		
		// --- SECOND SWAP
		// validate the current number of clusters in the solution
		if(k4 >= initial_nc) {  // Option 1 (easier): destcluster == nc, i.e., bestSrcVertex is going to a standalone cluster
			nc[0]++;  // just increment the number of clusters
		}
		// Option 2 (harder): bestSrcVertex's previous cluster no longer exists (removed its only element)
		// Check if there are other vertices in the previous cluster of bestSrcVertex
		found = false;
		for(int i = 0; i < n; i++) {
			if(clusterArray[i] == k2) {
				found = true;
				break;
			}
		}
		if(not found) {  // k1 is empty, has to be removed
			for(int i = 0; i < n; i++) {
				if(clusterArray[i] > k2) {
					clusterArray[i]--;
				}
			}
			nc[0]--;
		}
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArrays2(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int n, uint* ncArray) {
		uint nc = ncArray[0];
		for(int i = 0; i < n; i++) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[k + (nc+1) * i] = 0.0;
            	vertexClusterNegSumArray[k + (nc+1) * i] = 0.0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					vertexClusterPosSumArray[clusterArray[j] + (nc+1) * i] += fabs(weight);
				} else {
					vertexClusterNegSumArray[clusterArray[j] + (nc+1) * i] += fabs(weight);
				}
			}
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArrays(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int n, uint* ncArray) {
		
		// TODO: fazer cada thread da GPU processar uma aresta do grafo
		int i = blockDim.x*blockIdx.x + threadIdx.x;
		uint nc = ncArray[0];
		// shared memory
		extern __shared__ long sc[];       // n longs
		long *s_cluster = sc;              // n longs
		
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads();
	    if(i < n) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[k + (nc+1) * i] = 0.0;
            	vertexClusterNegSumArray[k + (nc+1) * i] = 0.0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					vertexClusterPosSumArray[s_cluster[j] + (nc+1) * i] += fabs(weight);
				} else {
					vertexClusterNegSumArray[s_cluster[j] + (nc+1) * i] += fabs(weight);
				}
			}
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArraysDelta(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int n, uint* ncArray, int moved_vertex, int old_cluster) {
		int i = blockDim.x*blockIdx.x + threadIdx.x;
		uint nc = ncArray[0];
		// shared memory
		extern __shared__ long sc[];       // n longs
		long *s_cluster = sc;              // n longs
		
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads();
	    	if(i == moved_vertex) {
        		for(int k = 0; k <= nc; k++) {
            			vertexClusterPosSumArray[i*(nc+1) + k] = 0.0;
            			vertexClusterNegSumArray[i*(nc+1) + k] = 0.0;
            		}
            		// in/out-edges of vertex i
            		int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					//if(i != j) {  // out-edges of i
						vertexClusterPosSumArray[s_cluster[j] + (nc+1) * i] += fabs(weight);
					//} else {  // in-edges of i
						vertexClusterPosSumArray[old_cluster + (nc+1) * j] -= fabs(weight);
						vertexClusterPosSumArray[s_cluster[i] + (nc+1) * j] += fabs(weight);
					//}
				} else {
					//if(i != j) {
						vertexClusterNegSumArray[s_cluster[j] + (nc+1) * i] += fabs(weight);
					//} else {
						vertexClusterNegSumArray[old_cluster + (nc+1) * j] -= fabs(weight);
						vertexClusterNegSumArray[s_cluster[i] + (nc+1) * j] += fabs(weight);	
					//}
				}
			}
		}
	}

	/// CUDA kernel for simple 1-opt search.
	__global__ void simpleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		float* destImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, uint* randomIndexArray, float* positiveSumArray,
		float* negativeSumArray, int threadsPerBlock) {

		ulong idx = blockIdx.x * blockDim.x + threadIdx.x;
		ulong i  = idx / (nc + 1);  // kernel executes local search for vertex i
		ulong k2 = idx % (nc + 1);  // kernel executes local search for vertex i, moving it to cluster k2
		
		//if(idx < numberOfChunks) {
		if(i < n && k2 <= nc) {
			// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
	                // calculates the cost of removing vertex i from cluster1 and inserting into cluster2		
			// positiveSum = +(positiveSumArray[nc] - positiveSumArray[k2]) - (positiveSumArray[nc] - positiveSumArray[k1]);
			// vertex i is in cluster(k1)
                	ulong k1 = clusterArray[i];			
			float positiveSum = +(- positiveSumArray[i*(nc+1)+k2]) - (- positiveSumArray[i*(nc+1)+k1]);
			float negativeSum = -negativeSumArray[i*(nc+1)+k1] + negativeSumArray[i*(nc+1)+k2];
			
			// numberOfTestedCombinations++;
			// updates thread idx / vertex i from k1 to k2 imbalance result
			destImbArray[idx] = funcArray[0] + positiveSum + negativeSum;
			//if(destImbArray[idx] <= 0) { // TODO corrigir isso aqui!!!
			//	destImbArray[idx] = funcArray[0];
			//}
		}
	}
	
	/// CUDA kernel for 2-opt search.
	__global__ void doubleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		ulong* destClusterArray1, ulong* destClusterArray2, float* destImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, uint* randomIndexArray, float* vertexClusterPosSumArray,
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
		// int baseI = i*(nc+1);
		// int baseJ = j*(nc+1);
		int tbase = threadIdx.x * (nc+1);
		int baseI = i, baseJ = j;
		for(ulong k = 0; k <= nc; k++, baseI += n, baseJ += n) {
			s_clusterPosSumI[tbase + k] = vertexClusterPosSumArray[baseI];
			s_clusterNegSumI[tbase + k] = vertexClusterNegSumArray[baseI];
			s_clusterPosSumJ[tbase + k] = vertexClusterPosSumArray[baseJ];
			s_clusterNegSumJ[tbase + k] = vertexClusterNegSumArray[baseJ];
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
			float originalImbalance = funcArray[0];
			float bestImbValue = originalImbalance;
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
			uint k3 = 0; // randomIndexArray[idx];
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
							if(originalImbalance + positiveSum2 + negativeSum2 < bestImbValue) {  // improvement in imbalance
								bestImbValue = positiveSum2 + negativeSum2 + originalImbalance;
								bestDestCluster1 = k3;
								bestDestCluster2 = k4;
								if(firstImprovement) {
									destImbArray[idx] = bestImbValue;
									destClusterArray1[idx] = bestDestCluster1;
									destClusterArray2[idx] = bestDestCluster2;
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
			destImbArray[idx] = bestImbValue;
			destClusterArray1[idx] = bestDestCluster1;
			destClusterArray2[idx] = bestDestCluster2;
		}
	}

        // C / CPU function to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
        void updateVertexClusterSumArraysCPU(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest, thrust::host_vector<int>& h_numedges,
                        thrust::host_vector<int>& h_offset, thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_VertexClusterPosSum, thrust::host_vector<float>& h_VertexClusterNegSum, int n, int nc) {
                for(int i = 0; i < n; i++) {
	            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
        	    for(int k = 0; k <= nc; k++) {
                	h_VertexClusterPosSum[k + (nc+1) * i] = 0.0;
                	h_VertexClusterNegSum[k + (nc+1) * i] = 0.0;
         	    }
		    // in/out-edges of vertex i
            	    int offset = h_offset[i];
                    int numedges = h_numedges[i];
                    ulong count = offset + numedges;
                    for (ulong edgenum = offset; edgenum < count; edgenum++) {
                                int j = h_dest[edgenum];
                                float weight = h_weights[edgenum];
                                if(weight > 0) {
                                        h_VertexClusterPosSum[ h_mycluster[j] + (nc+1) * i] += fabs(weight);
                                } else {
                                        h_VertexClusterNegSum[ h_mycluster[j] + (nc+1) * i] += fabs(weight);
                                }
                    }
        	}
        }

	/// Runs the kernel for variable neighborhood descent (VND).
	bool runVNDKernel(thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m, ushort threadsCount, ulong& nc, ulong numberOfChunks, bool firstImprovement,
			thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
			thrust::host_vector<float>& h_VertexClusterNegSum, std::vector<uint>& sourceVertexList,
			std::vector<uint>& destinationClusterList, 
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
		thrust::device_vector<uint> d_randomIndex(n * (nc+1));
		thrust::device_vector<float> d_VertexClusterPosSum(n * (nc+1));
		thrust::device_vector<float> d_VertexClusterNegSum(n * (nc+1));
		thrust::device_vector<uint> d_nc(1);
		thrust::device_vector<float> d_destFunctionValue(numberOfChunks);
		thrust::device_vector<unsigned long> d_destCluster1(numberOfChunks);
		thrust::device_vector<unsigned long> d_destCluster2(numberOfChunks);	
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		uint* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		
		// number of clusters - changes every iteration of VND
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = nc;
		d_nc = h_nc;

		// VND loop
		int r = 1, iteration = 0;
		float bestImbalance = h_functionValue[0];
		sourceVertexList.clear();
		destinationClusterList.clear();
		int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
		updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
			offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, ncArray);
		//updateVertexClusterSumArrays2<<<1, 1>>>(weightArray, destArray, numArray,
                //      offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, ncArray);
		checkCudaErrors(cudaDeviceSynchronize());
		// updateVertexClusterSumArraysCPU(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_VertexClusterPosSum, h_VertexClusterNegSum, n, h_nc[0]);
                // updates GPU arrays with CPU result
                //d_VertexClusterPosSum = h_VertexClusterPosSum;
                //d_VertexClusterNegSum = h_VertexClusterNegSum;
		
		while (r <= l /* && (timeSpentSoFar + timeSpentOnLocalSearch < timeLimit)*/) {
			printf("*** Local search iteration %d, r = %d, nc = %ld\n", iteration, r, h_nc[0]);
			if(r == 1) {
				numberOfChunks = (h_nc[0] + 1) * n;
			} else {
				numberOfChunks = n * n;
			}
			// result / destination vectors
			d_destFunctionValue.resize(numberOfChunks);
			d_destCluster1.resize(numberOfChunks);
			d_destCluster2.resize(numberOfChunks);
			unsigned long* destClusterArray1 = thrust::raw_pointer_cast( &d_destCluster1[0] );
			unsigned long* destClusterArray2 = thrust::raw_pointer_cast( &d_destCluster2[0] );
			float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
			
			printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
			blocksPerGrid = ((n * (h_nc[0] + 1)) + threadsCount - 1) / threadsCount;
			
			if(r == 1) {
		    	// printf("CUDA kernel launch with %d blocks of %d threads, shmem = %d\n", blocksPerGrid, threadsCount, threadsCount*2*h_nc[0]*sizeof(float));
				simpleSearchKernel<<<blocksPerGrid, threadsCount>>>(clusterArray, funcArray, 
						uint(n), uint(m), destImbArray, 
						ulong(h_nc[0]), ulong(numberOfChunks), firstImprovement, randomIndexArray,
						vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
			} else {		
		    	// printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsCount);
				doubleSearchKernel<<<blocksPerGrid, threadsCount, threadsCount*4*(h_nc[0]+1)*sizeof(float)>>>(clusterArray, funcArray, 
						uint(n), uint(m), destClusterArray1, destClusterArray2, destImbArray, 
						ulong(h_nc[0]), ulong(numberOfChunks), firstImprovement, randomIndexArray,
						vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
			}
			checkCudaErrors(cudaDeviceSynchronize());
			// printf("Begin reduce...\n");
			thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.begin()+numberOfChunks);
			float min_val = *iter;
						
			if(min_val < bestImbalance) {  // improvement in imbalance
				bestImbalance = min_val;
				// determines the position of the best improvement found in the result vector
				uint position = iter - d_destFunctionValue.begin();
				
				// post-processing of improved solution
				// printf("Begin post-process...\n");
				if(r == 1) {
					int resultIdx = position;
					thrust::device_vector<unsigned long>::iterator it = d_destCluster1.begin() + position;
					// for(int i = 0; i < position; i++, it++);
					ulong destcluster = resultIdx % (h_nc[0] + 1);  // ALTERADO AQUI o K2
					ulong bestSrcVertex = resultIdx / (h_nc[0] + 1);
					ulong sourceCluster = d_mycluster[bestSrcVertex];
					printf("Idx = %d: The best src vertex is %d to cluster %d with I(P) = %.2f\n", resultIdx, bestSrcVertex, destcluster, destFunctionValue);
					if(destFunctionValue < 0) {  printf("WARNING: I(P) < 0 !!!\n");  }
					updateClustering1opt<<< 1, 1 >>>(bestSrcVertex, destcluster, bestImbalance, clusterArray, funcArray, n, ncArray);
					sourceVertexList.push_back(bestSrcVertex);
					destinationClusterList.push_back(destcluster);
					int old_nc = h_nc[0];
	                                h_nc = d_nc;
        	                        // h_functionValue = d_functionValue;
                	                // printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
                                	if(old_nc != h_nc[0]) {  // the number of clusters in the solution changed
                                        	d_VertexClusterPosSum.resize(n * (h_nc[0]+1));
                                        	d_VertexClusterNegSum.resize(n * (h_nc[0]+1));
						h_VertexClusterPosSum.resize(n * (h_nc[0]+1));
                                                h_VertexClusterNegSum.resize(n * (h_nc[0]+1));
                                        	vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
                                        	vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
						printf("Number of clusters has changed.\n");
                                	}

					h_mycluster = d_mycluster; // retrieves new cluster configuration from GPU

					blocksPerGrid = (n + threadsCount - 1) / threadsCount;
					// updateVertexClusterSumArrays2<<<1, 1>>>(weightArray, destArray, numArray,
                      			//	offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, ncArray);
					updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray, offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, ncArray);
					// updateVertexClusterSumArraysDelta<<<blocksPerGrid, threadsCount, n*sizeof(long)+2*(nc+1)*sizeof(float)>>>(weightArray, destArray, numArray, offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, n, ncArray, bestSrcVertex, sourceCluster);
					checkCudaErrors(cudaDeviceSynchronize());
					// int moved_vertex, int old_cluster	
					// updateVertexClusterSumArraysCPU(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_VertexClusterPosSum, h_VertexClusterNegSum, n, h_nc[0]);
					// updates GPU arrays with CPU result
					// d_VertexClusterPosSum = h_VertexClusterPosSum;
					// d_VertexClusterNegSum = h_VertexClusterNegSum;
				} 
				checkCudaErrors(cudaDeviceSynchronize());
				// printf("Preparing new VND loop...\n");
				// h_functionValue = d_functionValue;
				// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
				r = 1;
				// if(bestImbalance < 0)  break;
			} else {  // no better result found in neighborhood
				r++;
			}
			iteration++;
		}
		h_mycluster = d_mycluster;
		nc = h_nc[0];
		destFunctionValue = bestImbalance;
		// printf("Exiting with nc = %ld\n", h_nc[0]);
		return true;
	}
	
	/// Runs a kernel for 2-opt search.
	bool run2optSearchKernel(thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue,
			ulong n, ulong m, thrust::host_vector<unsigned long>& h_destcluster1, 
			thrust::host_vector<unsigned long>& h_destcluster2, 
			thrust::host_vector<float>& h_destFunctionValue, ushort threadsCount, ulong nc, ulong numberOfChunks,
			bool firstImprovement, 
			thrust::host_vector<uint>& h_randomIndex, thrust::host_vector<float>& h_VertexClusterPosSum,
			thrust::host_vector<float>& h_VertexClusterNegSum) {

		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<uint> d_randomIndex = h_randomIndex;
		thrust::device_vector<float> d_VertexClusterPosSum = h_VertexClusterPosSum;
		thrust::device_vector<float> d_VertexClusterNegSum = h_VertexClusterNegSum;
		// destination vectors
		thrust::device_vector<float> d_destFunctionValue = h_destFunctionValue;
		thrust::device_vector<unsigned long> d_destCluster1 = h_destcluster1;
		thrust::device_vector<unsigned long> d_destCluster2 = h_destcluster2;
	
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		unsigned long* destClusterArray1 = thrust::raw_pointer_cast( &d_destCluster1[0] );
		unsigned long* destClusterArray2 = thrust::raw_pointer_cast( &d_destCluster2[0] );
		float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
		
		int blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
    	// printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsCount);
		doubleSearchKernel<<<blocksPerGrid, threadsCount, threadsCount*4*(nc+1)*sizeof(float)>>>(clusterArray, funcArray, 
				uint(n), uint(m), destClusterArray1, destClusterArray2, destImbArray, 
				ulong(nc), ulong(numberOfChunks), firstImprovement, randomIndexArray,
				vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
		
		checkCudaErrors(cudaDeviceSynchronize());

		h_destcluster1 = d_destCluster1;
		h_destcluster2 = d_destCluster2;
		h_destFunctionValue = d_destFunctionValue;
		
		return true;
	}
}
