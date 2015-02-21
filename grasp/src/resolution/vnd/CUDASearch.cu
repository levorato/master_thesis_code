#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include "problem/include/ClusteringProblem.h"
#include "graph/include/Graph.h"
#include "util/include/RandomUtil.h"
#include "../construction/include/VertexSet.h"

#include <boost/timer/timer.hpp>
#include <iostream>

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <stdio.h>
#include <curand_kernel.h>
#include <vector>

#define EPS 10e-6

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

using namespace problem;
using namespace resolution::construction;
using namespace clusteringgraph;
using namespace std;
using namespace boost;
	
	/// CUDA Kernel to update the cluster array after constructive phase iteration.
	__global__ void updateConstructClustering(const int bestSrcVertex, const int destcluster, const float destFunctionValue, 
			ulong* clusterArray, float* funcArray, ulong n, ulong nc) {
		funcArray[0] = destFunctionValue;
		clusterArray[bestSrcVertex] = destcluster;
		if(destcluster == nc) {  // if the number of clusters is to increase
			for(ulong v = 0; v < n; v++) {  // increments the cluster number of every vertex that does not have a cluster yet
				if((clusterArray[v] == nc) && (v != bestSrcVertex)) {
					clusterArray[v]++;
				}
			}
		}
	}
	
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
            	vertexClusterPosSumArray[k * n + i] = 0.0;
            	vertexClusterNegSumArray[k * n + i] = 0.0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(weight > 0) {
					vertexClusterPosSumArray[clusterArray[j] * n + i] += fabs(weight);
				} else {
					vertexClusterNegSumArray[clusterArray[j] * n + i] += fabs(weight);
				}
			}
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArrays(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int* isNeighborClusterArray, int n, uint* ncArray) {
		
		int i = blockDim.x*blockIdx.x + threadIdx.x;
		uint nc = ncArray[0];
		// shared memory
		extern __shared__ int sc[];       // n longs
		int *s_cluster = sc;              // n longs
		
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads();
	    if(i < n) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[k * n + i] = 0.0;
            	vertexClusterNegSumArray[k * n + i] = 0.0;
            	isNeighborClusterArray[k * n + i] = 0;
            }
            // in/out-edges of vertex i
            int offset = offsetArray[i];
			int numedges = numArray[i];
			ulong count = offset + numedges;
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				int j = destArray[edgenum];
				float weight = weightArray[edgenum];
				if(s_cluster[j] >= 0) {
					if(s_cluster[i] != s_cluster[j]) {  // different cluster
						isNeighborClusterArray[s_cluster[j] * n + i]++;  // vertex i now has external connection to cluster kj
						isNeighborClusterArray[s_cluster[i] * n + j]++;  // vertex j now has external connection to cluster ki
					}
					if(weight > 0) {
						vertexClusterPosSumArray[s_cluster[j] * n + i] += fabs(weight);
					} else {
						vertexClusterNegSumArray[s_cluster[j] * n + i] += fabs(weight);
					}
				}
			}
			// preenche a possibilidade de se mover o vertice i para um cluster novo (k = nc)
			isNeighborClusterArray[nc * n + i] = 1;
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArraysDelta(const float* weightArray, const int* destArray, const int* numArray,
			const int* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int* isNeighborClusterArray, int n, uint* old_ncArray, uint* ncArray, int i, int k1, int k2) {  // vertex i is being moved from cluster k1 to k2
		uint nc = old_ncArray[0];
		uint new_nc = ncArray[0];
		// shared memory
		extern __shared__ int sc[];       // n longs
		int *s_cluster = sc;              // n longs

		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads();

		/*  THIS PEACE OF CODE IS NEEEDED ONLY WHEN USING THE CUDA CONSTRUCTIVE KERNEL - DISABLED */
		if(new_nc > nc) {  // vertex i is being moved to a new cluster
			// move a fileira correspondente ao cluster k = nc na matriz de soma, shiftando os dados para a direita (nc + 1)
			for(int v = 0; v < n; v++) {
				vertexClusterPosSumArray[(new_nc) * n + v] = vertexClusterPosSumArray[(nc) * n + v];
				vertexClusterNegSumArray[(new_nc) * n + v] = vertexClusterNegSumArray[(nc) * n + v];
			}
			// zera a fileira movida anteriormente
			for(int v = 0; v < n; v++) {
				vertexClusterPosSumArray[(nc) * n + v] = 0.0;
				vertexClusterNegSumArray[(nc) * n + v] = 0.0;
			}
			// preenche a possibilidade de se mover todos os vertices para um cluster novo (k = new_nc)
			for(int v = 0; v < n; v++) {
				isNeighborClusterArray[v + new_nc * n] = 1;
				isNeighborClusterArray[v + nc * n] = 0;
			}
		}

		// in/out-edges of vertex i
		int offset = offsetArray[i];
		int numedges = numArray[i];
		ulong count = offset + numedges;
		// isNeighborClusterArray[i+k1*n] = 1;  // vertex i now has external connection to cluster k1
		// for(ulong k = 0; k < nc; k++) {
		//  	isNeighborClusterArray[i + k * n] = 0;
		// }
		for (ulong edgenum = offset; edgenum < count; edgenum++) {
			int j = destArray[edgenum];
			float weight = weightArray[edgenum];
			isNeighborClusterArray[j + k1 * n] -= 2;
			if(s_cluster[j] != k2) {
				isNeighborClusterArray[i + s_cluster[j] * n]++;  // vertex i now has external connection to cluster kj
				isNeighborClusterArray[j + k2 * n]++;  // vertex j now has external connection to cluster k2
			} else {  // cluster[i] == cluster[j] == k2
				isNeighborClusterArray[i + s_cluster[j] * n] = 0;  // vertex i now has NO external connections to cluster kj
				isNeighborClusterArray[j + k2 * n] = 0;  // vertex j now has NO external connections to cluster k2
			}
			if(weight > 0) {
				vertexClusterPosSumArray[j + k1 * n] -= fabs(weight);
				vertexClusterPosSumArray[j + k2 * n] += fabs(weight);
			} else {
				vertexClusterNegSumArray[j + k1 * n] -= fabs(weight);
				vertexClusterNegSumArray[j + k2 * n] += fabs(weight);
			}
		}

		if(new_nc < nc) {
			// remove a fileira correspondente ao cluster k1 na matriz de soma, shiftando os dados para a esquerda
			for(int k = k1 + 1; k <= nc; k++) {
				for(int v = 0; v < n; v++) {
					vertexClusterPosSumArray[(k - 1) * n + v] = vertexClusterPosSumArray[(k) * n + v];
					vertexClusterNegSumArray[(k - 1) * n + v] = vertexClusterNegSumArray[(k) * n + v];
					isNeighborClusterArray[(k - 1) * n + v] = isNeighborClusterArray[(k) * n + v];
				}
			}
		}
	}

	/// CUDA kernel for simple 1-opt search.
	__global__ void simpleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		float* destImbArray, ulong nc, ulong numberOfChunks, float* positiveSumArray,
		float* negativeSumArray, int* isNeighborClusterArray) {

		ulong idx = blockIdx.x * blockDim.x + threadIdx.x;
		ulong i  = idx % n;  // kernel executes local search for vertex i
		ulong k2 = idx / n;  // kernel executes local search for vertex i, moving it to cluster k2
		
		if(i < n && k2 <= nc) {
			ulong k1 = clusterArray[i];
			if( (k1 != k2) && (isNeighborClusterArray[i+k2*n] > 0) ) {
				// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
				// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
						// calculates the cost of removing vertex i from cluster1 and inserting into cluster2
				// positiveSum = +(positiveSumArray[nc] - positiveSumArray[k2]) - (positiveSumArray[nc] - positiveSumArray[k1]);
				// vertex i is in cluster(k1)

				float positiveSum = +(- positiveSumArray[i+k2*n]) - (- positiveSumArray[i+k1*n]);
				float negativeSum = -negativeSumArray[i+k1*n] + negativeSumArray[i+k2*n];

				// updates thread idx / vertex i from k1 to k2 imbalance result
				destImbArray[idx] = funcArray[0] + positiveSum + negativeSum;
			} else {
				destImbArray[idx] = funcArray[0];
			}
		}
	}

	/// CUDA kernel for simple construct clustering (addition of a vertex into all possible clusters).
	__global__ void simpleConstructKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		float* destImbArray, ulong nc, float* positiveSumArray, 
		float* negativeSumArray, int threadsPerBlock, ulong v) {

		ulong idx = blockIdx.x * blockDim.x + threadIdx.x;
		// kernel executes local search for vertex i
		ulong k2 = idx;  // kernel executes local search for vertex i, moving it to cluster k2
		
		if(k2 <= nc) {
			// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
			// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
			// vertex i is in no cluster
			/*
			 * Na inserção em cluster novo, contar apenas as relações externas entre o vértice i e os vértices
			 * que estão sem cluster. Não incluir as relações internas quando k = nc.
			 * Quando o vértice i for inserido em cluster existente, contar as relações internas a k2 negativas,
			 * bem como as relações externas a k2 positivas com i.
			 */
			destImbArray[k2] = funcArray[0];
			if(k2 < nc) {
				destImbArray[k2] += negativeSumArray[v + k2 * n];
			}
			for(int k = 0; k < nc; k++) {
				if(k != k2) {
					destImbArray[k2] += positiveSumArray[v + k * n];
				}
			}
			destImbArray[k2] += positiveSumArray[v + nc * n];
		}
	}
	
	__global__ void shuffleBestResult1opt(const float* destImbArray, uint imbArraySize, float bestImbalance, 
						uint randomInitialIndex, uint* resultIndexArray, float* resultValueArray) {
		// Starts the search for the bestImbalance element from the randomInitialIndex provided as parameter
		for(uint count = 0, i = randomInitialIndex; count < imbArraySize; i = (i + 1) % imbArraySize, count++) {
			if((destImbArray[i] - bestImbalance < EPS) && (fabs(destImbArray[i] - bestImbalance) > EPS)){   // (destImbArray[i] < bestImbalance)
				resultIndexArray[0] = i;
				resultValueArray[0] = destImbArray[i];
				break;
			}
		}
	}

	/// CUDA kernel for 2-opt search.
	__global__ void doubleSearchKernel(const ulong* clusterArray, const float* funcArray, uint n, uint m,
		ulong* destClusterArray1, ulong* destClusterArray2, float* destImbArray, ulong nc, ulong numberOfChunks,
		bool firstImprovement, float* vertexClusterPosSumArray,
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
                	h_VertexClusterPosSum[k * n + i] = 0.0;
                	h_VertexClusterNegSum[k * n + i] = 0.0;
         	    }
		    // in/out-edges of vertex i
            	    int offset = h_offset[i];
                    int numedges = h_numedges[i];
                    ulong count = offset + numedges;
                    for (ulong edgenum = offset; edgenum < count; edgenum++) {
                                int j = h_dest[edgenum];
                                float weight = h_weights[edgenum];
                                if(weight > 0) {
                                        h_VertexClusterPosSum[ h_mycluster[j] * n + i] += fabs(weight);
                                } else {
                                        h_VertexClusterNegSum[ h_mycluster[j] * n + i] += fabs(weight);
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
			thrust::host_vector<float>& h_VertexClusterNegSum, thrust::host_vector<uint> &h_neighbor_cluster,
			std::vector<uint>& sourceVertexList, std::vector<uint>& destinationClusterList,
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
		thrust::device_vector<uint> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue;
		// Result arrays for imbalance reduction
		thrust::device_vector<uint> d_result_index(1);
		thrust::device_vector<float> d_result_value(1);
		
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
		uint* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		uint* resultIndexArray = thrust::raw_pointer_cast( &d_result_index[0] );
		float* resultValueArray = thrust::raw_pointer_cast( &d_result_value[0] );
		
		// number of clusters - changes every iteration of VND
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = nc;
		d_nc = h_nc;

		// VND loop
		int r = 1, iteration = 0;
		float bestImbalance = h_functionValue[0];
		sourceVertexList.clear();
		destinationClusterList.clear();
		thrust::device_vector<int> d_neighbor_cluster(n * (h_nc[0]+1), 0);
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );

		int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
		updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(int)>>>(weightArray, destArray, numArray,
			offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
		checkCudaErrors(cudaDeviceSynchronize());
		
		while (r <= l) {
			// printf("*** Local search iteration %d, r = %d, nc = %ld\n", iteration, r, h_nc[0]);
			// (r == 1)
			numberOfChunks = (h_nc[0] + 1) * n;
			// result / destination vectors: every element is initialized with the current best imbalance
			d_destFunctionValue = thrust::device_vector<float>(numberOfChunks, h_functionValue[0]);
			float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
			
			// printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
			blocksPerGrid = ((n * (h_nc[0] + 1)) + threadsCount - 1) / threadsCount;
			
			// (r == 1)
			// printf("CUDA kernel launch with %d blocks of %d threads, shmem = 0\n", blocksPerGrid, threadsCount);
			simpleSearchKernel<<<blocksPerGrid, threadsCount>>>(clusterArray, funcArray, uint(n), uint(m), destImbArray,
					ulong(h_nc[0]), ulong(numberOfChunks), vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray);
			checkCudaErrors(cudaDeviceSynchronize());
			// printf("Begin reduce...\n");
			thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.begin()+numberOfChunks);
			float min_val = *iter;
			// printf("min element is %.2f\n", min_val);

			if((min_val - bestImbalance < EPS) && (fabs(min_val - bestImbalance) > EPS)) {  // (min_val < bestImbalance) => improvement in imbalance
				// As there may be more than one combination with a better imbalance,
				// chooses one better combination at random.
				// The algorithm starts the search with a random initial index.
				// uint randomInitialIndex = util::RandomUtil::next(0, numberOfChunks - 1);
				// printf("Best imbalance is %.2f and random index is %d\n", bestImbalance, randomInitialIndex);
				// shuffleBestResult1opt<<< 1,1 >>>(destImbArray, numberOfChunks, bestImbalance, randomInitialIndex, resultIndexArray, resultValueArray);
				// thrust::host_vector<uint> h_result_index = d_result_index;
				// thrust::host_vector<float> h_result_value = d_result_value;
				// printf("Result index is %d and result value is %.2f\n", h_result_index[0], h_result_value[0]);

				// determines the position of the best improvement found in the result vector - DISABLED
				uint position = iter - d_destFunctionValue.begin();
				//printf("Original position would be %d and original imbalance would be %.2f\n", position, min_val);
				int resultIdx = position;
				bestImbalance = *iter;
				// bestImbalance = h_result_value[0];
				ulong destCluster = resultIdx / n;
				ulong bestSrcVertex = resultIdx % n;
				ulong sourceCluster = d_mycluster[bestSrcVertex];
				// printf("Idx = %d: The best src vertex is %d to cluster %d (nc = %d) with I(P) = %.2f\n", resultIdx, bestSrcVertex, destCluster, h_nc[0], bestImbalance);
				if(bestImbalance < 0) {  printf("WARNING: I(P) < 0 !!!\n");  }
				updateClustering1opt<<< 1, 1 >>>(bestSrcVertex, destCluster, bestImbalance, clusterArray, funcArray, n, ncArray);
				thrust::host_vector<uint> h_old_nc(1);
				h_old_nc[0] = h_nc[0];
				d_old_nc = h_old_nc;
				h_nc = d_nc;
				sourceVertexList.push_back(bestSrcVertex);
				destinationClusterList.push_back(destCluster);
				// h_functionValue = d_functionValue;
				// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
				// CASO ESPECIAL 1 (um novo cluster k2 foi criado)
				if(h_nc[0] > h_old_nc[0]) {
					// acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
					// printf("New cluster. Growing vectors.\n");
					d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
					d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
					d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
					vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
					vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
					isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
				}

				updateVertexClusterSumArraysDelta<<<1, 1, n*sizeof(int)>>>(weightArray, destArray, numArray, offsetArray, clusterArray, vertexClusterPosSumArray,
						vertexClusterNegSumArray, isNeighborClusterArray, n, old_ncArray, ncArray, bestSrcVertex, sourceCluster, destCluster);
				/*
				blocksPerGrid = (n + threadsCount - 1) / threadsCount;
				d_neighbor_cluster = thrust::device_vector<uint>(n * (h_nc[0]+1), 0);
				uint* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
				isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
				updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
					offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
				*/
				checkCudaErrors(cudaDeviceSynchronize());
				// CASO ESPECIAL 2: o cluster k1 foi removido -> parcialmente tratado dentro do kernel anterior
				// h_mycluster = d_mycluster; // retrieves new cluster configuration from GPU
				// CASO ESPECIAL 2: cluster removido
				if(h_nc[0] < h_old_nc[0]) {
					// remove uma fileira correpondente a um cluster removido na matriz de soma
					// BOOST_LOG_TRIVIAL(debug) << "Deleted cluster. Shrinking vectors." << endl;
					d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
					d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
					d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
					vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
					vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
					isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
				}
				// printf("Preparing new VND loop...\n");
				if(bestImbalance < 0)  break;
			} else {  // no better result found in neighborhood
				// printf("Breaking VND loop...\n");
				break;
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
				ulong(nc), ulong(numberOfChunks), firstImprovement,
				vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount);
		
		checkCudaErrors(cudaDeviceSynchronize());

		h_destcluster1 = d_destCluster1;
		h_destcluster2 = d_destCluster2;
		h_destFunctionValue = d_destFunctionValue;
		
		return true;
	}

	/// Runs the kernel for GRASP with VND.
	bool runGRASPKernel(ClusteringProblem& problem, ConstructClustering &construct, 
			SignedGraph *g, int processRank, ulong timeLimit, int iter,
			thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			ulong n, ulong m, ushort threadsCount, bool firstImprovement,
			Clustering& result, int &totalIterations, double& timeSpentConstruct, double& timeSpentGRASP) {

		// 0. Triggers local processing time calculation
		double timeSpentInGRASP = 0;
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// declaration and initialization of variables
		// Graph data - read only, copy of host data
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		// 1. Construct clustering
		Clustering CStar = construct.constructClustering(g, problem, processRank);
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timer.resume();
		double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);

		Clustering previousCc = CStar;
		Clustering CBest = CStar;
		Clustering Cc = CStar;
		float bestGRASPValue = CStar.getImbalance().getValue();
		unsigned long nc = Cc.getNumberOfClusters();

		// current clustering data - changes every GRASP iteration
		thrust::device_vector<float> d_functionValue(1);
		thrust::device_vector<unsigned long> d_mycluster(1);
		thrust::device_vector<uint> d_randomIndex(n * (nc+1));
		thrust::device_vector<float> d_VertexClusterPosSum(n * (nc+1), 0.0);
		thrust::device_vector<float> d_VertexClusterNegSum(n * (nc+1), 0.0);
		thrust::device_vector<uint> d_nc(1);
		thrust::device_vector<uint> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue(1);
		thrust::device_vector<int> d_neighbor_cluster;
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		uint* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		uint* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
		
		int i = 0, totalIter = 0;
		while(i <= iter || iter < 0) {
			// cout << "CUDA GRASP iteration " << totalIter << endl;
			
			ClusterArray myCluster = Cc.getClusterArray();
			thrust::host_vector<unsigned long> h_mycluster(myCluster);
			thrust::host_vector<float> h_functionValue(1);
			h_functionValue[0] = Cc.getImbalance().getValue();
			// A -> Transfer to device
			// transfers the arrays to CUDA device
			d_mycluster.resize(nc);
			d_mycluster = h_mycluster;
			d_functionValue = h_functionValue;
			clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
			funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
			
			unsigned long numberOfChunks = n * (nc + 1);  // the search space for each vertex (dest cluster) will be split into n*(nc+1) chunks
			// 2. Execute local search algorithm: CUDA VND
			// number of clusters - changes every iteration of VND
			thrust::host_vector<ulong> h_nc(1);
			h_nc[0] = nc;
			d_nc = h_nc;

			// VND loop
			int iteration = 0;
			float bestImbalance = h_functionValue[0];
			int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
			updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
				offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
			checkCudaErrors(cudaDeviceSynchronize());
		
			while (true) {
				// printf("*** Local search iteration %d, nc = %ld, I(P) = %.2f\n", iteration, h_nc[0], bestImbalance);
				numberOfChunks = (h_nc[0] + 1) * n;
				// result / destination vector
				d_destFunctionValue.resize(numberOfChunks);
				float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
			
				// printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
				blocksPerGrid = ((n * (h_nc[0] + 1)) + threadsCount - 1) / threadsCount;
				simpleSearchKernel<<<blocksPerGrid, threadsCount>>>(clusterArray, funcArray,
						uint(n), uint(m), destImbArray, ulong(h_nc[0]), ulong(numberOfChunks),
						vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray);
				checkCudaErrors(cudaDeviceSynchronize());

				// printf("Begin reduce / post-process...\n");
				thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.begin()+numberOfChunks);
				float min_val = *iter;
				if(min_val < bestImbalance) {  // improvement in imbalance
					bestImbalance = min_val;
					// determines the position of the best improvement found in the result vector
					uint position = iter - d_destFunctionValue.begin();
					int resultIdx = position;
					ulong destcluster = resultIdx / n;
					ulong bestSrcVertex = resultIdx % n;
					ulong sourceCluster = d_mycluster[bestSrcVertex];
					// printf("Idx = %d: The best src vertex is %d to cluster %d with I(P) = %.2f\n", resultIdx, bestSrcVertex, destcluster, destFunctionValue);
					if(bestImbalance < 0) {  printf("WARNING: I(P) < 0 !!!\n");  }
					d_old_nc = d_nc;
					old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
					updateClustering1opt<<< 1, 1 >>>(bestSrcVertex, destcluster, bestImbalance, clusterArray, funcArray, n, ncArray);
					int old_nc = h_nc[0];
					h_nc = d_nc;
					nc = h_nc[0];
					// h_functionValue = d_functionValue;
					// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);

					updateVertexClusterSumArraysDelta<<<1, 1>>>(weightArray, destArray, numArray, offsetArray, clusterArray, vertexClusterPosSumArray,
							vertexClusterNegSumArray,isNeighborClusterArray, n, old_ncArray, ncArray, bestSrcVertex, sourceCluster, destcluster);
					checkCudaErrors(cudaDeviceSynchronize());
					// CASO ESPECIAL 1: o cluster k1 foi removido -> tratado dentro do kernel anterior
					/*
					if(nc < old_nc) {
						// verifica se a fileira inteira correspondente ao cluster removido ficou zerada: OK
						bool zerado = true;
						for(int v = 0; v < n; v++) {
							if( (vertexClusterPosSum[v + k1 * n] > 0) or (vertexClusterNegSum[v + k1 * n] > 0) ) {
								zerado = false;
								break;
							}
						}
						// BOOST_LOG_TRIVIAL(debug) << "*** Cluster removido. Zerado = " << zerado << endl;
						// remove a fileira correspondente ao cluster k1 na matriz de soma, shiftando os dados
						for(int k = k1 + 1; k <= nc; k++) {
							for(int v = 0; v < n; v++) {
								vertexClusterPosSum[(k - 1) * n + v] = vertexClusterPosSum[(k) * n + v];
								vertexClusterNegSum[(k - 1) * n + v] = vertexClusterNegSum[(k) * n + v];
							}
						}
						// remove fisicamente a ultima fileira da matriz - no proximo if
					} */
					h_mycluster = d_mycluster; // retrieves new cluster configuration from GPU
					// CASO ESPECIAL 2 (um novo cluster k2 foi criado) E CASO ESPECIAL 1
					if(nc != old_nc) {
						// acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
						// BOOST_LOG_TRIVIAL(debug) << "New cluster. Growing vectors." << endl;
						d_VertexClusterPosSum.resize(n * (nc+1), 0.0);
						d_VertexClusterNegSum.resize(n * (nc+1), 0.0);
						vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
						vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
						// BOOST_LOG_TRIVIAL(debug) << "New element value: " << vertexClusterPosSum[new_nc * n + 2] << endl;
					}
					// printf("Preparing new VND loop...\n");
					if(bestImbalance < 0)  break;
				} else {  // no better result found in neighborhood
					// printf("Breaking VND loop...\n");
					break;
				}
				iteration++;
			}
			// h_mycluster = d_mycluster;
			nc = h_nc[0];

			// 3. Select the best clustring so far
			// if Q(Cl) > Q(Cstar)
			// Imbalance newValue = Cl.getImbalance();
			if(bestImbalance < bestGRASPValue) {
				ClusterArray cArray;
				for(int x = 0; x < h_mycluster.size(); x++) {
					cArray.push_back(h_mycluster[x]);
				}
				bestGRASPValue = bestImbalance;
				// printf("Rebuilding clustering object.\n");
				/*  TODO FIXME consertar esse bug aqui
				Clustering Cl(cArray, *g, problem);
				printf("I(P) = %.2f\n", Cl.getImbalance().getValue());
				CStar = Cl;
				printf("Clustering object REBUILT.\n"); */
				// iterationValue = totalIter;
				i = 0;
				if(bestGRASPValue <= 0) break;
			}
			timer.stop();
			boost::timer::cpu_times end_time = timer.elapsed();
			timer.resume();
			timeSpentInGRASP = (end_time.wall - start_time.wall) / double(1000000000);
			// if elapsed time is bigger than timeLimit, break
			if(timeSpentInGRASP >= timeLimit) {
				cout << "Time limit exceeded." << endl;
				break;
			}

			// Increment to next loop
			i++, totalIter++, previousCc = Cc;

			// guarantees at least one execution of the GRASP when the number of iterations is smaller than one
			// used in vote/boem tests
			if(iter <= 0) {break;}

			// Avoids constructClustering if loop break condition is met
			if(i <= iter) {
				timer.stop();
				boost::timer::cpu_times start_timeC = timer.elapsed();
				timer.resume();

				// 1. Construct the next clustering
				// printf("New construct clustering\n");
				Cc = construct.constructClustering(g, problem, processRank);
				int old_nc = nc;
				nc = Cc.getNumberOfClusters();
				if(old_nc != nc) {  // the number of clusters in the solution changed
					d_VertexClusterPosSum.resize(n * (nc+1));
					d_VertexClusterNegSum.resize(n * (nc+1));
					vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
					vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
					// printf("Number of clusters has changed.\n");
				}

				timer.stop();
				end_time = timer.elapsed();
				timer.resume();
				timeSpentInConstruction += (end_time.wall - start_timeC.wall) / double(1000000000);
			}
			// printf("I(P) = %.2f\n", bestGRASPValue);
		}
		printf("I(P) = %.2f\n", bestGRASPValue);
		result = CStar;
		totalIterations = totalIter;
		timeSpentConstruct = timeSpentInConstruction;
		timeSpentGRASP = timeSpentInGRASP;
		return true;
	}

	/// Runs the kernel for metaheuristic construct clustering's gain function calculation.
	bool runConstructKernel(ulong randomSeed, thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
			thrust::host_vector<unsigned long>& h_mycluster, thrust::host_vector<float>& h_functionValue, 
			ulong n, ulong m, ulong nc, ushort threadsCount, thrust::host_vector<unsigned long>& h_newcluster, 
			double& imbalance) {

		// declaration and initialization of variables
		// Graph data - read only, copy of host data
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<int> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<int> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<int> d_offset = h_offset;  // initial edge number for vertex i
		// 1. Current (partial) clustering data
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<float> d_VertexClusterPosSum(n * (nc+1), 0.0);
		thrust::device_vector<float> d_VertexClusterNegSum(n * (nc+1), 0.0);
		thrust::device_vector<uint> d_nc(1);
		thrust::device_vector<uint> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue(1);
		thrust::device_vector<int> d_neighbor_cluster;
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		int* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		int* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		int* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		uint* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		uint* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] ); // TODO calcular valor aqui
		
		VertexSet lc(randomSeed, n); // L(Cc) = V(G)
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = nc;
		
		int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
		updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(long)>>>(weightArray, destArray, numArray,
			offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
		checkCudaErrors(cudaDeviceSynchronize());
		
		while (lc.size() > 0) { // lc != empty
			// alpha = 1.0 (completely random): no need to calculate all gains (saves time)
			int v = lc.chooseRandomVertex(lc.size()).vertex;
			
			// CALCULO DO GAIN PARA O VERTICE i
			d_nc = h_nc;
			
			// result / destination vector
			// cout << "The current number of clusters is " << h_nc[0] << " and there are " << lc.size() << " remaining vertices." << endl;
			ulong numberOfChunks = (h_nc[0] + 1);	
			d_destFunctionValue.resize(numberOfChunks);
			float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );		
			blocksPerGrid = ((h_nc[0] + 1) + threadsCount - 1) / threadsCount;
			simpleConstructKernel<<< blocksPerGrid, threadsCount >>>(clusterArray, funcArray, uint(n), uint(m), destImbArray,
				ulong(h_nc[0]), vertexClusterPosSumArray, vertexClusterNegSumArray, threadsCount, v);
			checkCudaErrors(cudaDeviceSynchronize());
	
			/*
			thrust::host_vector<float> h_destFunctionValue = d_destFunctionValue;
			cout << "destFunctionValue array: ";
			for(int x = 0; x < h_destFunctionValue.size(); x++) {
				cout << h_destFunctionValue[x] << " ";
			}
			cout << endl; */

			// printf("Begin reduce / post-process...\n");
			thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.begin()+numberOfChunks);
			float min_val = *iter;
			double bestImbalance = min_val;
			// determines the position of the best improvement found in the result vector
			uint position = iter - d_destFunctionValue.begin();
			int resultIdx = position;
			int clusterNumber = resultIdx;
			// invoca o kernel para atualizar o vetor de clustering (clusterArray) e o vetor de F.O. (funcArray)
			// printf("Inserting vertex %d into cluster %d with imbalance = %.2f\n", v, clusterNumber, bestImbalance);
			updateConstructClustering <<< 1,1 >>>(v, clusterNumber, bestImbalance, clusterArray, funcArray, n, h_nc[0]);
			int prevCluster = h_nc[0];  // before being moved to a cluster, vertex v is initially inside the cluster where k = nc.
			d_old_nc = h_nc;
			
			if(clusterNumber == h_nc[0]) {  // um novo cluster k2 sera criado
				// clusterNumber = Clustering::NEW_CLUSTER;
				h_nc[0]++;  // just increment the number of clusters
				d_nc = h_nc;
				// acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
				// printf("New cluster. Growing vectors.\n");
				d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
				d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
				vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
				vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
			}
			//printf("Idx = %d: The best src vertex is %ld to cluster %ld with I(P) = %.2f\n", resultIdx, v, clusterNumber, bestImbalance);
			if(bestImbalance < 0) {  printf("WARNING: Construct I(P) < 0 !!! I(P) = %.2f\n", bestImbalance);  }

			updateVertexClusterSumArraysDelta <<< 1,1 >>>(weightArray, destArray, numArray,
				offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, old_ncArray, ncArray,
				v, h_nc[0], clusterNumber);
			checkCudaErrors(cudaDeviceSynchronize());
			// 4. lc = lc - {i}
			// the choosing vertex i automatically removes it from the list
			// Removal already done by the chooseVertex methods above
		}
		h_newcluster = d_mycluster;
		h_functionValue = d_functionValue;
		imbalance = h_functionValue[0];
		
		return true;
	}

}

