#include "include/CUDASearch.h"
#include "include/CUDAHelper.h"

#include "problem/include/ClusteringProblem.h"
#include "../construction/include/ConstructClustering.h"
#include "graph/include/Graph.h"
#include "util/include/RandomUtil.h"
#include "../construction/include/VertexSet.h"
#include "graph/include/Perturbation.h"
#include "./include/VariableNeighborhoodDescent.h"

#include <boost/timer/timer.hpp>
#include <iostream>
#include <iomanip>

#include <assert.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>
#include <vector>
#include <utility>

// define precision for floating point single-precision in CUDA
#define EPS 0.00005

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
#define BLOCK_DIM_X 16
#define BLOCK_DIM_Y 16
#define BLOCK_SIZE 512


#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)

namespace clusteringgraph {

using namespace problem;
using namespace resolution::construction;
using namespace resolution::vnd;
using namespace clusteringgraph;
using namespace std;
	
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
	__global__ void updateClustering1opt(const int bestSrcVertex, const int destcluster, float positiveFunctionValue,
			float negativeFunctionValue, ulong* clusterArray, float* funcArray, long n, ulong* nc) {
		funcArray[0] = positiveFunctionValue;
		funcArray[1] = negativeFunctionValue;
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
	__global__ void updateVertexClusterSumArrays(const float* weightArray, const long* destArray, const long* numArray,
			const long* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int* isNeighborClusterArray, long n, ulong* ncArray) {
		
		long i = blockDim.x*blockIdx.x + threadIdx.x;
		ulong nc = ncArray[0];
		// shared memory
		/*
		extern __shared__ int sc[];       // n longs
		int *s_cluster = sc;              // n longs
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads(); */
	    if(i < n) {
            // For each vertex i, stores the sum of edge weights between vertex i and all clusters
            for(int k = 0; k <= nc; k++) {
            	vertexClusterPosSumArray[k * n + i] = 0.0;
            	vertexClusterNegSumArray[k * n + i] = 0.0;
            	isNeighborClusterArray[k * n + i] = 0;
            }
            // in/out-edges of vertex i
            ulong offset = offsetArray[i];
			ulong numedges = numArray[i];
			ulong count = offset + numedges;
			long ki = clusterArray[i];
			for (ulong edgenum = offset; edgenum < count; edgenum++) {   
				long j = destArray[edgenum];
				float weight = weightArray[edgenum];
				long kj = clusterArray[j];
				if(kj >= 0) {
					if(ki != kj) {  // different cluster
						isNeighborClusterArray[kj * n + i]++;  // vertex i now has external connection to cluster kj
						isNeighborClusterArray[ki * n + j]++;  // vertex j now has external connection to cluster ki
					}
					if(weight > 0) {
						vertexClusterPosSumArray[kj * n + i] += fabs(weight);
					} else {
						vertexClusterNegSumArray[kj * n + i] += fabs(weight);
					}
				}
			}
			// preenche a possibilidade de se mover o vertice i para um cluster novo (k = nc)
			isNeighborClusterArray[nc * n + i] = 1;
        }
	}
	
	/// CUDA Kernel to update the vertex-cluster edge-weight sum arrays after a change in the clustering.
	__global__ void updateVertexClusterSumArraysDelta(const float* weightArray, const long* destArray, const long* numArray,
			const long* offsetArray, const ulong* clusterArray, float* vertexClusterPosSumArray, float* vertexClusterNegSumArray, 
			int* isNeighborClusterArray, long n, ulong* old_ncArray, ulong* ncArray, int i, int k1, int k2) {  // vertex i is being moved from cluster k1 to k2
		ulong nc = old_ncArray[0];
		ulong new_nc = ncArray[0];
		// shared memory
		/* extern __shared__ int sc[];       // n longs
		int *s_cluster = sc;              // n longs
		for(ulong j = 0; j < n; j++) {
			s_cluster[j] = clusterArray[j];
		}
		__syncthreads(); */

		/*  THIS PEACE OF CODE IS NEEEDED ONLY WHEN USING THE CUDA CONSTRUCTIVE KERNEL - DISABLED */
		if(new_nc > nc) {  // vertex i is being moved to a new cluster
			// move a fileira correspondente ao cluster k = nc na matriz de soma, shiftando os dados para a direita (nc + 1)
			for(long v = 0; v < n; v++) {
				vertexClusterPosSumArray[(new_nc) * n + v] = vertexClusterPosSumArray[(nc) * n + v];
				vertexClusterNegSumArray[(new_nc) * n + v] = vertexClusterNegSumArray[(nc) * n + v];
			}
			// zera a fileira movida anteriormente
			for(long v = 0; v < n; v++) {
				vertexClusterPosSumArray[(nc) * n + v] = 0.0;
				vertexClusterNegSumArray[(nc) * n + v] = 0.0;
			}
			// preenche a possibilidade de se mover todos os vertices para um cluster novo (k = new_nc)
			for(long v = 0; v < n; v++) {
				isNeighborClusterArray[v + new_nc * n] = 1;
				isNeighborClusterArray[v + nc * n] = 0;
			}
		}

		// in/out-edges of vertex i
		ulong offset = offsetArray[i];
		ulong numedges = numArray[i];
		ulong count = offset + numedges;
		// isNeighborClusterArray[i+k1*n] = 1;  // vertex i now has external connection to cluster k1
		/* for(ulong k = 0; k < nc; k++) {
		 	isNeighborClusterArray[i + k * n] = 0;
		} */
		for (ulong edgenum = offset; edgenum < count; edgenum++) {
			ulong j = destArray[edgenum];
			float weight = weightArray[edgenum];
			long kj = clusterArray[j];

			isNeighborClusterArray[j + k1 * n] -= 2;
			if(kj != k2) {
				isNeighborClusterArray[i + kj * n]++;  // vertex i now has external connection to cluster kj
				isNeighborClusterArray[j + k2 * n]++;  // vertex j now has external connection to cluster k2
			} else {  // cluster[i] == cluster[j] == k2
				isNeighborClusterArray[i + kj * n] = 0;  // vertex i now has NO external connections to cluster kj
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
		/*
		if(new_nc < nc) {
			// remove a fileira correspondente ao cluster k1 na matriz de soma, shiftando os dados para a esquerda
			for(int k = k1 + 1; k <= nc; k++) {
				for(long v = 0; v < n; v++) {
					vertexClusterPosSumArray[(k - 1) * n + v] = vertexClusterPosSumArray[(k) * n + v];
					vertexClusterNegSumArray[(k - 1) * n + v] = vertexClusterNegSumArray[(k) * n + v];
					isNeighborClusterArray[(k - 1) * n + v] = isNeighborClusterArray[(k) * n + v];
				}
			}
		} */
	}

	/// CUDA kernel for simple 1-opt search.
	__global__ void simpleSearchKernel(const ulong* clusterArray, const float* funcArray, ulong n, ulong m,
		float* destImbArray, float* destPosImbArray, float* destNegImbArray, ulong nc, ulong numberOfChunks,
		float* positiveSumArray, float* negativeSumArray, int* isNeighborClusterArray) {

		// ulong idx = blockIdx.x * blockDim.x + threadIdx.x;
		int index_x = blockIdx.x * blockDim.x + threadIdx.x;
		int index_y = blockIdx.y * blockDim.y + threadIdx.y;

		// map the two 	2D indices to a single linear, 1D index
		int grid_width = gridDim.x * blockDim.x;
		int idx = index_y * grid_width + index_x;

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
				destImbArray[idx] = funcArray[0] + funcArray[1] + positiveSum + negativeSum;
				//destPosImbArray[idx] = funcArray[0] + positiveSum;
				//destNegImbArray[idx] = funcArray[1] + negativeSum;
			} else {
				destImbArray[idx] = funcArray[0] + funcArray[1];
				//destPosImbArray[idx] = funcArray[0];
                                //destNegImbArray[idx] = funcArray[1];
			}
		} else if(idx < numberOfChunks) {
			destImbArray[idx] = funcArray[0] + funcArray[1];
                        //destPosImbArray[idx] = funcArray[0];
                        //destNegImbArray[idx] = funcArray[1];
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
						int* resultIndexArray) {
		// Each thread will check if the result corresponding to its index is an improvement
		uint idx = blockIdx.x * blockDim.x + threadIdx.x;
		if(idx < imbArraySize) {
			if((destImbArray[idx] - bestImbalance < EPS) && (fabs(destImbArray[idx] - bestImbalance) > EPS)){   // (destImbArray[i] < bestImbalance)
				resultIndexArray[0] = idx;
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
	
	/// CUDA Kernel to update the cluster array and the objective function value after a permutation.
	__global__ void perturbationKernel(int node, int k2,  
			ulong* clusterArray, float* funcArray, int n, uint* nc, float* positiveSumArray,
			float* negativeSumArray) {
		
		if(k2 == clusterArray[node]) {  // Option 1 (easier): destcluster == cluster[node]
			k2 = nc[0];
		}
		int k1 = clusterArray[node];
		// updates objective function value
		float positiveSum = +(- positiveSumArray[node+k2*n]) - (- positiveSumArray[node+k1*n]);
		float negativeSum = -negativeSumArray[node+k1*n] + negativeSumArray[node+k2*n];
		// updates thread idx / vertex i from k1 to k2 imbalance result
		funcArray[0] += positiveSum;
		funcArray[1] += negativeSum;

		clusterArray[node] = k2;
		// validate the current number of clusters in the solution
		if(k2 >= nc[0]) {  // Option 1 (easier): destcluster == number of clusters
			nc[0]++;  // i.e., bestSrcVertex is going to a standalone cluster
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
		if(not found) {  // previousCluster is empty, has to be removed
			for(int i = 0; i < n; i++) {
				if(clusterArray[i] > k1) {
					clusterArray[i]--;
				}
			}
			nc[0]--;
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
		float& destPositiveImbalance, float& destNegativeImbalance, const long& timeSpentSoFar, const unsigned int& l) {

		// declaration and initialization of variables
		// Graph data - read only
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<long> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<long> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<long> d_offset = h_offset;  // initial edge number for vertex i
		// current clustering data
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<uint> d_randomIndex(n * (nc+1));
		thrust::device_vector<float> d_VertexClusterPosSum(n * (nc+1));
		thrust::device_vector<float> d_VertexClusterNegSum(n * (nc+1));
		thrust::device_vector<ulong> d_nc(1);
		thrust::device_vector<ulong> d_old_nc(1);
		thrust::device_vector<ulong> d_result_index(1);
		thrust::device_vector<float> d_destFunctionValue(numberOfChunks);
		thrust::device_vector<float> d_destPosFunctionValue(numberOfChunks);
		thrust::device_vector<float> d_destNegFunctionValue(numberOfChunks);
		thrust::device_vector<unsigned long> d_destCluster1(numberOfChunks);
		thrust::device_vector<unsigned long> d_destCluster2(numberOfChunks);
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		long* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		long* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		long* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		uint* randomIndexArray = thrust::raw_pointer_cast( &d_randomIndex[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		ulong* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		ulong* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		ulong* resultIndexArray = thrust::raw_pointer_cast( &d_result_index[0] );
		
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
		thrust::device_vector<float> d_destPosFunctionValue(numberOfChunks);
		thrust::device_vector<float> d_destNegFunctionValue(numberOfChunks);
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
			Clustering& result, int &totalIterations, double& timeSpentConstruct, double& timeSpentGRASP,
			stringstream &constructivePhaseResults, stringstream &iterationResults) {

		// 0. Triggers local processing time calculation
		double timeSpentInGRASP = 0;
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// declaration and initialization of variables
		// Graph data - read only, copy of host data
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<long> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<long> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<long> d_offset = h_offset;  // initial edge number for vertex i
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
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = Cc.getNumberOfClusters();

		// current clustering data - changes every GRASP iteration
		thrust::device_vector<float> d_functionValue(2);
		thrust::device_vector<unsigned long> d_mycluster(1);
		thrust::device_vector<uint> d_randomIndex(n * (h_nc[0]+1));
		thrust::device_vector<float> d_VertexClusterPosSum(n * (h_nc[0]+1), 0.0);
		thrust::device_vector<float> d_VertexClusterNegSum(n * (h_nc[0]+1), 0.0);
		thrust::device_vector<ulong> d_nc(1);
		thrust::device_vector<ulong> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue(1);
		thrust::device_vector<float> d_destPosFunctionValue(1);
		thrust::device_vector<float> d_destNegFunctionValue(1);
		thrust::device_vector<int> d_neighbor_cluster(n * (h_nc[0]+1), 0);
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		long* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		long* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		long* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		ulong* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		ulong* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
		
		int i = 0, totalIter = 0;
		while(i <= iter || iter < 0) {
			// cout << "CUDA GRASP iteration " << totalIter << endl;
			
			ClusterArray myCluster = Cc.getClusterArray();
			thrust::host_vector<unsigned long> h_mycluster(myCluster);
			thrust::host_vector<float> h_functionValue(2);
			float destPositiveImbalance = Cc.getImbalance().getPositiveValue();
			float destNegativeImbalance = Cc.getImbalance().getNegativeValue();
			h_nc[0] = Cc.getNumberOfClusters();
			// A -> Transfer to device
			// transfers the arrays to CUDA device
			d_mycluster.resize(h_nc[0]);
			d_mycluster = h_mycluster;
			d_functionValue = h_functionValue;
			clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
			funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
			
			unsigned long numberOfChunks = n * (h_nc[0] + 1);  // the search space for each vertex (dest cluster) will be split into n*(nc+1) chunks
			// 2. Execute local search algorithm: CUDA VND
			// number of clusters - changes every iteration of VND
			d_nc = h_nc;
			ncArray = thrust::raw_pointer_cast( &d_nc[0] );

			// VND loop
			int iteration = 0;
			float bestImbalance = destPositiveImbalance + destNegativeImbalance;
			long blocksPerGrid = (n + threadsCount - 1) / threadsCount;
			updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(int)>>>(weightArray, destArray, numArray,
				offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
			checkCudaErrors(cudaDeviceSynchronize());
		
			while (true) {
				// printf("*** Local search iteration %d, nc = %ld, I(P) = %.2f\n", iteration, h_nc[0], bestImbalance);
				numberOfChunks = (h_nc[0] + 1) * n;
				// result / destination vector
				d_destFunctionValue.resize(numberOfChunks);
				//d_destPosFunctionValue.resize(numberOfChunks);
				//d_destNegFunctionValue.resize(numberOfChunks);
				float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
				float* destPosImbArray = thrust::raw_pointer_cast( &d_destPosFunctionValue[0] );
				float* destNegImbArray = thrust::raw_pointer_cast( &d_destNegFunctionValue[0] );
				d_functionValue[0] = bestImbalance;
				d_functionValue[1] = 0.0;
			
				// printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
				                                        // create two dimensional 4x4 thread blocks
                                // dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
                                dim3 block(BLOCK_SIZE, 1, 1);

                                // configure a two dimensional grid as well
                                dim3 grid_size;
                                long num_elements_x = (numberOfChunks + (block.x) - 1) / block.x;
                                long num_elements_y = 1;
                                if(num_elements_x > 65535) {
                                        num_elements_y = (num_elements_x + 65535 - 1) / 65535;
                                        num_elements_x = 65535;
                                }
                                grid_size.x = num_elements_x;// / block.x;
                                grid_size.y = num_elements_y;// / block.y;

				simpleSearchKernel<<<grid_size, block>>>(clusterArray, funcArray,
						ulong(n), ulong(m), destImbArray, destPosImbArray, destNegImbArray, ulong(h_nc[0]), ulong(numberOfChunks),
						vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray);
				checkCudaErrors(cudaDeviceSynchronize());

				// printf("Begin reduce / post-process...\n");
				thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(),
						d_destFunctionValue.begin()+numberOfChunks);
				float min_val = *iter;
				// as there may be more than one combination with the same imbalance,
				// chooses one combination at random
				// the algorithm starts the search with a random initial index
				/*
				thrust::host_vector<int> h_result_index(1);
				h_result_index[0] = -1;
				thrust::device_vector<int> d_result_index = h_result_index;
				int* resultIndexArray = thrust::raw_pointer_cast( &d_result_index[0] );
				blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
				shuffleBestResult1opt<<< blocksPerGrid, threadsCount >>>(destImbArray, numberOfChunks, bestImbalance, resultIndexArray);
				h_result_index = d_result_index;
				int resultIdx = h_result_index[0];
				float min_val = bestImbalance;
				if(resultIdx >= 0) {
					min_val = d_destFunctionValue[resultIdx];
				} */

				if((min_val - bestImbalance < EPS) && (fabs(min_val - bestImbalance) > EPS)) {  // (min_val < bestImbalance) => improvement in imbalance
					bestImbalance = min_val;
					// determines the position of the best improvement found in the result vector
					ulong position = iter - d_destFunctionValue.begin();
					long resultIdx = position;
					ulong destCluster = resultIdx / n;
					ulong bestSrcVertex = resultIdx % n;
					ulong sourceCluster = d_mycluster[bestSrcVertex];
					//destPositiveImbalance = d_destPosFunctionValue[resultIdx];
					//destNegativeImbalance = d_destNegFunctionValue[resultIdx];
					destPositiveImbalance = bestImbalance;
                                        destNegativeImbalance = 0.0;

					// printf("Idx = %d: The best src vertex is %d to cluster %d with I(P) = %.2f\n", resultIdx, bestSrcVertex, destCluster, bestImbalance);
					if(bestImbalance < EPS) {  printf("WARNING: I(P) < 0 !!!\n");  }
					d_old_nc = d_nc;
					old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
					updateClustering1opt<<< 1, 1 >>>(bestSrcVertex, destCluster, destPositiveImbalance, destNegativeImbalance, clusterArray, funcArray, n, ncArray);
					checkCudaErrors(cudaDeviceSynchronize());
					thrust::host_vector<ulong> h_old_nc(1);
					h_old_nc = h_nc;
					h_nc = d_nc;
					// h_functionValue = d_functionValue;
					// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);
					// TODO check if full recalculation of the arrays is faster than delta calculation
					/*
                                        if(h_nc[0] > h_old_nc[0]) {
                                                // acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
                                                // printf("New cluster. Growing vectors.\n");
                                                d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
                                                d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
                                                d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
                                                vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
                                                vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
                                                isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
                                        } */
                                        if(h_nc[0] == h_old_nc[0]) {  // must not update the arrays if the number of clusters decreased
                                                updateVertexClusterSumArraysDelta<<<1, 1, n*sizeof(int)>>>(weightArray, destArray, numArray, offsetArray, clusterArray, vertexClusterPosSumArray,
                                                        vertexClusterNegSumArray, isNeighborClusterArray, n, old_ncArray, ncArray, bestSrcVertex, sourceCluster, destCluster);
                                                checkCudaErrors(cudaDeviceSynchronize());
                                        }
                                        // CASO ESPECIAL 2: o cluster k1 foi removido -> parcialmente tratado dentro do kernel anterior
                                        // h_mycluster = d_mycluster; // retrieves new cluster configuration from GPU
                                        // CASO ESPECIAL 2: cluster removido
                                        if(h_nc[0] != h_old_nc[0]) {  // number of cluster decreased => FULL recalculation
                                                        // recalculates sum matrices
                                                        d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
                                                        d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
                                                        d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
                                                        vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
                                                        vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
                                                        isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
                                                        long blocksPerGrid = (n + threadsCount - 1) / threadsCount;  // , n*sizeof(int)
                                                        updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount>>>(weightArray, destArray, numArray,
                                                                offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
                                                        checkCudaErrors(cudaDeviceSynchronize());
                                        }
					// printf("Preparing new VND loop...\n");
					if(bestImbalance <= 0)  break;
				} else {  // no better result found in neighborhood
					// printf("Breaking VND loop...\n");
					break;
				}
				iteration++;
			}
			// 3. Select the best clustring so far
			// if Q(Cl) > Q(Cstar)
			// Imbalance newValue = Cl.getImbalance();
			if((bestImbalance - bestGRASPValue < EPS) && (fabs(bestImbalance - bestGRASPValue) > EPS)) {  // (bestImbalance < bestGRASPValue) => improvement in imbalance
				// printf("Imbalance improved.\n");
				h_mycluster = d_mycluster;
				ClusterArray cArray(h_mycluster.size(), 0);
				for(int x = 0; x < h_mycluster.size(); x++) {
					cArray[x] = h_mycluster[x];
				}
				bestGRASPValue = bestImbalance;
				// printf("Rebuilding clustering object.\n");
				Clustering Cl(cArray, *g, problem);
				// printf("I(P) = %.2f\n", Cl.getImbalance().getValue());
				CStar = Cl;
				// printf("Clustering object REBUILT.\n");
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
				int old_nc = h_nc[0];
				h_nc[0] = Cc.getNumberOfClusters();
				if(old_nc != h_nc[0]) {  // the number of clusters in the solution changed
					d_VertexClusterPosSum.resize(n * (h_nc[0]+1));
					d_VertexClusterNegSum.resize(n * (h_nc[0]+1));
					d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
					vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
					vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
					isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
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
		thrust::device_vector<long> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<long> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<long> d_offset = h_offset;  // initial edge number for vertex i
		// 1. Current (partial) clustering data
		thrust::device_vector<float> d_functionValue = h_functionValue;
		thrust::device_vector<unsigned long> d_mycluster = h_mycluster;
		thrust::device_vector<float> d_VertexClusterPosSum(n * (nc+1), 0.0);
		thrust::device_vector<float> d_VertexClusterNegSum(n * (nc+1), 0.0);
		thrust::device_vector<ulong> d_nc(1);
		thrust::device_vector<ulong> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue(1);
		thrust::device_vector<int> d_neighbor_cluster;
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		long* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		long* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		long* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		ulong* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		ulong* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] ); // TODO calcular valor aqui
		
		VertexSet lc(randomSeed, n); // L(Cc) = V(G)
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = nc;
		
		int blocksPerGrid = (n + threadsCount - 1) / threadsCount;
		updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount, n*sizeof(int)>>>(weightArray, destArray, numArray,
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
			ulong position = iter - d_destFunctionValue.begin();
			long resultIdx = position;
			long clusterNumber = resultIdx;
			// invoca o kernel para atualizar o vetor de clustering (clusterArray) e o vetor de F.O. (funcArray)
			// printf("Inserting vertex %d into cluster %d with imbalance = %.2f\n", v, clusterNumber, bestImbalance);
			updateConstructClustering <<< 1,1 >>>(v, clusterNumber, bestImbalance, clusterArray, funcArray, n, h_nc[0]);
			long prevCluster = h_nc[0];  // before being moved to a cluster, vertex v is initially inside the cluster where k = nc.
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
			if(bestImbalance < -EPS) {  printf("WARNING: Construct I(P) < 0 !!! I(P) = %.2f\n", bestImbalance);  }

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

	 int generateRandomNumbers(size_t n, int lowerBound, int upperBound, std::vector<int>& randomNumbers) {
	    size_t i;
	    curandGenerator_t gen;
	    float *devData, *hostData;

	    /* Allocate n floats on host */
	    hostData = (float *)calloc(n, sizeof(float));

	    /* Allocate n floats on device */
	    CUDA_CALL(cudaMalloc((void **)&devData, n*sizeof(float)));

	    /* Create pseudo-random number generator */
	    CURAND_CALL(curandCreateGenerator(&gen,
	                CURAND_RNG_PSEUDO_DEFAULT));

	    /* Set seed */
	    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen,
	    		boost::random::random_device()()));

	    /* Generate n floats on device */
	    CURAND_CALL(curandGenerateUniform(gen, devData, n));

	    /* Copy device memory to host */
	    CUDA_CALL(cudaMemcpy(hostData, devData, n * sizeof(float),
	        cudaMemcpyDeviceToHost));

	    /* Show result */
	    randomNumbers.resize(n, 0);
	    for(i = 0; i < n; i++) {
	        int rnd_integer_range = lowerBound + hostData[i] * (upperBound - lowerBound);
	        // printf("%1.4f  %d  ", hostData[i], rnd_integer_range);
	        randomNumbers[i] = rnd_integer_range;
	    }
	    // printf("\n");

	    /* Cleanup */
	    CURAND_CALL(curandDestroyGenerator(gen));
	    CUDA_CALL(cudaFree(devData));
	    free(hostData);

	    return EXIT_SUCCESS;
	}

	bool runILSKernel(ClusteringProblem& problem, ConstructClustering &construct,
					SignedGraph *g, int processRank, ulong timeLimit,
					const int& iterMax, const int& iterMaxILS, const int& perturbationLevelMax,
					thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
					thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset,
					ulong n, ulong m, ushort threadsCount, bool firstImprovement,
					Clustering& result, int &totalIterations, double& timeSpentConstruct, double& timeSpentILS,
					stringstream &constructivePhaseResults, stringstream &iterationResults) {		

		double timeSpentInILS = 0;
		Perturbation perturbation(construct.getRandomSeed());
		util::RandomUtil randomUtil;
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();
		// declaration and initialization of variables
		// Graph data - read only, copy of host data
		thrust::device_vector<float> d_weights = h_weights;  // edge weights
		thrust::device_vector<long> d_dest = h_dest;  // edge destination (vertex j)
		thrust::device_vector<long> d_numedges = h_numedges;  // number of edges of each vertex i
		thrust::device_vector<long> d_offset = h_offset;  // initial edge number for vertex i
		// random number engine
		thrust::host_vector<float> h_random;
		thrust::host_vector<float> h_random2;
		thrust::device_vector<float> d_random;
		curandGenerator_t gen;
		curandGenerator_t gen2;
		/* Create pseudo-random number generator */
		curandCreateGenerator(&gen,	CURAND_RNG_PSEUDO_DEFAULT);
		curandCreateGenerator(&gen2,	CURAND_RNG_PSEUDO_DEFAULT);
		/* Set seed */
		curandSetPseudoRandomGeneratorSeed(gen,	boost::random::random_device()());
		curandSetPseudoRandomGeneratorSeed(gen2, boost::random::random_device()());

		// 1. Initial clustering (construct)
		Clustering Cc = construct.constructClustering(g, problem, processRank);
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
		timeSpentInILS += timeSpentInConstruction;
		timer.resume();
		start_time = timer.elapsed();
		
		Clustering previousCc = Cc;
		Clustering CBest = Cc;
		Clustering CStar = Cc;
		int iterationValue = 0;
		double timeSpentOnBestSolution = 0.0;
		double initialImbalanceSum = 0.0;
		
		thrust::host_vector<ulong> h_nc(1);
		h_nc[0] = CStar.getNumberOfClusters();

		// current clustering data - changes every ILS iteration
		thrust::device_vector<float> d_functionValue(2);
		thrust::device_vector<unsigned long> d_mycluster(1);
		thrust::device_vector<uint> d_randomIndex(n * (h_nc[0]+1));
		thrust::device_vector<float> d_VertexClusterPosSum(n * (h_nc[0]+1), 0.0);
		thrust::device_vector<float> d_VertexClusterNegSum(n * (h_nc[0]+1), 0.0);
		thrust::device_vector<ulong> d_nc(1);
		thrust::device_vector<ulong> d_old_nc(1);
		thrust::device_vector<float> d_destFunctionValue(1);
		thrust::device_vector<float> d_destPosFunctionValue(1);
		thrust::device_vector<float> d_destNegFunctionValue(1);
		thrust::device_vector<int> d_neighbor_cluster(n * (h_nc[0]+1), 0);
		
		float* weightArray = thrust::raw_pointer_cast( &d_weights[0] );
		long* destArray = thrust::raw_pointer_cast( &d_dest[0] );
		long* numArray = thrust::raw_pointer_cast( &d_numedges[0] );
		long* offsetArray = thrust::raw_pointer_cast( &d_offset[0] );
		
		unsigned long* clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
		float* funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
		float* vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
		float* vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
		ulong* ncArray = thrust::raw_pointer_cast( &d_nc[0] );
		ulong* old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
		int* isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
		
		// Multi-start ILS
		for (int i = 0; i < iterMax || iterMax < 0 ; i++, previousCc = Cc) {
			double totalTime = 0.0;
			 printf("ILS iteration %d, best solution so far: %.2f\n", i, CBest.getImbalance().getValue());
			// cout << "Best solution so far: I(P) = " << fixed << setprecision(0) << bestValue.getValue() << endl;
			//    Store initial solution value in corresponding results file
			constructivePhaseResults << (i+1) << "," << Cc.getImbalance().getValue() << ","
					<< Cc.getImbalance().getPositiveValue()
					<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
			initialImbalanceSum += Cc.getImbalance().getValue();
			// Registers the Cc result
			// notifyNewValue(Cc, 0.0, i);

			CStar = Cc;
			Clustering Cl = Cc;
			float bestImbalance = CStar.getImbalance().getValue();
			float destPositiveImbalance = Cc.getImbalance().getPositiveValue();
			float destNegativeImbalance = Cc.getImbalance().getNegativeValue();
			// printf("Constructive phase I(P) = %.2f\n", bestImbalance);
			int perturbationLevel = 1;
			bool perturbated = false;
			
			for(int j = 1, total = 0; j <= iterMaxILS; total++) {  // internal ILS loop
				// printf("ILS internal loop iteration %d\n", j);
				ClusterArray cArrayCl = Cl.getClusterArray();
				thrust::host_vector<unsigned long> h_myclusterCl(cArrayCl);
				thrust::host_vector<float> h_functionValue(2);
				h_functionValue[0] = Cl.getImbalance().getPositiveValue();
				h_functionValue[1] = Cl.getImbalance().getNegativeValue();
				destPositiveImbalance = Cl.getImbalance().getPositiveValue();
				destNegativeImbalance = Cl.getImbalance().getNegativeValue();
				h_nc[0] = Cl.getNumberOfClusters();
				unsigned long numberOfChunks = n * (h_nc[0] + 1);  // the search space for each vertex (dest cluster) will be split into n*(nc+1) chunks
				// A -> Transfer to device
				// transfers the arrays to CUDA device
				d_mycluster.resize(h_nc[0]);
				d_mycluster = h_myclusterCl;
				d_functionValue = h_functionValue;
				clusterArray = thrust::raw_pointer_cast( &d_mycluster[0] );
				funcArray = thrust::raw_pointer_cast( &d_functionValue[0] );
				d_nc = h_nc;
				ncArray = thrust::raw_pointer_cast( &d_nc[0] );
				// recalculates sum matrices
				d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
				d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
				d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
				vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
				vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
				isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
				long blocksPerGrid = (n + threadsCount - 1) / threadsCount;  // , n*sizeof(int)
				updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount>>>(weightArray, destArray, numArray,
					offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
				checkCudaErrors(cudaDeviceSynchronize());

				// VND loop
				int iteration = 0;
				bestImbalance = destPositiveImbalance + destNegativeImbalance;
				while (true) {
					// printf("*** Local search iteration %d, nc = %ld, best I(P) = %.2f\n", iteration, h_nc[0], bestImbalance);
					ulong numberOfChunks = (h_nc[0] + 1) * n;
					// result / destination vector
					d_destFunctionValue.resize(numberOfChunks, 0.0);
					//d_destPosFunctionValue.resize(numberOfChunks, 0.0);
					//d_destNegFunctionValue.resize(numberOfChunks, 0.0);
					float* destImbArray = thrust::raw_pointer_cast( &d_destFunctionValue[0] );
					float* destPosImbArray = thrust::raw_pointer_cast( &d_destPosFunctionValue[0] );
					float* destNegImbArray = thrust::raw_pointer_cast( &d_destNegFunctionValue[0] );
					d_functionValue[0] = destPositiveImbalance;
					d_functionValue[1] = destNegativeImbalance;
					
					// create two dimensional 4x4 thread blocks
					// dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
					dim3 block(BLOCK_SIZE, 1, 1);

					// configure a two dimensional grid as well
					dim3 grid_size;
					long num_elements_x = (numberOfChunks + (block.x) - 1) / block.x;
					long num_elements_y = 1;
					if(num_elements_x > 65535) {
						num_elements_y = (num_elements_x + 65535 - 1) / 65535;
						num_elements_x = 65535;
					}
					grid_size.x = num_elements_x;// / block.x;
					grid_size.y = num_elements_y;// / block.y;
					//printf("grid_size.x = %d, grid_size.y = %d\n", grid_size.x,  grid_size.y);
					
					// printf("The current number of clusters is %ld and bestImbalance = %.2f\n", h_nc[0], bestImbalance);
					//blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
					//if(numberOfChunks != blocksPerGrid * threadsCount) {
					//	printf("numberOfChunks = %ld, numberOfThreads = %ld\n", numberOfChunks, blocksPerGrid * threadsCount);
					//	if(numberOfChunks > blocksPerGrid * threadsCount) {  printf("ERROR!\n");  }
					//}
					//printf("Invoking kernel with %ld chunks, %ld threads, %ld clusters and %ld vertices \n", ulong(numberOfChunks), ulong(grid_size.x*grid_size.y*block.x*block.y), ulong(h_nc[0]), ulong(n));
					simpleSearchKernel<<<grid_size, block>>>(clusterArray, funcArray,
							ulong(n), ulong(m), destImbArray, destPosImbArray, destNegImbArray, ulong(h_nc[0]), ulong(numberOfChunks),
							vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray);
					//cudaCheckErrors("kernel fail");
					checkCudaErrors(cudaDeviceSynchronize());

					// printf("Begin reduce / post-process...\n");

					thrust::device_vector<float>::iterator iter = thrust::min_element(d_destFunctionValue.begin(), d_destFunctionValue.begin() + numberOfChunks);
					float min_val = *iter;
					// as there may be more than one combination with the same imbalance,
					// chooses one combination at random
					// the algorithm starts the search with a random initial index
					/*
					thrust::host_vector<int> h_result_index(1);
					h_result_index[0] = -1;
					thrust::device_vector<int> d_result_index = h_result_index;
					int* resultIndexArray = thrust::raw_pointer_cast( &d_result_index[0] );
					blocksPerGrid = (numberOfChunks + threadsCount - 1) / threadsCount;
					shuffleBestResult1opt<<< blocksPerGrid, threadsCount >>>(destImbArray, numberOfChunks, bestImbalance, resultIndexArray);
					h_result_index = d_result_index;
					int resultIdx = h_result_index[0];
					float min_val = bestImbalance;
					if(resultIdx >= 0) {
						min_val = d_destFunctionValue[resultIdx];
					} */
					// printf("min val is %.2f\n", min_val);

					if((min_val - bestImbalance < EPS) && (fabs(min_val - bestImbalance) > EPS)) {  // (min_val < bestImbalance) => improvement in imbalance
						bestImbalance = min_val;
						// printf("min val is %.2f\n", min_val);
						// determines the position of the best improvement found in the result vector
						 ulong position = iter - d_destFunctionValue.begin();
						 long resultIdx = position;
						ulong destCluster = resultIdx / n;
						ulong bestSrcVertex = resultIdx % n;
						ulong sourceCluster = d_mycluster[bestSrcVertex];
						//destPositiveImbalance = d_destPosFunctionValue[resultIdx];
						//destNegativeImbalance = d_destNegFunctionValue[resultIdx];
						destPositiveImbalance = bestImbalance;
                                                destNegativeImbalance = 0.0;

						// printf("Max idx is %ld\n", numberOfChunks - 1);
						// printf("Idx = %ld: The best src vertex is %ld to cluster %ld with I(P) = %.2f\n", resultIdx, bestSrcVertex, destCluster, min_val);
						if(bestImbalance < EPS) {  /* printf("WARNING: I(P) < 0 !!!\n"); */  break;  }
						d_old_nc = d_nc;
						old_ncArray = thrust::raw_pointer_cast( &d_old_nc[0] );
						updateClustering1opt<<< 1, 1 >>>(bestSrcVertex, destCluster, destPositiveImbalance, destNegativeImbalance, clusterArray, funcArray, n, ncArray);
						checkCudaErrors(cudaDeviceSynchronize());
						thrust::host_vector<ulong> h_old_nc(1);
						h_old_nc = h_nc;
						h_nc = d_nc;
						// h_functionValue = d_functionValue;
						// printf("After: nc = %d, cluster = %d, imbalance = %.2f\n", h_nc[0], h_mycluster[bestSrcVertex], h_functionValue[0]);

						// TEST IMBALANCE CALCULATION FIXME
						/*
 						thrust::host_vector<unsigned long> h_mycluster = d_mycluster;
                        ClusterArray cArray(h_mycluster.size(), 0);
                        for(long x = 0; x < h_mycluster.size(); x++) {
                                cArray[x] = h_mycluster[x];
                        }
   		                // iterationValue = total;
   		                 // printf("Rebuilding clustering object.\n");
           		         // TODO remontar o clustering sem recalcular a FO
	
                   		Clustering Cltest(cArray, *g, problem); // , destPositiveImbalance, destNegativeImbalance);
                		    // assert(bestImbalance == destPositiveImbalance + destNegativeImbalance);
            	        if(destPositiveImbalance + destNegativeImbalance != Cltest.getImbalance().getValue()) {
            	       	         printf("Warning: imbalance on VND does not match. CUDA I(P) = %.2f, CPU: %.2f\n", bestImbalance, Cltest.getImbalance().getValue());
                    	} */

						// FIXME
                        // recalculates sum matrices
						/*
                                d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
                                d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
                                d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
                                vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
                                vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
                                isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
                                long blocksPerGrid = (n + threadsCount - 1) / threadsCount;  // , n*sizeof(int)
                                updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount>>>(weightArray, destArray, numArray,
                                        offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
                                checkCudaErrors(cudaDeviceSynchronize());
						*/
						
						if(h_nc[0] > h_old_nc[0]) {
							// acrescenta uma nova fileira correpondente a um novo cluster na matriz de soma
							//printf("New cluster. nc = %d. Growing vectors to size %2.f MB.\n", h_nc[0], (n * (h_nc[0]+1) * sizeof(float)) / (1024.0 * 1024.0));
							d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
							d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
							d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
							vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
							vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
							isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
						}
						if(h_nc[0] < h_old_nc[0]) {
							//printf("Will shrink vectors.\n");
							//printf("Idx = %ld: The best src vertex is %ld from cluster %ld to cluster %ld with I(P) = %.2f\n", resultIdx, bestSrcVertex, sourceCluster, destCluster, min_val);						
						}
						// CASO ESPECIAL 2: o cluster k1 foi removido -> parcialmente tratado dentro do kernel anterior
						// h_mycluster = d_mycluster; // retrieves new cluster configuration from GPU
						// CASO ESPECIAL 2: cluster removido
						if(h_nc[0] < h_old_nc[0]) {
							// recalculates sum matrices
                            d_VertexClusterPosSum.resize(n * (h_nc[0]+1), 0.0);
        			        d_VertexClusterNegSum.resize(n * (h_nc[0]+1), 0.0);
                			d_neighbor_cluster.resize(n * (h_nc[0]+1), 0);
                            vertexClusterPosSumArray = thrust::raw_pointer_cast( &d_VertexClusterPosSum[0] );
                            vertexClusterNegSumArray = thrust::raw_pointer_cast( &d_VertexClusterNegSum[0] );
                            isNeighborClusterArray = thrust::raw_pointer_cast( &d_neighbor_cluster[0] );
						}

                        long blocksPerGrid = (n + threadsCount - 1) / threadsCount;  // , n*sizeof(int)
						start_time = timer.elapsed();
                        updateVertexClusterSumArrays<<<blocksPerGrid, threadsCount>>>(weightArray, destArray, numArray,
                                offsetArray, clusterArray, vertexClusterPosSumArray, vertexClusterNegSumArray, isNeighborClusterArray, n, ncArray);
                        checkCudaErrors(cudaDeviceSynchronize());
						timer.stop();
                        end_time = timer.elapsed();
                        double timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
						totalTime += timeSpent;
                        //printf("Time spent on matrix full recalculation: %.2f\n", timeSpent);
                        timer.resume();
						
						// printf("Preparing new VND loop...\n");
						if(bestImbalance <= -EPS)   {  /* printf("WARNING: I(P) <= 0 !!!\n"); */  break;  }
					} else {  // no better result found in neighborhood
						// printf("Breaking VND loop...\n");
						break;
					}
					iteration++;
				}
				/*
				SequentialNeighborhoodSearch neighborhoodSearchSeq;
				VariableNeighborhoodDescent vnd(neighborhoodSearchSeq, construct.getRandomSeed(), 1, true, timeLimit);
				Clustering Clnew = vnd.localSearch(g, Cl, i, problem, timeSpentInILS, 0);
				float bestImbalance2 = Clnew.getImbalance().getValue();
				if(bestImbalance != bestImbalance2)   printf("CUDA: %.2f  CPU: %.2f\n", bestImbalance, bestImbalance2);
				*/

				// 3. Select the best clustring so far
				if((bestImbalance - CStar.getImbalance().getValue() < EPS) && (fabs(bestImbalance - CStar.getImbalance().getValue()) > EPS)) {  // (bestImbalance < CStar.getImbalance().getValue()) => improvement in imbalance
					// printf("Imbalance improved in the current iteration.\n");
					/*
					if(perturbated) {
						printf("Improved because of perturbation: j = %d, perturb = %d\n", j, perturbationLevel);
					} else {
						printf("Improved because of constructive phase.\n");
					} */

					thrust::host_vector<unsigned long> h_mycluster = d_mycluster;
					ClusterArray cArray(h_mycluster.size(), 0);
					for(long x = 0; x < h_mycluster.size(); x++) {
						cArray[x] = h_mycluster[x];
					}
					iterationValue = total;
					// printf("Rebuilding clustering object.\n");
					// TODO remontar o clustering sem recalcular a FO
					Clustering Clnew(cArray, *g, problem); // , destPositiveImbalance, destNegativeImbalance);
					// assert(bestImbalance == destPositiveImbalance + destNegativeImbalance);
					
					if(fabs(destPositiveImbalance + destNegativeImbalance - Clnew.getImbalance().getValue()) > EPS) {  // (destPositiveImbalance + destNegativeImbalance != Clnew.getImbalance().getValue()) => different values
						//printf("Warning: imbalance does not match. CUDA I(P) = %.6f, CPU: %.6f\n", bestImbalance, Clnew.getImbalance().getValue());
						destPositiveImbalance = Clnew.getImbalance().getPositiveValue();
						destNegativeImbalance = Clnew.getImbalance().getNegativeValue();
						bestImbalance = Clnew.getImbalance().getValue();
					}
					// assert( destPositiveImbalance + destNegativeImbalance == Clnew.getImbalance().getValue());
					// printf("I(P) = %.2f\n", Clnew.getImbalance().getValue());
					CStar = Clnew;
					// printf("Clustering object REBUILT.\n");
					// restarts the internal ILS loop and the perturbation
					if((bestImbalance - CStar.getImbalance().getValue() < EPS) && (fabs(bestImbalance - CStar.getImbalance().getValue()) > EPS)) {  // (bestImbalance < CStar.getImbalance().getValue()) => improvement in imbalance
						j = 1;
						perturbationLevel = 1;
						timeSpentOnBestSolution = timeSpentInILS;
					} else {
						j++;
                    	if(j > iterMaxILS) {
                            // printf("Increasing perturbation to level %d...\n", perturbationLevel + 1);
                            perturbationLevel++;
                            j = 1;
                            if(perturbationLevel > perturbationLevelMax) {
                                    break;
                            }
                    	}
					}
					if(bestImbalance < EPS)  break;
				} else {  // did not improve solution
					j++;
					if(j > iterMaxILS) {
						// printf("Increasing perturbation to level %d...\n", perturbationLevel + 1);
						perturbationLevel++;
						j = 1;
						if(perturbationLevel > perturbationLevelMax) {
							break;
						}
					}
				}
				// 4. Generate perturbation over C* -> kernel perturbation
				// Cl = perturbation.randomMove(g, CStar, problem, perturbationLevel);
				Cl = CStar;

				long n = g->getN();
				long nc = Cl.getNumberOfClusters();
				/*
				std::vector<int> nodeList(n, 0);
				for(int i = 0; i < n; i++) {
					nodeList[i] = i;
				}
				std::random_shuffle(nodeList.begin(), nodeList.end());
				std::vector<int> nodeList; */
				//generateRandomNumbers(perturbationLevel, 0, n - 1, nodeList);
				/* Allocate n floats on device */
				d_random.resize(perturbationLevel, 0);
				float* devData = thrust::raw_pointer_cast( &d_random[0] );
				/* Generate n floats on device */
				curandGenerateUniform(gen, devData, perturbationLevel);
				h_random = d_random;
				curandGenerateUniform(gen2, devData, perturbationLevel);
				h_random2 = d_random;

				for(int i = 0; i < perturbationLevel; i++) {
					nc = Cl.getNumberOfClusters();
					// int k2 = randomUtil.next(0, nc - 1);
					// int idx = randomUtil.next(0, n - 1);
					// cout << "k2 = " << k2 << endl;
					int rnd_integer_range = 0 + h_random[i] * (n - 1 - 0);
					int k2 = 0 + h_random2[i] * (nc - 1 - 0);
					Cl = perturbation.move1optCCProblem(g, Cl, problem, rnd_integer_range, k2);
				}
				// printf("Perturbated I(P) = %.2f\n", Cl.getImbalance().getValue());
				// printf("After perturb: nc = %d, imbalance = %.2f\n", h_nc[0], h_functionValue[0]);
				perturbated = true;
				
				// 5. Stops the timer and stores the elapsed time
				timer.stop();
				end_time = timer.elapsed();

				// 6. Write the results into ostream os, using csv format
				// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
				timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
				timer.resume();
				start_time = timer.elapsed();
				// if elapsed time is bigger than timeLimit, break
				if(timeSpentInILS >= timeLimit) {
					// BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
					break;
				}
			}
			Imbalance newValue = CStar.getImbalance();
			Imbalance bestValue = CBest.getImbalance();
			double newImbValue = newValue.getValue();
			double bestImbValue = bestValue.getValue();
			if((newImbValue - bestImbValue < EPS) && (fabs(newImbValue - bestImbValue) > EPS)) {  // (newValue < bestValue) => improvement in imbalance
				//printf("A better global solution was found: I(P) = %.5f\n", newValue.getValue());
				CBest = CStar;
				bestValue = newValue;
				iterationValue = i;
				timeSpentOnBestSolution = timeSpentInILS;
				if(newValue.getValue() == 0)  break;
			}
			timer.stop();
			end_time = timer.elapsed();
			timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);

			iterationResults << (i+1) << "," << bestValue.getValue() << "," << bestValue.getPositiveValue()
					<< "," << bestValue.getNegativeValue() << "," << CBest.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInILS << "\n";
			printf("total time spent on matrix recalc: %.2f\n", totalTime);
			// if elapsed time is bigger than timeLimit, break
			if(timeSpentInILS >= timeLimit) {
				// BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
				break;
			}
			timer.resume();
			start_time = timer.elapsed();
			// guarantees at least one execution of the ILS when the number of iterations is smaller than one
			if(iterMax <= 0) {  break;  }

			// Avoids constructClustering if loop break condition is met
			if((i + 1) < iterMax) {
				// 0. Triggers local processing time calculation
				start_time = timer.elapsed();
				// 1. Construct the next clustering
				Cc = construct.constructClustering(g, problem, processRank);
				
				timer.stop();
				end_time = timer.elapsed();
				timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
				timeSpentInILS += timeSpentInConstruction;
				timer.resume();
				start_time = timer.elapsed();
			}
		}

		curandDestroyGenerator(gen);
		curandDestroyGenerator(gen2);

		printf("I(P) = %.5f\n", CBest.getImbalance().getValue());
		result = CBest;
		totalIterations = iterMax;
		timeSpentConstruct = timeSpentInConstruction;
		timeSpentILS = timeSpentInILS;
		return true;
	}
}
