/*
 * CUDANeighborhoodSearch.cpp
 *
 *  Created on: May 25, 2014
 *      Author: mlevorato
 */

#include "include/CUDANeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include "include/CUDASearch.h"
#include "util/include/RandomUtil.h"

#include <boost/log/trivial.hpp>
#include <boost/timer/timer.hpp>

#include <limits>
#include <thrust/host_vector.h>


namespace clusteringgraph {

using namespace thrust;
using namespace util;


CUDANeighborhoodSearch::CUDANeighborhoodSearch() {
	// TODO Auto-generated constructor stub

}

CUDANeighborhoodSearch::~CUDANeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

Clustering CUDANeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
        Clustering* clustering, ClusteringProblem& problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, bool firstImprovementOnOneNeig) {

	// Resets the number of combinations tested on neighborhood search
	numberOfTestedCombinations = 0;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
	}

	unsigned long threadsCount = 120;

	// Splits the processing in (n / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = 0;
	unsigned long remainingVertices = 0;
	if(g->getN() > threadsCount) {
		sizeOfChunk = (unsigned long)((double)g->getN() / threadsCount);
		remainingVertices = g->getN() % threadsCount;
	} else {
		sizeOfChunk = 1;
		threadsCount = g->getN();
		remainingVertices = 0;
	}
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	// invoke the 1-opt CUDA method
	BOOST_LOG_TRIVIAL(debug) << "Waiting for CUDA...\n";
	Clustering bestClustering;
	double bestValue = numeric_limits<double>::infinity();
	bestClustering = this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed, myRank, 0, threadsCount*sizeOfChunk - 1,
			firstImprovementOnOneNeig, k);

	// the leader (me) does its part of the work too
	BOOST_LOG_TRIVIAL(trace) << "Total number of vertices is " << g->getN() << endl;
	BOOST_LOG_TRIVIAL(trace) << "CUDA Parallelization status: " << threadsCount <<
			" search slaves (threads) will process " << sizeOfChunk << " vertice(s) each one." << endl;
	BOOST_LOG_TRIVIAL(trace) << "RemainingVertices is " << remainingVertices << endl;

	if(remainingVertices > 0) {
		bestClustering = NeighborhoodSearch::search1opt(g, clustering,
				problem, timeSpentSoFar, timeLimit, randomSeed, myRank,
				threadsCount * sizeOfChunk, threadsCount * sizeOfChunk + remainingVertices - 1, firstImprovementOnOneNeig, k);
		bestValue = bestClustering.getImbalance().getValue();
	}
	// the leader receives the processing results, if that is the case
	// TODO implement global first improvement (l = 2-opt) in parallel VND here!
	if(sizeOfChunk > 0) {
		/*
		OutputMessage omsg;
		for(i = 0; i < numberOfSearchSlaves; i++) {
			mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_VND_TAG, omsg);
			BOOST_LOG_TRIVIAL(debug) << "Message received from process " << stat.source() << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			BOOST_LOG_TRIVIAL(trace) << omsg.clustering.toString();
			// processes the result of the execution of process p(i)
			// sums the number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// checks if the obj value improved
			if(omsg.clustering.getImbalance().getValue() < bestValue) {
				Clustering clustering = omsg.clustering;
				bestClustering = clustering;
				bestValue = omsg.clustering.getImbalance().getValue();
				BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VND] Better value found for objective function in node "
						<< stat.source() << ": " <<
						omsg.clustering.getImbalance().getValue() << endl;
				// IMPORTANT: Terminate other slaves' VND search if 2-opt and CC Problem
				if( (l == 2) and (problem.getType() == ClusteringProblem::CC_PROBLEM) ) {
					InputMessageParallelVND imsg;
					BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VND] First improvement on 2-opt: interrupting other VND slaves.";
					for(int j = 0; j < numberOfSearchSlaves; j++) {
						world.send(slaveList[j], MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG, imsg);
					}
				}
			}
		}
		*/
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel Search CUDA] Best solution found: Obj = " << bestClustering.getImbalance().getValue() << endl;
	bestClustering.printClustering();
	return bestClustering;
}

Clustering CUDANeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialSearchIndex,
		unsigned long finalSearchIndex, bool firstImprovementOnOneNeig, unsigned long k) {
	// Unused method
	assert(initialSearchIndex < g->getN());
	assert(finalSearchIndex < g->getN());

	if (l == 1) {  // 1-opt
		// Parallel search always does best improvement in 1-opt
		// Therefore, parameter firstImprovementOnOneNeig is ignored
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, false, k);
	} else {  // 2-opt
		// TODO implement global first improvement in parallel 2-opt
		// first improvement found in one process must break all other processes' loop
		// IMPORTANT: Parallel VND does first improvement on 2-opt if CC Problem
		// Or best improvement on 2-opt if RCC Problem
		bool firstImprovementOn2Opt = true;
		if(problem.getType() == ClusteringProblem::CC_PROBLEM) {
			firstImprovementOn2Opt = true;
		} else if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			// TODO check if this is still necessary with GRASP RCC
			// firstImprovementOn2Opt = false;
		}
		return this->search2opt(g, clustering, &problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, firstImprovementOn2Opt, k);
	}
}

Clustering CUDANeighborhoodSearch::search1opt(SignedGraph* g,
                Clustering* clustering, ClusteringProblem& problem,
                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                int myRank, unsigned long initialSearchIndex,
        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k) {
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	unsigned long n = g->getN();
	unsigned long m = g->getM();
	unsigned long numberOfVerticesInInterval = finalSearchIndex - initialSearchIndex + 1;
	unsigned long totalNumberOfClusters = clustering->getNumberOfClusters();
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	// stores the best clustering combination generated (minimum imbalance) - used by 1-opt neighborhood
	Clustering cBest = *clustering;
	// pre-calculates, in an array, to which cluster each vertex belongs
	boost::shared_ptr<unsigned long[]> myCluster(new unsigned long[n]);
	for(unsigned long k = 0; k < totalNumberOfClusters; k++) {  // for each cluster k
		BoolArray clusterK = clustering->getCluster(k);
		for(unsigned long i = 0; i < n; i++) {  // for each vertex i
			if(clusterK[i]) {  // vertex i is in cluster k
				myCluster[i] = k;
			}
		}
	}

	BOOST_LOG_TRIVIAL(info) << "[CUDA] Begin transfer to device...";
	// A -> Transfer to device
	// transfers the myClusters array to CUDA device
	int numberOfThreads = finalSearchIndex - initialSearchIndex;
	unsigned long* cluster = myCluster.get();
	thrust::host_vector<unsigned long> h_mycluster(cluster, cluster + n);
	thrust::host_vector<unsigned long> h_destcluster(cluster, cluster + n);
	// objective function value
	thrust::host_vector<float> h_functionValue(2);
	h_functionValue[0] = clustering->getImbalance().getPositiveValue();
	h_functionValue[1] = clustering->getImbalance().getNegativeValue();
	thrust::host_vector<float> h_destPosFunctionValue(numberOfThreads);
	thrust::host_vector<float> h_destNegFunctionValue(numberOfThreads);
	// graph structure (adapted adjacency list)
	thrust::host_vector<float> h_weights(m);  // edge weights
	thrust::host_vector<int> h_dest(m);  // edge destination (vertex j)
	thrust::host_vector<int> h_numedges(n);  // number of edges of each vertex i
	thrust::host_vector<int> h_offset(n);  // initial edge number for vertex i
	// For each vertex i
	for(int i = 0, offset = 0, edge = 0; i < n; i++) {
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		h_offset[i] = offset;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		h_numedges[i] = count;
		offset += count;
	}
	// TODO transform into class constant
	unsigned short threadsCount = 512;

	// Pass raw array and its size to kernel
	runSimpleSearchKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, n, m,
			h_destcluster, h_destPosFunctionValue, h_destNegFunctionValue, threadsCount, totalNumberOfClusters);

	BOOST_LOG_TRIVIAL(info) << "[CUDA] Before: " << h_functionValue[0] << "; " << h_functionValue[1];
	BOOST_LOG_TRIVIAL(info) << "[CUDA] Processing complete. " << h_destPosFunctionValue[2] << "; " << h_destPosFunctionValue[4];

	// Process result vectors


	/*
	// for each vertex i, tries to move i to another cluster in myNeighborClusterList[i]
	// For each node i in cluster(k1)
	for (unsigned long i = randomUtil.next(initialSearchIndex, finalSearchIndex), cont2 = 0; cont2 < numberOfVerticesInInterval; cont2++) {
		// vertex i is in cluster(k1)
		unsigned long k1 = myCluster[i];
		// RCC Problem: must have exactly k clusters -> Cannot remove node from standalone cluster
		unsigned long s = clustering->getClusterSize(k1);
		if((problem.getType() == ClusteringProblem::RCC_PROBLEM) && (s == 1)) {
			continue;
		}
		// Option 1: node i is moved from k1 to another existing cluster k2 != k1
		for (unordered_set<unsigned long>::iterator itr = myNeighborClusterList[i].begin(); itr != myNeighborClusterList[i].end(); ++itr) {
			// cluster(k2)
			unsigned long k2 = *itr;
			// removes node i from cluster1 and inserts in cluster2
			Clustering cTemp = *clustering;
			BoolArray cluster2 = cTemp.getCluster(k2);

			//BOOST_LOG_TRIVIAL(trace) << "Option 1: Taking node " << i << " from cluster " << k1 << " to cluster " << k2;
			int nc = cTemp.getNumberOfClusters();
			cTemp.removeNodeFromCluster(*g, problem, i, k1);
			// recalculates the number of clusters, as one of them may have been removed
			if((cTemp.getNumberOfClusters() < nc) && (k2 >= k1)) {
				// cluster k1 has been removed
				cTemp.addNodeToCluster(*g, problem, i, k2 - 1);
			} else {
				cTemp.addNodeToCluster(*g, problem, i, k2);
			}
			numberOfTestedCombinations++;
			// cTemp->printClustering();
			Imbalance newImbalance = cTemp.getImbalance();
			Imbalance bestImbalance = cBest.getImbalance();
			if (newImbalance < bestImbalance) {
				//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 1-neighborhood: " << setprecision(2) << newImbalance.getValue();
				// First improvement for 1-opt neighborhood
				cBest = cTemp;
				if(firstImprovement) {
					return cBest;
				}
			}
			// return if time limit is exceeded
			boost::timer::cpu_times end_time = timer.elapsed();
			double localTimeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			// std::cout << timeSpentSoFar + localTimeSpent << endl;
			if(timeSpentSoFar + localTimeSpent >= timeLimit)  return cBest;
		}
		// Option 2: node i is moved to a new cluster, alone
		// removes node i from cluster k1 and inserts in newCluster
		// cout << "New clustering combination generated." << endl;
		if((problem.getType() == ClusteringProblem::CC_PROBLEM) && (s > 1)) {
			// this code is not executed for RCC Problem, as it must have exactly k clusters
			// this code is not executed if node i is to be moved from a standalone cluster to another standalone cluster (s == 1)
			Clustering cTemp = *clustering;
			//BOOST_LOG_TRIVIAL(trace) << "Option 2: Taking node " << i << " from " << k1 << " to new cluster.";
			cTemp.removeNodeFromCluster(*g, problem, i, k1);
			BoolArray cluster2 = cTemp.addCluster(*g, problem, i);
			numberOfTestedCombinations++;
			// cTemp->printClustering();
			Imbalance newImbalance = cTemp.getImbalance();
			Imbalance bestImbalance = cBest.getImbalance();
			if (newImbalance < bestImbalance) {
				//BOOST_LOG_TRIVIAL(trace) << "Better solution found in 1-neighborhood: " << setprecision(2) << newImbalance.getValue();
				// First improvement for 1-opt neighborhood
				cBest = cTemp;
				if(firstImprovement) {
					return cBest;
				}
			}
		}
		// increment rule
		i++;
		if(i > finalSearchIndex) {
			i = initialSearchIndex;
		}
	}*/
	// returns the best combination found in 1-opt
	return cBest;
}

} /* namespace clusteringgraph */
