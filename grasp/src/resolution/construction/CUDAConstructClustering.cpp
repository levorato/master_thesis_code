/*
 * CUDAConstructClustering.cpp
 *
 *  Created on: 09/02/2015
 *      Author: czt0
 */

#include "include/CUDAConstructClustering.h"
#include "include/CUDAImbalanceGainFunction.h"
#include "../vnd/include/CUDASearch.h"
#include "../construction/include/VertexSet.h"

#include <boost/log/trivial.hpp>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>

#include <thrust/host_vector.h>

using namespace thrust;

namespace resolution {
namespace construction {

CUDAConstructClustering::CUDAConstructClustering(CUDAImbalanceGainFunction *f, const unsigned long& seed) :
		ConstructClustering(f, seed, double(1.0)) {

}

CUDAConstructClustering::~CUDAConstructClustering() {

}


Clustering CUDAConstructClustering::constructClustering(SignedGraph *g,
		ClusteringProblem& problem, const int& myRank) {
	Clustering Cc(*g); // Cc = empty
	BOOST_LOG_TRIVIAL(debug)<< "CUDA parallel construct clustering...\n";
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// CUDA variables
	unsigned long n = g->getN();
	unsigned long m = g->getM();
	// CUDA graph structure (adapted adjacency list)
	thrust::host_vector<float> h_weights(2 * m);  // in/out edge weights
	thrust::host_vector<int> h_dest(2 * m);  // edge destination (vertex j)
	thrust::host_vector<int> h_numedges(n);  // number of edges of each vertex i
	thrust::host_vector<int> h_offset(n);  // initial edge number for vertex i
	// For each vertex, creates a list of in and out edges
	int i = 0, offset = 0;
	for(int edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		h_offset[i] = offset;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		h_numedges[i] = count;
		offset += count;
	}

	// objective function value
	thrust::host_vector<float> h_functionValue(1);
	h_functionValue[0] = c.getImbalance().getValue();
	thrust::host_vector<unsigned long> h_mycluster(c.getClusterArray());
	ulong nc = c.getNumberOfClusters();
	for(int e = 0; e < h_mycluster.size(); e++) {
		if (h_mycluster[e] == Clustering::NO_CLUSTER) {
			h_mycluster[e] = nc;
		}
	}
	thrust::host_vector<unsigned long> h_newcluster;
	double imbalance = 0.0;
	runConstructKernel(randomSeed, h_weights, h_dest, h_numedges, h_offset, h_mycluster, h_functionValue, graph->getN(),
			graph->getM(), c.getNumberOfClusters(), threadsCount, h_newcluster, imbalance);

	// recreates clustering object based on cluster array
	ClusterArray cArray;
	for(int x = 0; x < h_newcluster.size(); x++) {
		cArray.push_back(h_newcluster[x]);
	}
	Clustering newCc(cArray, *g, problem);
	// Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	// => Finally: Stops the timer and stores the elapsed time
	newCc.setImbalance(problem.objectiveFunction(*g, newCc));
	assert(imbalance == newCc.getImbalance().getValue());
	BOOST_LOG_TRIVIAL(debug)<< "Initial clustering completed. Obj = " << Cc.getImbalance().getValue();
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpent = (end_time.wall - start_time.wall)
					 / double(1000000000);

	// Cc.printClustering();
	return Cc;
}


} /* namespace construction */
} /* namespace resolution */
