/*
 * CUDAConstructClustering.cpp
 *
 *  Created on: 09/02/2015
 *      Author: czt0
 */

#include "include/CUDAConstructClustering.h"
#include "include/ConstructClustering.h"

#include "../construction/include/VertexSet.h"
#include <boost/log/trivial.hpp>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>

#include <thrust/host_vector.h>

using namespace thrust;

namespace resolution {
namespace construction {

CUDAConstructClustering::CUDAConstructClustering(SignedGraph *g, const unsigned long& seed) :
		ConstructClustering(CUDAImbalanceGainFunction(g), seed, double(1.0)) {

}

CUDAConstructClustering::~CUDAConstructClustering() {

}


Clustering ConstructClustering::constructClustering(SignedGraph *g,
		ClusteringProblem& problem, const int& myRank) {
	Clustering Cc(*g); // Cc = empty
	VertexSet lc(randomSeed, g->getN()); // L(Cc) = V(G)
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

	while (lc.size() > 0) { // lc != empty
		GainCalculation gainCalculation;
		if (_alpha == 1.0) {
			// alpha = 1.0 (completely random): no need to calculate all gains (saves time)
			int i = lc.chooseRandomVertex(lc.size()).vertex;
			gainCalculation = gainFunction->calculateIndividualGain(problem, Cc,
					i, h_weights, h_dest, h_numedges, h_offset);
		} else {
			// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
			// according to the value of the gain function
			gainFunction->calculateGainList(problem, Cc, lc.getVertexList(),
					h_weights, h_dest, h_numedges, h_offset);

			// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
			// (alpha x |lc|) is a rounded number
			// IMPORTANT: if alpha < 0, the constructClustering method will always
			//             choose the first vertex in the gainFunction list, that is,
			//             the one that minimizes the objective (VOTE algorithm).
			if (_alpha <= 0.0) {
				// chooses first element (min) without ordering (saves time)
				gainCalculation = lc.chooseFirstElement(gainFunction);
			} else {
				lc.sort(gainFunction);
				gainCalculation = lc.chooseRandomVertex(
						boost::math::iround(_alpha * lc.size()));
			}
			// std::cout << "Random vertex between 0 and " << boost::math::iround(alpha * lc.size()) << " is " << i << std::endl;
		}
		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		// cout << "Selected vertex is " << i << endl;
		// The gain function ensures no clusters of size > k are created if RCC Problem
		if (gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			Cc.addCluster(*g, problem, gainCalculation.vertex);
		} else {
			// inserts i into existing cluster
			Cc.addNodeToCluster(*g, problem, gainCalculation.vertex,
					gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the choosing vertex i automatically removes it from the list
		// Removal already done by the chooseVertex methods above
		// Cc->printClustering();
	}
	// Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	// => Finally: Stops the timer and stores the elapsed time
	Cc.setImbalance(problem.objectiveFunction(*g, Cc));
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
