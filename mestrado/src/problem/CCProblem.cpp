/*
 * CCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/CCProblem.h"
#include <cmath>
#include <boost/log/trivial.hpp>

namespace problem {

using namespace std;
using namespace clusteringgraph;
using namespace boost;

CCProblem::CCProblem() {
	// TODO Auto-generated constructor stub

}

CCProblem::~CCProblem() {
	// TODO Auto-generated destructor stub
}

int CCProblem::getType() const {
	return ClusteringProblem::CC_PROBLEM;
}

/**
 * The Imbalance of a partition P (I(P)) is defined as the total
 * weight of negative uncut arcs and positive cut arcs.
 */
Imbalance CCProblem::objectiveFunction(SignedGraph& g, const ClusterList& c) const {
	double positiveSum = 0, negativeSum = 0;
	int nc = c.size();
	int n = g.getN();

	BOOST_LOG_TRIVIAL(trace) << "[CCProblem] Disparando calculo da funcao objetivo." << endl;

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		// Find out to which cluster vertex i belongs
		for(unsigned long k = 0; k < nc; k++) {
			cluster = c.at(k);
			if(cluster[i]) {
				break;
			}
		}
		DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g.graph);
			bool sameCluster = cluster[j];
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (positiveSum + negativeSum) << endl;
	return Imbalance(positiveSum, negativeSum);
}

// Calculates the delta of the objective function
Imbalance CCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, const ClusterList& c,
			const unsigned long& k, const unsigned long& i) const {
	double negativeSum = 0, positiveSum = 0;
	BoolArray cluster = c.at(k);
	// unsigned long n = g.getN();

	// iterates over out-edges of vertex i
	//typedef graph_traits<Graph> GraphTraits;
	DirectedGraph::out_edge_iterator out_i, out_end;
	DirectedGraph::edge_descriptor e;

	// std::cout << "out-edges of " << i << ": ";
	for (tie(out_i, out_end) = out_edges(i, g.graph); out_i != out_end; ++out_i) {
		e = *out_i;
		Vertex src = source(e, g.graph), targ = target(e, g.graph);
		double weight = ((Edge*)out_i->get_property())->weight;
		if(cluster[targ.id]) {
			if(weight < 0) {
				negativeSum += fabs(weight);
			}
		} else {
			if(weight > 0) {
				positiveSum += weight;
			}
		}
		// std::cout << "(" << src.id << "," << targ.id << ") ";
	}

	// iterates over in-edges of vertex i
	DirectedGraph::in_edge_iterator in_i, in_end;
	// std::cout << "in-edges of " << i << ": ";
	for (tie(in_i, in_end) = in_edges(i, g.graph); in_i != in_end; ++in_i) {
		e = *in_i;
		Vertex src = source(e, g.graph), targ = target(e, g.graph);
		double weight = ((Edge*)in_i->get_property())->weight;
		if(cluster[src.id]) {
			if(weight < 0) {
				negativeSum += fabs(weight);
			}
		} else {
			if(weight > 0) {
				positiveSum += weight;
			}
		}
		// std::cout << "(" << src.id << "," << targ.id << ") ";
	}
	return Imbalance(positiveSum, negativeSum);
}

} /* namespace problem */
