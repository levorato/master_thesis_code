/*
 * CCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/CCProblem.h"
#include <cmath>
#include <boost/log/trivial.hpp>
#include "../graph/include/Clustering.h"

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

string CCProblem::getName() {
	return "CC";
}

/**
 * The Imbalance of a partition P (I(P)) is defined as the total
 * weight of negative uncut arcs and positive cut arcs.
 */
Imbalance CCProblem::objectiveFunction(SignedGraph& g, Clustering& c) {
	double positiveSum = 0, negativeSum = 0;
	int nc = c.getNumberOfClusters();
	int n = g.getN();

	BOOST_LOG_TRIVIAL(trace) << "[CCProblem] Starting obj function calculation.";

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		// Find out to which cluster vertex i belongs
		for(unsigned long k = 0; k < nc; k++) {
			cluster = c.getCluster(k);
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
	BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (positiveSum + negativeSum) << endl;
	return Imbalance(positiveSum, negativeSum);
}

Imbalance CCProblem::calculateDeltaPlusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i) {
	return calculateDeltaObjectiveFunction(g, c, k, i);
}

Imbalance CCProblem::calculateDeltaMinusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i) {
	return calculateDeltaObjectiveFunction(g, c, k, i);
}

// Calculates the delta of the objective function
Imbalance CCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i) {
	double negativeSum = 0, positiveSum = 0;
	BoolArray cluster = c.getCluster(k);
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

/**
  * Contabilizar, para cada vertice, o total de arestas que estao contribuindo para o 
  * desequilibrio (positivo + negativo). Assim a gente vai poder comparar quem tem 
  * mais relacoes em desequilibrio e quem sabe tirar alguma conclusao disso.
  * Esse script será executado sobre os resultados de todos os grafos de intâncias reais.
*/
string CCProblem::analyzeImbalance(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	stringstream ss1, ss2;
	DirectedGraph::edge_descriptor e;

	BOOST_LOG_TRIVIAL(info) << "[CCProblem] Starting imbalance analysis.";
	ss1 << endl << "Imbalance analysis (out edges contribution):" << endl;
	ss1 << "Vertex,PositiveSum,NegativeSum" << endl;
	ss2 << "Imbalance analysis (in edges contribution):" << endl;
	ss2 << "Vertex,PositiveSum,NegativeSum" << endl;

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		// Find out to which cluster vertex i belongs
		for(unsigned long k = 0; k < nc; k++) {
			cluster = c.getCluster(k);
			if(cluster[i]) {
				break;
			}
		}
		DirectedGraph::out_edge_iterator f, l;
		double positiveSum = 0, negativeSum = 0;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ((Edge*)f->get_property())->weight;
			bool sameCluster = cluster[targ.id];
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
			}
		}
		ss1 << i << "," << positiveSum << "," << negativeSum << endl;
		
		DirectedGraph::in_edge_iterator f2, l2;
		positiveSum = 0, negativeSum = 0;
		// For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {
			e = *f2;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ((Edge*)f2->get_property())->weight;
			bool sameCluster = cluster[src.id];
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
			}
		}
		ss2 << i << "," << positiveSum << "," << negativeSum << endl;
	}
	
	BOOST_LOG_TRIVIAL(info) << "[CCProblem] Graph analysis done." << endl;
	return ss1.str() + ss2.str();
}
} /* namespace problem */
