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
#include <boost/numeric/ublas/matrix.hpp>

namespace problem {

using namespace std;
using namespace clusteringgraph;
using namespace boost;
using namespace boost::numeric::ublas;

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
	ClusterArray myCluster = c.getClusterArray();

	BOOST_LOG_TRIVIAL(trace) << "[CCProblem] Starting obj function calculation.";

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			double weight = ((Edge*)f->get_property())->weight;
			long j = target(*f, g.graph);
			long kj = myCluster[j];
			bool sameCluster = (ki == kj);
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
	ClusterArray myCluster = c.getClusterArray();
	long ki = myCluster[i];
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
		if(myCluster[targ.id] == ki) {
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
		if(myCluster[src.id] == ki) {
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
	ClusterArray myCluster = c.getClusterArray();
	stringstream ss1, ss2, ss3;
	DirectedGraph::edge_descriptor e;
	// Cluster to cluster matrix containing positive / negative contribution to imbalance
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);

	BOOST_LOG_TRIVIAL(info) << "[CCProblem] Starting imbalance analysis.";
	ss1 << endl << "Imbalance analysis (out edges contribution):" << endl;
	ss1 << "Vertex,PositiveSum,NegativeSum" << endl;
	ss2 << "Imbalance analysis (in edges contribution):" << endl;
	ss2 << "Vertex,PositiveSum,NegativeSum" << endl;

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		DirectedGraph::out_edge_iterator f, l;
		double positiveSum = 0, negativeSum = 0;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ((Edge*)f->get_property())->weight;
			bool sameCluster = (myCluster[targ.id] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
				clusterImbMatrix(ki, myCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
				clusterImbMatrix(ki, myCluster[targ.id]) += fabs(weight);
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
			bool sameCluster = (myCluster[src.id] == ki);
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

	ss3 << endl << "Cluster contribution to imbalance analysis (cluster-cluster matrix):" << endl;
	for(int i = 0; i < nc; i++) {
		for(int j = 0; j < nc; j++) {
			ss3 << clusterImbMatrix(i, j) << ", ";
		}
		ss3 << endl;
	}
	
	BOOST_LOG_TRIVIAL(info) << "[CCProblem] Graph analysis done." << endl;
	return ss1.str() + ss2.str() + ss3.str();
}
} /* namespace problem */
