/*
 * RCCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/RCCProblem.h"

#include <boost/scoped_array.hpp>
#include <cmath>
#include <algorithm>
#include <boost/log/trivial.hpp>

using namespace boost;

namespace problem {

RCCProblem::RCCProblem() {
	// TODO Auto-generated constructor stub

}

RCCProblem::~RCCProblem() {
	// TODO Auto-generated destructor stub
}

int RCCProblem::getType() const {
	return ClusteringProblem::RCC_PROBLEM;
}

// TODO: implement this
/**
 * Returns he Relaxed Imbalance of a partition P (RI(P)).
 */
Imbalance RCCProblem::objectiveFunction(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	c.positiveSum.resize(nc, nc);
	c.negativeSum.resize(nc, nc);

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo da funcao objetivo." << endl;

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		unsigned long ki;
		// Find out to which cluster vertex i belongs
		for(ki = 0; ki < nc; ki++) {
			cluster = c.getCluster(ki);
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
			if(sameCluster) {
				if(weight > 0) { // positive edge
					c.positiveSum(ki, ki) += weight;
				} else { // negative edge
					c.negativeSum(ki, ki) += fabs(weight);
				}
			} else {
			    // Find out to which cluster vertex j belongs
				unsigned long kj = 0;
				for(kj = 0; kj < nc; kj++) {
					cluster = c.getCluster(kj);
					if(cluster[j]) {
						break;
					}
				}
				if(weight > 0) { // positive edge
					c.positiveSum(ki, kj) += weight;
				} else { // negative edge
					c.negativeSum(ki, kj) += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < k1; k2++) {
			externalSum += min(c.positiveSum(k1, k2), c.negativeSum(k1, k2));
		}
	}
	
	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	return Imbalance(internalSum, externalSum);
}

// Calculates the delta of the objective function
// TODO: Separar os casos de inclusao e exclusao de vertice atraves de parametro.
Imbalance RCCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, Clustering& c,
		const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo do delta da funcao objetivo." << endl;

	BoolArray currentCluster;
	unsigned long ki;
	// Find out to which cluster vertex i belongs (currentCluster)
	for(ki = 0; ki < nc; ki++) {
		currentCluster = c.getCluster(ki);
		if(currentCluster[i]) {
			break;
		}
	}
	// gets vertex i's new cluster
	BoolArray newCluster = c.getCluster(k);

	DirectedGraph::out_edge_iterator f, l;
	// For each out edge of i
	for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
		double weight = ((Edge*)f->get_property())->weight;
		int j = target(*f, g.graph);
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = currentCluster[j];
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) -= weight;
			} else { // negative edge
				c.negativeSum(ki, kj) -= fabs(weight);
			}
		}
		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		sameCluster = newCluster[j];
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) += weight;
			} else { // negative edge
				c.negativeSum(ki, kj) += fabs(weight);
			}
		}
	}
	// recalculates the imbalance based on the matrices
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < k1; k2++) {
			externalSum += min(c.positiveSum(k1, k2), c.negativeSum(k1, k2));
		}
	}

	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	return Imbalance(internalSum, externalSum);
	// TODO marcar diferenças entre antes e depois para calcular apenas o que mudou nas matrizes
}

} /* namespace problem */
