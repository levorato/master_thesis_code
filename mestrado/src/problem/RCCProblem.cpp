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
Imbalance RCCProblem::objectiveFunction(SignedGraph& g, const ClusterList& c) const {
	int nc = c.size();
	int n = g.getN();
	scoped_array<double> internalPositiveSum(new double[nc]);
	scoped_array<double> internalNegativeSum(new double[nc]);
	// os valores de soma entre clusters devem compor uma matriz
	scoped_array<double> externalPositiveSum(new double[nc]);
	scoped_array<double> externalNegativeSum(new double[nc]);
	std::fill_n(internalPositiveSum.get(), nc, 0.0);
	std::fill_n(internalNegativeSum.get(), nc, 0.0);
	std::fill_n(externalPositiveSum.get(), nc, 0.0);
	std::fill_n(externalNegativeSum.get(), nc, 0.0);

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo da funcao objetivo." << endl;

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		unsigned long k;
		// Find out to which cluster vertex i belongs
		for(k = 0; k < nc; k++) {
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
			if(sameCluster) {
				if(weight > 0) { // positive edge
					externalPositiveSum[k][k] += weight;
				} else { // negative edge
					externalNegativeSum[k][k] += fabs(weight);
				}
			} else {
			    // Find out to which cluster vertex j belongs
				int k2 = 0;
				for(k2 = 0; k2 < nc; k2++) {
					cluster = c.at(k2);
					if(cluster[j]) {
						break;
					}
				}
				if(weight > 0) { // positive edge
					externalPositiveSum[k][k2] += weight;
				} else { // negative edge
					externalNegativeSum[k][k2] += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	for(int k = 0; k < nc; k++) {
		internalSum += min(externalPositiveSum[k][k], externalNegativeSum[k][k]);
		for(int k2 = 0; k2 < k; k2++) {
			externalSum += min(externalPositiveSum[k][k2], externalNegativeSum[k][k2]);
		}
	}
	
	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	return Imbalance(internalSum, externalSum);
}

// Calculates the delta of the objective function
Imbalance RCCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, const ClusterList& c,
		const unsigned long& k, const unsigned long& i) const {
	double negativeSum = 0, positiveSum = 0;
	// Not implemented.
	// TODO marcar diferenças entre antes e depois para calcular apenas o que mudou nas matrizes
	return Imbalance(0.0, 0.0);
}

} /* namespace problem */
