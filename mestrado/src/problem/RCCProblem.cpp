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
					internalPositiveSum[k] += weight;
				} else { // negative edge
					internalNegativeSum[k] += fabs(weight);
				}
			} else {
				if(weight > 0) { // positive edge
					externalPositiveSum[k] += weight;
				} else { // negative edge
					externalNegativeSum[k] += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	for(int k = 0; k < nc; k++) {
		internalSum += min(internalPositiveSum[k], internalNegativeSum[k]);
		externalSum += min(externalPositiveSum[k], externalNegativeSum[k]);
	}
	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	return Imbalance(internalSum, externalSum);
}

// Calculates the delta of the objective function
Imbalance RCCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, const ClusterList& c,
		const unsigned long& k, const unsigned long& i) const {
	double negativeSum = 0, positiveSum = 0;
	// Not implemented.
	return Imbalance(0.0, 0.0);
}

} /* namespace problem */
