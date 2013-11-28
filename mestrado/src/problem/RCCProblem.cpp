/*
 * RCCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/RCCProblem.h"

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
	cerr << "Unimplemented function called!!!n";
	return Imbalance(0, 0);
}

// Calculates the delta of the objective function
Imbalance RCCProblem::calculateDeltaObjectiveFunction(SignedGraph& g, const ClusterList& c,
		const unsigned long& k, const unsigned long& i) const {
	cerr << "Unimplemented function called!!!n";
	return Imbalance(0, 0);
}

} /* namespace problem */
