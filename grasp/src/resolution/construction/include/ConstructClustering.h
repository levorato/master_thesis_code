/*
 * ConstructClustering.h
 *
 *  Created on: Apr 23, 2014
 *      Author: mlevorato
 */

#ifndef CONSTRUCTCLUSTERING_H_
#define CONSTRUCTCLUSTERING_H_

#include "../../construction/include/GainFunction.h"
#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"

namespace resolution {
namespace construction {

class ConstructClustering {
public:
	ConstructClustering(GainFunction* f, const unsigned long& seed, const double& alpha);
	virtual ~ConstructClustering();

	/**
	 * Constructs a clustering in a greedy randomized fashion,
	 * starting from the empty set.
	 * This is the first phase of the metaheuristic.
	 * @param g graph to be used as the base
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @return Clustering Cc
	 */
	Clustering constructClustering(SignedGraph *g, ClusteringProblem& problem,
			const int& myRank);

	double getAlpha() {
		return _alpha;
	}

	int getGainFunctionType() {
		return gainFunction->getType();
	}

private:
	GainFunction* gainFunction;
	unsigned long randomSeed;
	/**
	 * alpha parameter belonging to the interval [0, 1]
	 * if alpha < 0, the constructClustering method will always choose the first vertex
	 * in the gainFunction list, that is, the one that minimizes the objective (VOTE algorithm).
	 */
	double _alpha;
};

} /* namespace construction */
} /* namespace resolution */

#endif /* CONSTRUCTCLUSTERING_H_ */
