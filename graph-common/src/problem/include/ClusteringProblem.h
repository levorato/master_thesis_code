/*
 * ClusteringProblem.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEM_H_
#define CLUSTERINGPROBLEM_H_

#include "../../graph/include/Graph.h"
#include "../../graph/include/ParallelBGLSignedGraph.h"
#include "../../graph/include/Imbalance.h"

using namespace clusteringgraph;

namespace clusteringgraph{
class Clustering;
}

namespace problem {

class ClusteringProblem {
public:
	static const int CC_PROBLEM = 0, RCC_PROBLEM = 1;

	ClusteringProblem();
	virtual ~ClusteringProblem();

	virtual Imbalance objectiveFunction(SignedGraph& g, Clustering& c) = 0;

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	virtual Imbalance calculateDeltaPlusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i) = 0;

	/**
	 * Calculates the delta of the objective function caused by the
	 * removal of node i in cluster k.
	 */
	virtual Imbalance calculateDeltaMinusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i) = 0;

	virtual int getType() const = 0;

	virtual string getName() = 0;

	int getNumberOfIterations() const {
		return numberOfIterations;
	}

	void setNumberOfIterations(int numberOfIterations) {
		this->numberOfIterations = numberOfIterations;
	}


private:
	int numberOfIterations;

};

} /* namespace problem */
#endif /* CLUSTERINGPROBLEM_H_ */
