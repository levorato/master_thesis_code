/*
 * ClusteringProblem.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEM_VALIDATION_H_
#define CLUSTERINGPROBLEM_VALIDATION_H_

#include "./Graph.h"
#include "graph/include/Imbalance.h"

namespace clusteringgraph{
namespace validation {
class Clustering;
}
}

namespace problem {
namespace validation {

using namespace clusteringgraph::validation;

class ClusteringProblem {
public:
	static const int CC_PROBLEM = 0, RCC_PROBLEM = 1;

	ClusteringProblem();
	virtual ~ClusteringProblem();

	virtual clusteringgraph::Imbalance objectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c) = 0;

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	virtual clusteringgraph::Imbalance calculateDeltaPlusObjectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c,
			const unsigned long& k, const unsigned long& i) = 0;

	/**
	 * Calculates the delta of the objective function caused by the
	 * removal of node i in cluster k.
	 */
	virtual clusteringgraph::Imbalance calculateDeltaMinusObjectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c,
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

}
} /* namespace problem */
#endif /* CLUSTERINGPROBLEM_H_ */
