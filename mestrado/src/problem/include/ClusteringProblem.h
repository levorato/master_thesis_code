/*
 * ClusteringProblem.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEM_H_
#define CLUSTERINGPROBLEM_H_

#include "../../graph/include/Clustering.h"
#include "../../graph/include/Imbalance.h"

using namespace clusteringgraph;

namespace problem {

class ClusteringProblem {
public:
	static const int CC_PROBLEM = 0, RCC_PROBLEM = 1;

	ClusteringProblem();
	virtual ~ClusteringProblem();

	virtual Imbalance objectiveFunction(SignedGraph* g, Clustering* c) const = 0;

	virtual int getType() const = 0;

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
