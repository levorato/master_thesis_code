/*
 * ClusteringProblem.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEM_H_
#define CLUSTERINGPROBLEM_H_

#include "../../graph/include/Clustering.h"

using namespace clusteringgraph;

namespace problem {

class ClusteringProblem {
public:
	ClusteringProblem();
	virtual ~ClusteringProblem();

	virtual int objectiveFunction(Clustering *c);

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
