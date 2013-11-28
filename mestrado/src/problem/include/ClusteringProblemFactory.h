/*
 * ClusteringProblemFactory.h
 *
 *  Created on: May 30, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEMFACTORY_H_
#define CLUSTERINGPROBLEMFACTORY_H_

#include "CCProblem.h"
#include "RCCProblem.h"

namespace problem {

class ClusteringProblemFactory {
public:
	CCProblem ccProblem;
	RCCProblem rccProblem;

	ClusteringProblemFactory();
	virtual ~ClusteringProblemFactory();

	ClusteringProblem& build(int problemType) {
		if(problemType == ClusteringProblem::CC_PROBLEM) {
			return ccProblem;
		} else if(problemType == ClusteringProblem::RCC_PROBLEM) {
			return rccProblem;
		}
		cerr << "Unknown problem type to build!!! Returning CC problem instance.\n";
		return ccProblem;
	}
};

} /* namespace problem */
#endif /* CLUSTERINGPROBLEMFACTORY_H_ */
