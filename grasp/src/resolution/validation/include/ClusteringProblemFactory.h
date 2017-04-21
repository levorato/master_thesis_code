/*
 * ClusteringProblemFactory.h
 *
 *  Created on: May 30, 2013
 *      Author: mario
 */

#ifndef CLUSTERINGPROBLEMFACTORY_VALIDATION_H_
#define CLUSTERINGPROBLEMFACTORY_VALIDATION_H_

#include "./CCProblem.h"

namespace problem {
namespace validation {

class ClusteringProblemFactory {
public:
	problem::validation::CCProblem ccProblem;

	ClusteringProblemFactory();
	virtual ~ClusteringProblemFactory();

	ClusteringProblem& build(int problemType, long k = 0) {
		if(problemType == ClusteringProblem::CC_PROBLEM)
			return ccProblem;
		cerr << "Unknown problem type to build!!! Returning CC problem instance.\n";
		return ccProblem;
	}
};

}
} /* namespace problem */
#endif /* CLUSTERINGPROBLEMFACTORY_H_ */
