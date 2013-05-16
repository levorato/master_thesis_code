/*
 * RCCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef RCCPROBLEM_H_
#define RCCPROBLEM_H_

#include "ClusteringProblem.h"

namespace problem {

class RCCProblem: public problem::ClusteringProblem {
public:
	RCCProblem();
	virtual ~RCCProblem();

	int objectiveFunction(SignedGraph* g, Clustering *c);
};

} /* namespace problem */
#endif /* RCCPROBLEM_H_ */
