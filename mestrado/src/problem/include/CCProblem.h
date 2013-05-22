/*
 * CCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef CCPROBLEM_H_
#define CCPROBLEM_H_

#include "ClusteringProblem.h"

using namespace clusteringgraph;

namespace problem {

class CCProblem: public problem::ClusteringProblem {
public:
	CCProblem();
	virtual ~CCProblem();

	virtual float objectiveFunction(SignedGraph* g, Clustering* c) const;
};

} /* namespace problem */
#endif /* CCPROBLEM_H_ */
