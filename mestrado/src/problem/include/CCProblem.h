/*
 * CCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef CCPROBLEM_H_
#define CCPROBLEM_H_

#include "ClusteringProblem.h"
#include "../../graph/include/Imbalance.h"

using namespace clusteringgraph;

namespace problem {

class CCProblem: public problem::ClusteringProblem {
public:
	CCProblem();
	virtual ~CCProblem();

	virtual Imbalance objectiveFunction(SignedGraph* g, Clustering* c) const;

	virtual int getType() const;
};

} /* namespace problem */
#endif /* CCPROBLEM_H_ */
