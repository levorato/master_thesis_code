/*
 * RCCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef RCCPROBLEM_H_
#define RCCPROBLEM_H_

#include "ClusteringProblem.h"
#include "../../graph/include/Imbalance.h"

using namespace clusteringgraph;

namespace problem {

class RCCProblem: public problem::ClusteringProblem {
public:
	RCCProblem();
	virtual ~RCCProblem();

	virtual Imbalance objectiveFunction(SignedGraph* g, Clustering *c) const;

	virtual int getType() const;
};

} /* namespace problem */
#endif /* RCCPROBLEM_H_ */
