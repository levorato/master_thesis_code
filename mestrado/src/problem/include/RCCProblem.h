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
#include "../../graph/include/Graph.h"

using namespace clusteringgraph;

namespace problem {

class RCCProblem: public problem::ClusteringProblem {
public:
	RCCProblem();
	virtual ~RCCProblem();

	virtual Imbalance objectiveFunction(SignedGraph& g, const ClusterList& c) const;

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	virtual Imbalance calculateDeltaObjectiveFunction(SignedGraph& g, const ClusterList& c,
			const unsigned long& k, const unsigned long& i) const;

	virtual int getType() const;
};

} /* namespace problem */
#endif /* RCCPROBLEM_H_ */
