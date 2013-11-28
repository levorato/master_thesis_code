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
#include "../../graph/include/Graph.h"

using namespace clusteringgraph;

namespace problem {

class CCProblem: public problem::ClusteringProblem {
public:
	CCProblem();
	virtual ~CCProblem();

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
#endif /* CCPROBLEM_H_ */
