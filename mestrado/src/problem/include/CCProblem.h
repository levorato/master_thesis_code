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

	virtual Imbalance objectiveFunction(SignedGraph& g, Clustering& c);

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	virtual Imbalance calculateDeltaPlusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i);

	/**
	 * Calculates the delta of the objective function caused by the
	 * removal of node i in cluster k.
	 */
	virtual Imbalance calculateDeltaMinusObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i);

	Imbalance calculateDeltaObjectiveFunction(SignedGraph& g, Clustering& c,
			const unsigned long& k, const unsigned long& i);

	virtual int getType() const;
};

} /* namespace problem */
#endif /* CCPROBLEM_H_ */
