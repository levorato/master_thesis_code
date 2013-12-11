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
#include "../../graph/include/Clustering.h"

using namespace clusteringgraph;

namespace problem {

class RCCProblem: public problem::ClusteringProblem {
public:
	RCCProblem();
	virtual ~RCCProblem();

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


	virtual int getType() const;
};

} /* namespace problem */
#endif /* RCCPROBLEM_H_ */
