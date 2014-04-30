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

class EdgeContribution {
public:
	EdgeContribution() : i(0), j(0), value(0) {

	}

	EdgeContribution(int _i, int _j, double _value) : i(_i), j(_j), value(_value) {

	}
	virtual ~EdgeContribution() {

	}

	int i, j;
	double value;
};

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

	string analyzeImbalance(SignedGraph& g, Clustering& c);

	virtual int getType() const;

	virtual string getName();

	void setK(long _k) {
		k = _k;
	}

	long getK() {
		return k;
	}

private:
	unsigned long k;
	static const int POSITIVE_EDGE = -1;
	static const int NEGATIVE_EDGE = 1;

	list<EdgeContribution> computeEdges(SignedGraph& g, Clustering& c, int c1, int c2, int edgeType);
};

} /* namespace problem */
#endif /* RCCPROBLEM_H_ */
