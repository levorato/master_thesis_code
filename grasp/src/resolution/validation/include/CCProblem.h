/*
 * CCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef CCPROBLEM_VALIDATION_H_
#define CCPROBLEM_VALIDATION_H_

#include "./ClusteringProblem.h"
#include "graph/include/Imbalance.h"
#include "./Graph.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace clusteringgraph::validation;
using namespace boost::numeric::ublas;

namespace problem {
namespace validation {

class CCProblem: public problem::validation::ClusteringProblem {
public:
	CCProblem();
	virtual ~CCProblem();

	virtual clusteringgraph::Imbalance objectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c);

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	virtual clusteringgraph::Imbalance calculateDeltaPlusObjectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c,
			const unsigned long& k, const unsigned long& i);

	/**
	 * Calculates the delta of the objective function caused by the
	 * removal of node i in cluster k.
	 */
	virtual clusteringgraph::Imbalance calculateDeltaMinusObjectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c,
			const unsigned long& k, const unsigned long& i);

	clusteringgraph::Imbalance calculateDeltaObjectiveFunction(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c,
			const unsigned long& k, const unsigned long& i);

	virtual int getType() const;

	virtual string getName();
	
	string analyzeImbalance(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c);

	matrix<double> calculateClusterToClusterImbalanceMatrix(clusteringgraph::validation::SignedGraph& g, clusteringgraph::validation::Clustering& c);
};

}
} /* namespace problem */
#endif /* CCPROBLEM_H_ */
