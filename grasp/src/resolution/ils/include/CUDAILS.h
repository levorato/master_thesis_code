/*
 * CUDAILS.h
 *
 *  Created on: Mar 2, 2015
 *      Author: mlevorato
 */

#ifndef CUDAILS_H_
#define CUDAILS_H_

#include "ILS.h"

namespace resolution {
namespace ils {

using namespace problem;
using namespace resolution::construction;

class CUDAILS : public ILS {
public:
	CUDAILS();
	virtual ~CUDAILS();

	/**
	 * Executes the ILS algorithm. Returns the local optimum
	 * solution C(l) in the l-neighborhood of the current solution C.
	 * This ILS algorithm consists of two phases: constructClustering
	 * and localSearch.
	 * @param g the graph to be used as the base
	 * @param iterMax maximum number of iterations of multistart ILS
	 * @param iterMaxILS maximum number of iterations of internal ILS loop
	 * @param perturbationLevelMax maximum perturbation level for ILS
	 * @param problem the ClusteringProblem (objective function) to be used
	 * @param executionInfo auxiliary data about execution
	 */
	Clustering executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			SignedGraph *g, const int& iterMax, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem,	ExecutionInfo& info);

};

} /* namespace ils */
} /* namespace resolution */

#endif /* CUDAILS_H_ */
