/*
 * CUDAVariableNeighborhoodDescent.h
 *
 *  Created on: 5/6/2014
 *      Author: Mario Levorato
 */

#ifndef CUDAVND_H_
#define CUDAVND_H_

#include "../../include/ResolutionStrategy.h"
#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"
#include "graph/include/NeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
// #include "../../include/LocalSearch.h"
#include "VariableNeighborhoodDescent.h"


using namespace clusteringgraph;
using namespace problem;

namespace resolution {
namespace vnd {

class CUDAVariableNeighborhoodDescent : public VariableNeighborhoodDescent {
public:
	/**
	 * @param seed random seed
	 * @param lsize the size of the neighborhood
	 * @param firstImprovement1Opt true if first-improvement on 1-opt is enabled
	 * @param tlimit time limit in seconds
	 */
	CUDAVariableNeighborhoodDescent(NeighborhoodSearch &neighborhoodSearch, unsigned long seed,
			const int &lsize, const bool& firstImprovement1Opt, const long &tlimit);
	virtual ~CUDAVariableNeighborhoodDescent();

	/**
	 * Executes the local search algorithm. Repeatedly derives
	 * the local optimum solution C(l) in the l-neighborhood of
	 * the current solution C.
	 * This is the local search phase of the metaheuristic. It uses the
	 * Variable Neighborhood Descent (VND) algorithm.
	 * @param g the graph to be used as the base
	 * @param Cc the clustering given by the construct clustering
	 * phase of GRASP.
	 * @param graspIteration the number of the grasp iteration this LS belongs to
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @return Clustering C(l), the local optimum solution
	 */
	virtual Clustering localSearch(SignedGraph *g, Clustering& Cc, const int& graspIteration,
			ClusteringProblem& problem, const long& timeSpentSoFar, const int& myRank);

private:
	/**
	 * Computes the best known search result at each time interval.
	 */
	void measureTimeResults(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	void notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search.
	 */
	stringstream timeResults;

	/**
	 * Time measure interval for search results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;

};

} /* namespace vnd */
} /* namespace resolution */
#endif /* CUDAVND_H_ */
