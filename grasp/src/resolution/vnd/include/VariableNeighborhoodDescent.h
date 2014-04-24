/*
 * RCCLocalSearch.h
 *
 *  Created on: 28/11/2013
 *      Author: Mario Levorato
 */

#ifndef RCCLOCALSEARCH_H_
#define RCCLOCALSEARCH_H_

#include "../../include/ResolutionStrategy.h"
#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"
#include "../../../problem/include/ClusteringProblem.h"
#include "../../../graph/include/NeighborhoodSearch.h"

using namespace clusteringgraph;
using namespace problem;

namespace resolution {
namespace vnd {

class VariableNeighborhoodDescent : public ResolutionStrategy {
public:
	VariableNeighborhoodDescent(unsigned long seed);
	virtual ~VariableNeighborhoodDescent();

	/**
	 * Executes the local search algorithm. Repeatedly derives
	 * the local optimum solution C(l) in the l-neighborhood of
	 * the current solution C.
	 * This is the local search phase of the metaheuristic. It uses the
	 * Variable Neighborhood Descent (VND) algorithm.
	 * @param g the graph to be used as the base
	 * @param Cc the clustering given by the construct clustering
	 * phase of GRASP.
	 * @param l the size of the neighborhood
	 * @param graspIteration the number of the grasp iteration this LS belongs to
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @param timeLimit Maximum processing time in seconds.
	 * @return Clustering C(l), the local optinum solution
	 */
	Clustering localSearch(SignedGraph *g, Clustering& Cc, const int &l, const int& graspIteration,
			const bool& firstImprovementOnOneNeig, ClusteringProblem& problem, NeighborhoodSearch &neig,
			const long& timeLimit, const long& timeSpentSoFar, const int &numberOfSlaves,
			const int& myRank, const int& numberOfSearchSlaves);

	unsigned long getNumberOfTestedCombinations();

private:
	/**
	 * Computes the best known search result at each time interval.
	 */
	void measureTimeResults(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	void notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	/**
	 * Time spent so far in the Local Search in seconds.
	 */
	double timeSpentInSearch;

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search.
	 */
	stringstream timeResults;

	/**
	 * Random seed.
	 */
	unsigned long randomSeed;

	/**
	 * Time measure interval for search results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;
	/**
	 * Total number of tested combinations in local search.
	 */
	unsigned long numberOfTestedCombinations;
};

} /* namespace vnd */
} /* namespace resolution */
#endif /* RCCLOCALSEARCH_H_ */
