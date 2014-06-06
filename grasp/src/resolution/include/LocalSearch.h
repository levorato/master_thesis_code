/*
 * LocalSearch.h
 *
 *  Created on: Jun 5, 2014
 *      Author: mlevorato
 */

#ifndef LOCALSEARCH_H_
#define LOCALSEARCH_H_

#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"

namespace clusteringgraph {

class LocalSearch {
public:
	LocalSearch(unsigned long seed, const int &lsize, const bool& firstImprovement1Opt, const long &tlimit,
			const double &timeOnLS, const unsigned long &numcomb, const double &tsearch) : randomSeed(seed), l(lsize),
			firstImprovementOnOneNeig(firstImprovement1Opt),
			timeLimit(tlimit), timeSpentOnLocalSearch(timeOnLS), numberOfTestedCombinations(numcomb), timeSpentInSearch(tsearch)
	{

	}
	virtual ~LocalSearch();

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
			ClusteringProblem& problem, const long& timeSpentSoFar, const int& myRank) = 0;

	unsigned int getNeighborhoodSize() {
		return l;
	}

	double getTimeSpentOnLocalSearch() {
		return timeSpentOnLocalSearch;
	}

	unsigned long getTimeLimit() {
		return timeLimit;
	}

	bool isFirstImprovementOnOneNeig() {
		return firstImprovementOnOneNeig;
	}

	unsigned long getNumberOfTestedCombinations() {
		return numberOfTestedCombinations;
	}

	unsigned long getRandomSeed() {
		return randomSeed;
	}

protected:
	/**
	 * Neighborhood size of local search.
	 */
	int l;

	/**
	 * Controls if first-improvement is enabled on 1-opt neighborhood.
	 */
	bool firstImprovementOnOneNeig;

	/**
	 * Time limit of local search execution in seconds.
	 */
	long timeLimit;
	/**
	 * Time spent in the local search.
	 */
	double timeSpentOnLocalSearch;

	/**
	 * Total number of tested combinations in local search.
	 */
	unsigned long numberOfTestedCombinations;
	/**
	 * Time spent so far in the Local Search in seconds.
	 */
	double timeSpentInSearch;

	/**
	 * Random seed.
	 */
	unsigned long randomSeed;

};

} /* namespace clusteringgraph */

#endif /* LOCALSEARCH_H_ */
