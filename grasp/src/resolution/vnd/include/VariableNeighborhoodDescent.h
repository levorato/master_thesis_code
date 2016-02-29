/*
 * RCCLocalSearch.h
 *
 *  Created on: 28/11/2013
 *      Author: Mario Levorato
 */

#ifndef RCCLOCALSEARCH_H_
#define RCCLOCALSEARCH_H_

#include "../../include/ResolutionStrategy.h"
#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"
#include "graph/include/NeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "graph/include/ParallelNeighborhoodSearch.h"


using namespace clusteringgraph;
using namespace problem;

namespace resolution {
namespace vnd {

class VariableNeighborhoodDescent : public ResolutionStrategy {
public:

	VariableNeighborhoodDescent(const VariableNeighborhoodDescent &vnd);
	/**
	 * @param neighborhoodSearch neighborhood search algorithm (sequential, parallel)
	 * @param seed random seed
	 * @param lsize the size of the neighborhood
	 * @param firstImprovement1Opt true if first-improvement on 1-opt is enabled
	 * @param tlimit time limit in seconds
	 */
	VariableNeighborhoodDescent(NeighborhoodSearch &neigborhoodSearch, unsigned long seed,
			const int &lsize, const bool& firstImprovement1Opt, const long &tlimit);
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
	 * @param graspIteration the number of the grasp iteration this LS belongs to
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @return Clustering C(l), the local optimum solution
	 */
	virtual Clustering localSearch(SignedGraph *g, Clustering& Cc, const int& graspIteration,
			ClusteringProblem& problem, const long& timeSpentSoFar, const int& myRank);

	unsigned long getNumberOfTestedCombinations();

	unsigned int getNeighborhoodSize() {
		return l;
	}

	unsigned long getTimeLimit() {
		return timeLimit;
	}

	bool isFirstImprovementOnOneNeig() {
		return firstImprovementOnOneNeig;
	}

	unsigned long getRandomSeed() {
		return randomSeed;
	}

	double getTimeSpentOnLocalSearch() {
		return timeSpentOnLocalSearch;
	}

	void setTimeLimit(double time) {
		this->timeLimit = time;
	}

	/**
	 * Resets internal auxiliary matrices (used between Metaheuristic iterations)
	 * and updates the best clustering based on the constructive phase result.
	 */
	void reset(Clustering coClustering) {
		this->_neighborhoodSearch.reset(coClustering);
	}

protected:
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

	/**
	 * Associated NeigborhoodSearch class.
	 */
	NeighborhoodSearch& _neighborhoodSearch;

	/**
	 * Controls if first-improvement is enabled on 1-opt neighborhood.
	 */
	bool firstImprovementOnOneNeig;

	/**
	 * Neighborhood size of local search.
	 */
	unsigned int l;

	/**
	 * Time limit of local search execution in seconds.
	 */
	long timeLimit;

	/**
	 * Time spent on local search.
	 */
	double timeSpentOnLocalSearch;
};

} /* namespace vnd */
} /* namespace resolution */
#endif /* RCCLOCALSEARCH_H_ */
