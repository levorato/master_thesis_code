/*
 * RCCLocalSearch.h
 *
 *  Created on: 28/11/2013
 *      Author: czt0
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

class RCCLocalSearch : public ResolutionStrategy {
public:
	RCCLocalSearch();
	virtual ~RCCLocalSearch();

	/**
	 * Executes the local search algorithm. Repeatedly derives
	 * the local optimum solution C(l) in the l-neighborhood of
	 * the current solution C.
	 * This is the second phase of the GRASP algorithm. It uses the
	 * Variable Neighborhood Descent (VND) strategy.
	 * @param g the graph to be used as the base
	 * @param Cc the clustering given by the construct clustering
	 * phase of GRASP.
	 * @param l the size of the neighborhood
	 * @param graspIteration the number of the grasp iteration this LS belongs to
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @param timeLimit Maximum processing time in seconds.
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr localSearch(SignedGraph *g, Clustering& Cc, const int &l, const int& graspIteration,
			const bool& firstImprovementOnOneNeig, const ClusteringProblem& problem, NeighborhoodSearch &neig,
			const long& timeLimit, const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);

	long getNumberOfTestedCombinations();

private:
	/**
	 * TODO document this method
	 */
	void generateOutputFile(stringstream& fileContents, const string& rootFolder, const string& fileId, const string& timestamp,
			const int &processNumber, const string& fileSuffix, const double& alpha, const int& l, const int& numberOfIterations);

	/**
	 * Computes the best known GRASP result at each time interval.
	 */
	void measureTimeResults(const double& timeSpentOnLocalSearch, const int& graspIteration);

	void notifyNewValue(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& graspIteration);

	/**
	 * Time spent so far in the GRASP in seconds.
	 */
	double timeSpentInSearch;

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search phase of GRASP.
	 */
	stringstream timeResults;

	/**
	 * Random seed.
	 */
	unsigned long randomSeed;

	/**
	 * Time measure interval for GRASP results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;
	/**
	 * Total number of tested combinations in GRASP local search.
	 */
	long numberOfTestedCombinations;
	/** The best clustering found in all iterations. */
	ClusteringPtr CBest, CBefore;
};

} /* namespace vnd */
} /* namespace resolution */
#endif /* RCCLOCALSEARCH_H_ */
