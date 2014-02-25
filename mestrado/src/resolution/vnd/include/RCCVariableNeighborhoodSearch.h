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

class RCCVariableNeighborhoodSearch : public ResolutionStrategy {
public:
	RCCVariableNeighborhoodSearch(unsigned long seed);
	virtual ~RCCVariableNeighborhoodSearch();

	/**
	 * Executes the local search algorithm for RCC. Repeatedly derives
	 * the local optimum solution C(l) in the l-neighborhood of
	 * the current solution C assuming that number of clusters <= k.
	 * It is based on Variable Neighborhood Descent (VND) algorithm.
	 * @param g the graph to be used as base
	 * @param Cc the best clustering provided by the CC GRASP algorithm.
	 * @param l the size of the neighborhood
	 * @param k the maximum number of clusters in the solution
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @param timeLimit Maximum processing time in seconds.
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr executeSearch(SignedGraph *g, Clustering& Cc, const int &l,
			const long unsigned int& k, const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
			NeighborhoodSearch &neig, string& executionId, string& fileId, string& outputFolder,
			const long& timeLimit, const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);

	unsigned long getNumberOfTestedCombinations();

private:
	/**
	 * Generates CSV output file for RCC Search.
	 */
	void generateOutputFile(stringstream& fileContents, const string& rootFolder, const string& fileId, const string& timestamp,
			const int &processNumber, const string& fileSuffix, const unsigned int& k, const int& l);

	/**
	 * Computes the best known RCC Search result at each time interval.
	 */
	void measureTimeResults(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	void notifyNewValue(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& iteration);

	/**
	 * Time spent so far in the RCC Search in seconds.
	 */
	double timeSpentInSearch;

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search of RCC.
	 */
	stringstream timeResults;

	/**
	 * Random seed.
	 */
	unsigned long randomSeed;

	/**
	 * Time measure interval for RCC results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;
	/**
	 * Total number of tested combinations in RCC local search.
	 */
	unsigned long numberOfTestedCombinations;
};

} /* namespace vnd */
} /* namespace resolution */
#endif /* RCCLOCALSEARCH_H_ */
