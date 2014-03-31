/*
 * Grasp.h
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#ifndef GRASP_H_
#define GRASP_H_

#include "GainFunction.h"
#include "../../include/ResolutionStrategy.h"
#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"
#include "../../../problem/include/ClusteringProblem.h"
#include "../../../graph/include/NeighborhoodSearch.h"

#include <iostream>

using namespace clusteringgraph;
using namespace problem;

namespace resolution {
namespace grasp {

class Grasp: public ResolutionStrategy {
public:
	Grasp(GainFunction& f, unsigned long seed);
	virtual ~Grasp();

	/**
	 * Executes the GRASP algorithm. Returns the local optimum
	 * solution C(l) in the l-neighborhood of the current solution C.
	 * This GRASP algorithm consists of two phases: constructClustering
	 * and localSearch.
	 * @param g the graph to be used as the base
	 * @param iter maximum number of iterations
	 * @param alpha randomness factor belonging to the interval (0, 1); if alpha < 0, the
	 * constructClustering method will always choose the first vertex in the gainFunction list,
	 * that is, the one that minimizes the objective (VOTE algorithm).
	 * @param l the size of the neighborhood
	 * @param problem the ClusteringProblem (objective function) to be used
	 * @param fileId string representing the identification of the input graph file
	 * @param timeLimit Maximum processing time in seconds
	 * @param numberOfSlaves number of slaves used for parallel GRASP processing
	 * @param myRank processor rank (when running with MPI)
	 * @param numberOfSearchSlaves number of slaves used for parallel VNS processing
	 */
	Clustering executeGRASP(SignedGraph *g, const int& iter, const double& alpha, const int& l,
			const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const long& timeLimit,
			const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);

	unsigned long getNumberOfTestedCombinations();

protected:
	/**
	 * Constructs a clustering in a greedy ramdomized fashion,
	 * starting from the empty set.
	 * This is the first phase of the GRASP algorithm.
	 * @param g graph to be used as the base
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @param alpha ramdom seed belonging to the interval (0, 1)
	 * @param ramdomSeed seed to be used in ramdom number generation
	 * @return Clustering C(c)
	 */
	Clustering constructClustering(SignedGraph *g, ClusteringProblem& problem,
			double alpha, int myRank);

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
	Clustering localSearch(SignedGraph *g, Clustering& Cc, const int &l, const int& graspIteration,
			const bool& firstImprovementOnOneNeig, ClusteringProblem& problem, NeighborhoodSearch &neig,
			const long& timeLimit, const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);

	/**
	 * Generates CSV output file for GRASP local Search.
	 */
	void generateOutputFile(ClusteringProblem& problem, stringstream& fileContents, const string& rootFolder,
			const string& fileId, const string& timestamp, const int &processNumber, const string& fileSuffix,
			const double& alpha, const int& l, const int& numberOfIterations);

	/**
	 * Computes the best known GRASP result at each time interval.
	 */
	void measureTimeResults(const double& timeSpentOnLocalSearch, const int& graspIteration);

	void notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& graspIteration) {
		Imbalance i1 = CStar.getImbalance();
		Imbalance i2 = CBest.getImbalance();
		if(i1 < i2) {
			CBest = CStar;
			measureTimeResults(timeSpentOnLocalSearch, graspIteration);
		}
	}

	/**
	 * Time spent so far in the GRASP in seconds.
	 */
	double timeSpentInGRASP;

	/**
	 * The gain function to be used in the construction phase.
	 */
	GainFunction& gainFunction;

	/**
	 * Random seed.
	 */
	unsigned long randomSeed;

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search phase of GRASP.
	 */
	stringstream timeResults;

	/**
	 * Time measure interval for GRASP results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;
	/**
	 * Total number of tested combinations in GRASP local search.
	 */
	unsigned long numberOfTestedCombinations;
	/** The best clustering found in all GRASP iterations. */
	Clustering CBest, CBefore;
};

} /* namespace grasp */
} /* namespace resolution */

#endif /* GRASP_H_ */
