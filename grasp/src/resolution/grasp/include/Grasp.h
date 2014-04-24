/*
 * Grasp.h
 *
 *  Created on: 30/04/2013
 *      Author: Mario Levorato
 */

#ifndef GRASP_H_
#define GRASP_H_

#include "../../construction/include/ConstructClustering.h"
#include "../../vnd/include/VariableNeighborhoodDescent.h"
#include "../../include/ResolutionStrategy.h"
#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"
#include "../../../problem/include/ClusteringProblem.h"
#include "../../../graph/include/NeighborhoodSearch.h"

#include <iostream>

using namespace clusteringgraph;
using namespace problem;
using namespace resolution::construction;
using namespace resolution::vnd;

namespace resolution {
namespace grasp {

class Grasp: public ResolutionStrategy {
public:
	Grasp();
	virtual ~Grasp();

	/**
	 * Executes the GRASP algorithm. Returns the local optimum
	 * solution C(l) in the l-neighborhood of the current solution C.
	 * This GRASP algorithm consists of two phases: constructClustering
	 * and localSearch.
	 * @param ConstructClustering constructClustering - the construct clustering phase
	 * @param VariableNeighborhoodDescent vnd the local search phase (VND)
	 * @param g the graph to be used as the base
	 * @param iter maximum number of iterations
	 * constructClustering method will always choose the first vertex in the gainFunction list,
	 * that is, the one that minimizes the objective (VOTE algorithm).
	 * @param problem the ClusteringProblem (objective function) to be used
	 * @param fileId string representing the identification of the input graph file
	 * @param myRank processor rank (when running with MPI)
	 */
	Clustering executeGRASP(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
			SignedGraph *g, const int& iter, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const int& myRank);

	unsigned long getNumberOfTestedCombinations();

protected:

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
