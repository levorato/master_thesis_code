/*
 * ILS.h
 *
 *  Created on: 31/03/2014
 *      Author: Mario Levorato
 */

#ifndef ILS_H_
#define ILS_H_

#include "../../construction/include/GainFunction.h"
#include "../../include/ResolutionStrategy.h"
#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"
#include "../../../problem/include/ClusteringProblem.h"
#include "../../../graph/include/NeighborhoodSearch.h"
#include "../../construction/include/ConstructClustering.h"
#include "../../vnd/include/VariableNeighborhoodDescent.h"

#include <iostream>

using namespace clusteringgraph;
using namespace problem;
using namespace resolution::construction;
using namespace resolution::vnd;

namespace resolution {
namespace ils {

class ILS: public ResolutionStrategy {
public:
	ILS();
	virtual ~ILS();

	/**
	 * Executes the ILS algorithm. Returns the local optimum
	 * solution C(l) in the l-neighborhood of the current solution C.
	 * This ILS algorithm consists of two phases: constructClustering
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
	 * @param numberOfSlaves number of slaves used for parallel ILS processing
	 * @param myRank processor rank (when running with MPI)
	 * @param numberOfSearchSlaves number of slaves used for parallel VND processing
	 */
	Clustering executeILS(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
			SignedGraph *g, const int& iter, const int& l,
			const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
			string& executionId, string& fileId, string& outputFolder, const long& timeLimit,
			const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves);

	unsigned long getNumberOfTestedCombinations();

protected:

	/**
	 * Generates CSV output file for ILS local Search.
	 */
	void generateOutputFile(ClusteringProblem& problem, stringstream& fileContents, const string& rootFolder,
			const string& fileId, const string& timestamp, const int &processNumber, const string& fileSuffix,
			const double& alpha, const int& l, const int& numberOfIterations);

	/**
	 * Computes the best known ILS result at each time interval.
	 */
	void measureTimeResults(const double& timeSpentOnLocalSearch, const int& ILSIteration);

	void notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& ILSIteration) {
		Imbalance i1 = CStar.getImbalance();
		Imbalance i2 = CBest.getImbalance();
		if(i1 < i2) {
			CBest = CStar;
			measureTimeResults(timeSpentOnLocalSearch, ILSIteration);
		}
	}

	/**
	 * Time spent so far in the ILS in seconds.
	 */
	double timeSpentInILS;

	/**
	 * Stringstream containing the best result found at each moment of time.
	 * Results are collected in local search phase of ILS.
	 */
	stringstream timeResults;

	/**
	 * Time measure interval for ILS results in seconds.
	 */
	static const double timeMeasureInterval = 10.0;
	double timeSum;
	/**
	 * Total number of tested combinations in ILS local search.
	 */
	unsigned long numberOfTestedCombinations;
	/** The best clustering found in all ILS iterations. */
	Clustering CBest, CBefore;
};

} /* namespace ils */
} /* namespace resolution */

#endif /* ILS_H_ */
