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
#include "../../../graph/include/Neighborhood.h"

using namespace clusteringgraph;
using namespace problem;

namespace resolution {
namespace grasp {

class Grasp: public ResolutionStrategy {
public:
	Grasp(GainFunction* f);
	virtual ~Grasp();

	/**
	 * Executes the GRASP algorithm. Returns the local optimum
	 * solution C(l) in the l-neighborhood of the current solution C.
	 * This GRASP algorithm consists of two phases: constructClustering
	 * and localSearch.
	 * @param g the graph to be used as the base
	 * @param iter maximum number of iterations
	 * @param alpha ramdom seed belonging to the interval (0, 1)
	 * @param l the size of the neighborhood
	 * @param problem the ClusteringProblem (objective function) to be used
	 * @param fileId string representing the identification of the input graph file
	 * @param timeLimit Maximum processing time in seconds
	 * @param myRank processor rank (when running with MPI)
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, const int& iter, const double& alpha, const int& l,
			const ClusteringProblem& problem, string& timestamp, string& fileId, 
			string& outputFolder, const long& timeLimit, const int& myRank);

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
	ClusteringPtr constructClustering(SignedGraph *g, const ClusteringProblem& problem,
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
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @param timeLimit Maximum processing time in seconds.
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr localSearch(SignedGraph *g, Clustering& Cc, const int &l,
			const ClusteringProblem& problem, NeighborhoodListGenerator &neig,
			const long& timeLimit, const int& myRank);

	/**
	 * TODO document this method
	 */
	void generateOutputFile(stringstream& fileContents, string& outputFolder, string& fileId, string& timestamp,
			const int &processNumber, const double& alpha, const int& l, const int& numberOfIterations);

	/**
	 * Time spent so far in the GRASP.
	 */
	double timeSpentSoFar;

	/**
	 * The gain function to be used in the construction phase.
	 */
	GainFunction* gainFunction;
};

} /* namespace grasp */
} /* namespace resolution */

#endif /* GRASP_H_ */
