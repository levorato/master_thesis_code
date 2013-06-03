/*
 * Grasp.h
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#ifndef GRASP_H_
#define GRASP_H_

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
	Grasp();
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
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, const int& iter, const float& alpha, const int& l,
			const ClusteringProblem& problem, string& timestamp, string& fileId, const int& myRank);

private:
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
			float alpha, unsigned int ramdomSeed);

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
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr localSearch(SignedGraph *g, Clustering& Cc, const int &l,
			const ClusteringProblem& problem, NeighborhoodListGenerator &neig);

	/**
	 * TODO document this method
	 */
	void generateOutputFile(stringstream& fileContents, string& fileId, string& timestamp,
			const int &processNumber, const float& alpha, const int& l, const int& numberOfIterations);
};

} /* namespace grasp */
} /* namespace resolution */

#endif /* GRASP_H_ */
