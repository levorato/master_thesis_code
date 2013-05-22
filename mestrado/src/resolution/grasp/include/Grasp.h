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
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, int iter, float alpha, int l,
			ClusteringProblem& problem);

private:
	/**
	 * Constructs a clustering in a greedy ramdomized fashion,
	 * starting from the empty set.
	 * This is the first phase of the GRASP algorithm.
	 * @param g graph to be used as the base
	 * @param alpha ramdom seed belonging to the interval (0, 1)
	 * @return Clustering C(c)
	 */
	ClusteringPtr constructClustering(SignedGraph *g, float alpha,
			unsigned int ramdomSeed);

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
	 * @return Clustering C(l), the local optinum solution
	 */
	ClusteringPtr localSearch(SignedGraph *g, Clustering& Cc, int l,
			ClusteringProblem& problem);

};

} /* namespace grasp */
} /* namespace resolution */

#endif /* GRASP_H_ */
