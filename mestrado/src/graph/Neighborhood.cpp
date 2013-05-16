/*
 * Neighborhood.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/Neighborhood.h"
#include <limits>


namespace clusteringgraph {

NeighborhoodList::NeighborhoodList(Clustering *c, int n) :
		numberOfNodes(n), clusteringPtr(c) {

}

Clustering* NeighborhoodList::findLocalOptimum(SignedGraph *g, ClusteringProblem* problem) {
	Clustering* Cl = NULL;
	int bestValue = std::numeric_limits<int>::max();

	// for all C in N do
	for(unsigned int i = 0; i < size(); i++) {
		Clustering* C = new Clustering(at(i).get(), numberOfNodes);
		// if Q(C) > Q(Cl) then
		if(problem->objectiveFunction(g, C) < bestValue)
			Cl = C;
		// end if
		// N = N \ {c}
	}
	return Cl;
}


// TODO: Implementar de acordo com o especificado pelo Yuri
// para 1-opt e 2-opt
NeighborhoodList* NeighborhoodList::generateNeighborhood(int l) {
	return NULL;
}

} /* namespace clusteringgraph */
