/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#include "include/Grasp.h"
#include "include/VertexSet.h"
#include "../../graph/include/Neighborhood.h"
#include "../../problem/include/ClusteringProblem.h"
#include "../../problem/include/CCProblem.h"

#include <algorithm>
#include <boost/math/special_functions/round.hpp>

using namespace problem;

namespace resolution {
namespace grasp {

Grasp::Grasp() {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

Clustering* Grasp::executeGRASP(SignedGraph* g, int iter, float alpha, int l,
		ClusteringProblem* problem) {
	Clustering* CStar = NULL;
	std::cout << "Initializing GRASP procedure...\n";
	unsigned int ramdomSeed = 0;

	for (int i = 0; i < iter; i++) {
		// 1. Construct the clustering
		Clustering* Cc = constructClustering(g, alpha, ramdomSeed);
		// 2. Execute local search algorithm
		Clustering* Cl = localSearch(g, Cc, l, problem);
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		if(problem->objectiveFunction(g, Cl) < problem->objectiveFunction(g, CStar)) {
			CStar = Cl;
		}
	}
	return CStar;
}

Clustering* Grasp::constructClustering(SignedGraph* g, float alpha, unsigned int ramdomSeed) {
	Clustering *Cc = new Clustering(g->getN()); // Cc = empty
	VertexSet *lc = new VertexSet(g->getN()); // L(Cc) = V(G)
	std::cout << "GRASP construct clustering...\n";

	while(lc->size() > 0) { // lc != empty
		// 1. Compute L(Cc): this is done automatically by the VertexSet class (lc)

		// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		int i = lc->chooseRandomVertex(boost::math::iround(alpha * lc->size()));
		std::cout << "Random vertex is " << i << std::endl;

		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		// TODO: Colocar o vertice i em um cluster a parte ou em um cluster
		// existente, dependendo do valor calculado da funcao gain gc(i)
		int vertexList[1] = {i};
		Cc->addCluster(vertexList, 1);

		// 4. lc = lc - {i}
		// the removal of vertex i automatically recalculates the gain function
		lc->removeVertex(i);

		Cc->printClustering();
	}
	std::cout << "\nInitial clustering completed.\n";
	Cc->printClustering();
	return Cc;
}

Clustering* Grasp::localSearch(SignedGraph* g, Clustering* Cc, int l,
		ClusteringProblem* problem) {
	// k is the current neighborhood distance in the local search
	int k = 1;
	Clustering* cStar = Cc; // C* := Cc
	std::cout << "GRASP local search...\n";

	while(k <= l) {
		// N := Nl(C*)
		// apply a local search in cStar using the k-neighborhood
		NeighborhoodList neig(cStar, g->getN());
		ClusteringProblem* problem = new problem::CCProblem();
		neig.generateNeighborhood(k);
		Clustering* cl = neig.findLocalOptimum(g, problem);
		if(problem->objectiveFunction(g, cl) < problem->objectiveFunction(g, cStar)) {
			cStar = cl;
			k = 1;
		} else {
			k++;
		}
	}
	std::cout << "GRASP local search done.\n";
	return cStar;
}

} /* namespace grasp */
} /* namespace resolution */
