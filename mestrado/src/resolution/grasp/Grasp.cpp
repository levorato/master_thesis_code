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

ClusteringPtr Grasp::executeGRASP(SignedGraph* g, int iter, float alpha, int l,
		ClusteringProblem* problem) {
	std::cout << "Initializing GRASP procedure...\n";
	unsigned int ramdomSeed = 0;
	ClusteringPtr CStar(constructClustering(g, alpha, ramdomSeed));

	for (int i = 0; i < iter; i++) {
		cout << "GRASP iteration " << i << endl;
		// 1. Construct the clustering
		ClusteringPtr Cc(constructClustering(g, alpha, ramdomSeed));
		// 2. Execute local search algorithm
		ClusteringPtr Cl(localSearch(g, Cc.get(), l, problem));
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		if(problem->objectiveFunction(g, Cl.get()) < problem->objectiveFunction(g, CStar.get())) {
			cout << "A better solution was found." << endl;
			CStar = Cl;
		}
	}
	return CStar;
}

ClusteringPtr Grasp::constructClustering(SignedGraph* g, float alpha, unsigned int ramdomSeed) {
	ClusteringPtr Cc(new Clustering(g->getN())); // Cc = empty
	VertexSet lc(g->getN()); // L(Cc) = V(G)
	std::cout << "GRASP construct clustering...\n";
	// It is important to calculate the modularity matrix first (used by vertex sorting)
	g->calculateModularityMatrix();

	while(lc.size() > 0) { // lc != empty
		cout << "Vertex list size is " << lc.size() << endl;

		// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
		lc.sort(g, Cc.get());

		// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		int i = lc.chooseRandomVertex(boost::math::iround(alpha * lc.size()));
		std::cout << "Random vertex is " << i << std::endl;

		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		GainCalculation gainCalculation = Cc->gain(g, i);
		if(gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			int vertexList[1] = {i};
			Cc->addCluster(vertexList, 1);
		} else {
			// inserts i into existing cluster
			Cc->addNodeToCluster(i, gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the removal of vertex i automatically recalculates the gain function
		lc.removeVertex(i);

		Cc->printClustering();
	}
	std::cout << "\nInitial clustering completed.\n";
	Cc->printClustering();
	return Cc;
}

ClusteringPtr Grasp::localSearch(SignedGraph* g, Clustering* Cc, int l,
		ClusteringProblem* problem) {
	// k is the current neighborhood distance in the local search
	int k = 1, iteration = 0;
	ClusteringPtr cStar(Cc); // C* := Cc
	std::cout << "GRASP local search...\n";

	while(k <= l) {
		cout << "Local search iteration " << iteration << endl;
		// N := Nl(C*)
		// apply a local search in cStar using the k-neighborhood
		NeighborhoodList neig(cStar.get(), g->getN());
		cout << "Generating neighborhood of size l = " << k << endl;
		ClusteringPtr cl = neig.generateNeighborhood(k, g, problem);
		if(problem->objectiveFunction(g, cl.get()) < problem->objectiveFunction(g, cStar.get())) {
			cStar = cl;
			k = 1;
		} else {
			k++;
		}
		iteration++;
	}
	std::cout << "GRASP local search done.\n";
	return cStar;
}

} /* namespace grasp */
} /* namespace resolution */
