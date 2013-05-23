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
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace problem;
using namespace boost;

namespace resolution {
namespace grasp {

Grasp::Grasp() {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr Grasp::executeGRASP(SignedGraph *g, int iter, float alpha, int l,
		ClusteringProblem& problem) {
	std::cout << "Initializing GRASP procedure...\n";
	unsigned int ramdomSeed = 0;
	ClusteringPtr CStar = constructClustering(g, alpha, ramdomSeed);

	for (int i = 0; i < iter; i++) {
		cout << "GRASP iteration " << i << endl;
		// 1. Construct the clustering
		ClusteringPtr Cc = constructClustering(g, alpha, ramdomSeed);
		// 2. Execute local search algorithm
		ClusteringPtr Cl = localSearch(g, *Cc, l, problem);
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		float newValue = problem.objectiveFunction(g, Cl.get());
		if(newValue < problem.objectiveFunction(g, CStar.get())) {
			cout << "A better solution was found." << endl;
			CStar.reset();
			CStar = Cl;
			// TODO validar se essa saida eh valida: nao ha valor de FO menor que zero
			if(newValue == 0)  break;
		}
	}
	cout << "GRASP procedure done." << endl;
	return CStar;
}

ClusteringPtr Grasp::constructClustering(SignedGraph *g, float alpha, unsigned int ramdomSeed) {
	ClusteringPtr Cc = make_shared<Clustering>(g->getN()); // Cc = empty
	VertexSet lc(g->getN()); // L(Cc) = V(G)
	std::cout << "GRASP construct clustering...\n";

	while(lc.size() > 0) { // lc != empty
		// cout << "Vertex list size is " << lc.size() << endl;

		// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
		// It is important to calculate the modularity matrix first (used by vertex sorting)
		g->calculateModularityMatrix();
		lc.sort(g, Cc.get());

		// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		int i = lc.chooseRandomVertex(boost::math::iround(alpha * lc.size()));
		// std::cout << "Random vertex is " << i << std::endl;

		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		GainCalculation gainCalculation = Cc->gain(*g, i);
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

		// Cc->printClustering();
	}
	std::cout << "\nInitial clustering completed.\n";
	Cc->printClustering();
	return Cc;
}

ClusteringPtr Grasp::localSearch(SignedGraph *g, Clustering& Cc, int l,
		ClusteringProblem& problem) {
	// k is the current neighborhood distance in the local search
	int k = 1, iteration = 0;
	ClusteringPtr CStar = make_shared<Clustering>(Cc, g->getN()); // C* := Cc
	std::cout << "GRASP local search...\n";

	while(k <= l) {
		cout << "Local search iteration " << iteration << endl;
		// N := Nl(C*)
		// apply a local search in CStar using the k-neighborhood
		NeighborhoodList neig(g->getN());
		cout << "Generating neighborhood of size l = " << k << endl;
		ClusteringPtr Cl = neig.generateNeighborhood(k, g, CStar.get(), problem);
		cout << "Comparing local solution value." << endl;
		if(Cl.get() != NULL && CStar != NULL) {
			if(problem.objectiveFunction(g, Cl.get()) < problem.objectiveFunction(g, CStar.get())) {
				cout << "New local solution found." << endl;
				Cl->printClustering();
				CStar.reset();
				CStar = Cl;
				k = 1;
			} else {
				k++;
			}
		} else {  // no better result found in neighborhood
			k++;
		}
		iteration++;
	}
	std::cout << "GRASP local search done.\n";
	return CStar;
}

} /* namespace grasp */
} /* namespace resolution */
