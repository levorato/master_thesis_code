/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#include "include/Grasp.h"

namespace resolution {
namespace grasp {

Grasp::Grasp() {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

Clustering* Grasp::executeGRASP(SignedGraph* g, int iter, float alpha, int l) {

}

Clustering* Grasp::constructClustering(SignedGraph* g, float alpha) {
	Clustering *c = NULL; // c = empty
	Clustering *lc = new Clustering(1, g->getN()); // lc = VG

	while(lc->getN() > 0) { // lc != empty
		// compute lc

		// choose i ramdomly among the first alpha elements of lc
		// and recalculates the objective function

		// c = c + {i}

		// lc = lc - {i}

	}
	return c;
}

Clustering* localSearch(SignedGraph* g, Clustering* c, int l) {
	Clustering* cl = c;
	Clustering* cStar = NULL;
	std::vector<Clustering*> neighborhood;
	do {
		cStar = cl;
		neighborhood = Clustering::generateNeighborhood(cStar, l);
		for each Clustering* c in neighborhood do {
			if(Q(c) > Q(cl))
				cl = c;
			// N = N - {c}

		}
	} while(cl == cStar); // TODO sobecarregar funcao de igual

	return cl;
}

} /* namespace grasp */
} /* namespace resolution */
