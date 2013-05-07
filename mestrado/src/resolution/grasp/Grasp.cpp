/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#include "include/Grasp.h"
#include "../../graph/include/Clustering.h"
#include "include/VertexSet.h"

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
	VertexSet *lc = new VertexSet(g->getN()); // lc = VG

	while(lc->size() > 0) { // lc != empty
		// compute lc

		// choose i randomly among the first alpha elements of lc
		int i = 0;

		// c = c + {i}

		// lc = lc - {i}
		// the removal of vertex i recalculates the gain function
		lc->removeVertex(i);
	}
	return c;
}

Clustering* localSearch(SignedGraph* g, Clustering* c, int l) {
	Clustering* cl = c;
	Clustering* cStar = NULL;
	NeighborhoodList* neighborhood;
	do {
		cStar = cl;
		neighborhood = cStar->generateNeighborhood(l);
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
