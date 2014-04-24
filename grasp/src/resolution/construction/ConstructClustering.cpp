/*
 * ConstructClustering.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: mlevorato
 */

#include "include/ConstructClustering.h"
#include "../construction/include/VertexSet.h"
#include <boost/log/trivial.hpp>
#include <boost/math/special_functions/round.hpp>


namespace resolution {
namespace construction {

ConstructClustering::ConstructClustering(GainFunction* f, const unsigned long& seed,
		const double& alpha) : gainFunction(f), randomSeed(seed), _alpha(alpha) {

}

ConstructClustering::~ConstructClustering() {

}

Clustering ConstructClustering::constructClustering(SignedGraph *g,
		ClusteringProblem& problem, const int& myRank) {
	Clustering Cc; // Cc = empty
	VertexSet lc(randomSeed, g->getN()); // L(Cc) = V(G)
	BOOST_LOG_TRIVIAL(debug)<< "GRASP construct clustering...\n";

	while (lc.size() > 0) { // lc != empty
		GainCalculation gainCalculation;
		if (_alpha == 1.0) {
			// alpha = 1.0 (completely random): no need to calculate all gains (saves time)
			int i = lc.chooseRandomVertex(lc.size()).vertex;
			gainCalculation = gainFunction->calculateIndividualGain(problem, Cc,	i);
		} else {
			// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
			// according to the value of the gain function
			gainFunction->calculateGainList(problem, Cc, lc.getVertexList());

			// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
			// (alpha x |lc|) is a rounded number
			// IMPORTANT: if alpha < 0, the constructClustering method will always
			//             choose the first vertex in the gainFunction list, that is,
			//             the one that minimizes the objective (VOTE algorithm).
			if (_alpha <= 0.0) {
				// chooses first element (min) without ordering (saves time)
				gainCalculation = lc.chooseFirstElement(gainFunction);
			} else {
				lc.sort(gainFunction);
				gainCalculation = lc.chooseRandomVertex(
						boost::math::iround(_alpha * lc.size()));
			}
			// std::cout << "Random vertex between 0 and " << boost::math::iround(alpha * lc.size()) << " is " << i << std::endl;
		}
		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		// cout << "Selected vertex is " << i << endl;
		if (gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			Cc.addCluster(*g, problem, gainCalculation.vertex);
		} else {
			// inserts i into existing cluster
			Cc.addNodeToCluster(*g, problem, gainCalculation.vertex, gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the choosing vertex i automatically removes it from the list
		// Removal already done by the chooseVertex methods above
		// Cc->printClustering();
	}
	BOOST_LOG_TRIVIAL(debug)<< myRank << ": Initial clustering completed.\n";
	Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	// Cc.printClustering();
	return Cc;
}

} /* namespace construction */
} /* namespace resolution */
