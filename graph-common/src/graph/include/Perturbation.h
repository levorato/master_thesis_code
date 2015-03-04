/*
 * Perturbation.h
 *
 *  Created on: Apr 25, 2014
 *      Author: mlevorato
 */

#ifndef PERTURBATION_H_
#define PERTURBATION_H_

#include "Clustering.h"
#include "Graph.h"
#include "../../problem/include/ClusteringProblem.h"
#include "../../problem/include/RCCProblem.h"

#include <utility>
#include <vector>

namespace clusteringgraph {

class Perturbation {
public:
	Perturbation(unsigned long randomSeed) : _randomSeed(randomSeed) {

	}

	virtual ~Perturbation();

	Clustering randomMove(SignedGraph* g, Clustering clustering, ClusteringProblem& p,
			unsigned long numberOfMoves);

	std::pair< std::vector<int>, std::vector<int> > generateRandomMoveListCCProblem(int n,
			int nc, unsigned long numberOfMoves);

	std::pair<int, int> randomMove1optCCProblem(int n, int nc);

private:
    Clustering randomMove1opt(SignedGraph* g, Clustering clustering, ClusteringProblem& p);

    unsigned long _randomSeed;
    static const int NEW_CLUSTER = -1;
};

} /* namespace clusteringgraph */

#endif /* PERTURBATION_H_ */
