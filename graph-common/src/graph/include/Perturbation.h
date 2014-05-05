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

namespace clusteringgraph {

class Perturbation {
public:
	Perturbation(unsigned long randomSeed) : _randomSeed(randomSeed) {

	}

	virtual ~Perturbation();

	Clustering randomMove(SignedGraph* g, Clustering clustering, ClusteringProblem& p,
			unsigned long numberOfMoves);

private:
    Clustering move1opt(SignedGraph* g, Clustering clustering, ClusteringProblem& p);

    unsigned long _randomSeed;
};

} /* namespace clusteringgraph */

#endif /* PERTURBATION_H_ */
