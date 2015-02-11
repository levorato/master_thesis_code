/*
 * CUDAConstructClustering.h
 *
 *  Created on: 09/02/2015
 *      Author: czt0
 */

#ifndef CUDACONSTRUCTCLUSTERING_H_
#define CUDACONSTRUCTCLUSTERING_H_

#include "ConstructClustering.h"
#include "CUDAImbalanceGainFunction.h"

namespace resolution {
namespace construction {

class CUDAConstructClustering : public ConstructClustering {
public:
	CUDAConstructClustering(CUDAImbalanceGainFunction *f, const unsigned long& seed);
	virtual ~CUDAConstructClustering();

	/**
	 * Constructs a clustering in a greedy randomized fashion,
	 * starting from the empty set. CUDAConstructClustering is based on alpha = 1.0, that is,
	 * it chooses a random vertex to insert in the set and selects the destination cluster
	 * based on the best imbalance (calculated with CUDA kernel).
	 * WARNING: CUDAConstructClustering is ONLY for the CC problem (not RCC)!
	 * This is the first phase of the metaheuristic.
	 * @param g graph to be used as the base
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @return Clustering Cc
	 */
	virtual Clustering constructClustering(SignedGraph *g, ClusteringProblem& problem,
			const int& myRank);

	/**
	 * CUDAConstructClustering uses the imbalance gain function.
	 */
	int getGainFunctionType() {
		return GainFunction::IMBALANCE;
	}

};

} /* namespace construction */
} /* namespace resolution */
#endif /* CUDACONSTRUCTCLUSTERING_H_ */
