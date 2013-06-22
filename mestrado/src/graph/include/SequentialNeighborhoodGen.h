/*
 * SequentialNeighborhoodGen.h
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#ifndef SEQUENTIALNEIGHBORHOODGEN_H_
#define SEQUENTIALNEIGHBORHOODGEN_H_

#include "Neighborhood.h"

namespace clusteringgraph {

class SequentialNeighborhoodGenerator: public clusteringgraph::NeighborhoodListGenerator {
public:
	SequentialNeighborhoodGenerator(int n);
	virtual ~SequentialNeighborhoodGenerator();

	virtual ClusteringPtr generateNeighborhood(int l, SignedGraph* g,
				Clustering* clustering, const ClusteringProblem& problem,
				double timeSpentSoFar, double timeLimit, unsigned long randomSeed);
};

} /* namespace clusteringgraph */
#endif /* SEQUENTIALNEIGHBORHOODGEN_H_ */
