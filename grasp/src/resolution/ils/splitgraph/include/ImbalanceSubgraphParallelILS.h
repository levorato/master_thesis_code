/*
 * ImbalanceSubgraphParallelILS.h
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_
#define SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_

#include "../../include/ParallelILS.h"

namespace resolution {
namespace ils {

class ImbalanceSubgraphParallelILS: public ILS {
public:
	ImbalanceSubgraphParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
			const bool& split = false, const bool& cuda = true);
	virtual ~ImbalanceSubgraphParallelILS();

	/**
	 * Triggers the parallel execution of the ILS algorithm using Imbalance ratio / SplitGraph.
	 */
	virtual Clustering executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

protected:
	int machineProcessAllocationStrategy;
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
	Clustering CCclustering;
	bool splitGraph;
	bool cudaEnabled;

};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_ */
