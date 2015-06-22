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

struct VertexDegree {
    long id;
    double positiveDegree;
    double negativeDegree;
    VertexDegree() : id(0), positiveDegree(0.0), negativeDegree(0.0) { }
    VertexDegree(long vid, double posdeg, double negdeg) : id(vid), positiveDegree(posdeg), negativeDegree(negdeg) { }
};

struct neg_degree_ordering_asc
{
    inline bool operator() (const VertexDegree& struct1, const VertexDegree& struct2)
    {
        return (struct1.negativeDegree < struct2.negativeDegree);
    }
};

struct pos_degree_ordering_asc
{
    inline bool operator() (const VertexDegree& struct1, const VertexDegree& struct2)
    {
        return (struct1.positiveDegree < struct2.positiveDegree);
    }
};

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
