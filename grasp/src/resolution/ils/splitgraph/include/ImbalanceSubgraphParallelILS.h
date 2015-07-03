/*
 * ImbalanceSubgraphParallelILS.h
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_
#define SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_

#include "../../include/ParallelILS.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace resolution {
namespace ils {

struct Coordinate {
	int x, y;
	double value;
	Coordinate() : x(0), y(0), value(0.0) { }
	Coordinate(int a, int b, double vl = 0.0) : x(a), y(b), value(vl) { }
};

struct coordinate_ordering_asc
{
    inline bool operator() (const Coordinate& struct1, const Coordinate& struct2)
    {
        return (struct1.value < struct2.value);
    }
};

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

	Clustering preProcessSplitgraphPartitioning(SignedGraph *g, ClusteringProblem& problem);

	Clustering distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
			ClusterArray& splitgraphClusterArray);

	void calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster);

	Coordinate findMaximumElementInMatrix(matrix<double> &mat);

	std::vector< Coordinate > getMatrixElementsAsList(matrix<double> &mat);

	long findMostImbalancedVertexInProcessPair(SignedGraph& g, ClusterArray& myCluster, Coordinate processPair);


protected:
	int machineProcessAllocationStrategy;
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
	Clustering CCclustering;
	bool splitGraph;
	bool cudaEnabled;

	// data structures containing the imbalance contribution of each vertex and between processes
	std::vector< pair<long, double> > vertexImbalance;
	matrix<double> clusterImbMatrix;

};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_ */
