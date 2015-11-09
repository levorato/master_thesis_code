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
#include <list>

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
        if(struct1.value == struct2.value) {
        	return struct1.y < struct2.y;
        }
    	return (struct1.value < struct2.value);
    }
};

struct coordinate_ordering_desc
{
    inline bool operator() (const Coordinate& struct1, const Coordinate& struct2)
    {
        if(struct1.value == struct2.value) {
        	return struct1.y > struct2.y;
        }
    	return (struct1.value > struct2.value);
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

struct pos_degree_ordering_desc
{
    inline bool operator() (const VertexDegree& struct1, const VertexDegree& struct2)
    {
        return (struct1.positiveDegree > struct2.positiveDegree);
    }
};

struct ImbalanceMatrix {
	matrix<double> pos, neg;
	ImbalanceMatrix() : pos(), neg() { }
	ImbalanceMatrix(int nc) : pos(zero_matrix<double>(nc, nc)), neg(zero_matrix<double>(nc, nc)) { }
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

	Clustering preProcessSplitgraphPartitioning(SignedGraph *g, ClusteringProblem& problem, bool partitionByVertex);

	Clustering distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
			ClusterArray& splitgraphClusterArray, ImbalanceMatrix& processClusterImbMatrix);

	ImbalanceMatrix calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster);

	void updateProcessToProcessImbalanceMatrix(SignedGraph& g, const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix);

	bool moveVertex1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses, ImbalanceMatrix& processClusterImbMatrix,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	bool moveCluster1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses, ImbalanceMatrix& processClusterImbMatrix,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * 	swap-cluster: in this neighborhood the priority is swapping a bigger cluster in an overloaded process
	 * 	with a smaller cluster from a less-loaded process.
	 */
	bool swapCluster1opt(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses, ImbalanceMatrix& processClusterImbMatrix,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * 2-move-cluster: try to move 2 clusters from an overloaded process to 2 less-loaded processes.
	 */
	bool twoMoveCluster(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses, ImbalanceMatrix& processClusterImbMatrix,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 *
	 */
	bool positiveCliqueMove(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses, ImbalanceMatrix& processClusterImbMatrix,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	long variableNeighborhoodDescent(SignedGraph* g, Clustering& bestSplitgraphClustering,
			Clustering& bestClustering, const int& numberOfProcesses,
			ImbalanceMatrix& processClusterImbMatrix, ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info, const double& timeSpentSoFar);

protected:
	int machineProcessAllocationStrategy;
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
	Clustering CCclustering;
	bool splitGraph;
	bool cudaEnabled;

	// data structures containing the imbalance contribution of each vertex and between processes
	std::vector< pair<long, double> > vertexImbalance;

	// time spent by each process in each iteration
	std::vector< std::vector<double> > timeSpentAtIteration;

	Coordinate findMaximumElementInMatrix(ImbalanceMatrix &mat);

	std::vector< Coordinate > getMatrixElementsAsList(ImbalanceMatrix &mat);

	long findMostImbalancedVertexInProcessPair(SignedGraph& g, ClusterArray& splitGraphCluster,
			ClusterArray& globalCluster, Coordinate processPair);

	std::vector<Coordinate> obtainListOfImbalancedClusters(SignedGraph& g,
			ClusterArray& splitGraphCluster, Clustering& globalClustering);

	/**
	 * A vertex-overloaded process is a process with more than (n / numberOfProcesses) vertices.
	 */
	std::vector<Coordinate> obtainListOfOverloadedProcesses(SignedGraph& g,
				Clustering& splitGraphClustering);

	std::vector<Coordinate> obtainListOfClustersFromProcess(SignedGraph& g,
				Clustering& globalClustering, int processNumber);

	std::vector<long> getListOfVeticesInCluster(SignedGraph& g, Clustering& globalClustering,
			long clusterNumber);

	/**
	  *  Returns a list containing the vertices that belong to the positive clique C+,
	  *  found inside global cluster X. This is a greedy heuristic.
	*/
	std::vector<long> findPositiveCliqueC(SignedGraph *g, Clustering& globalClustering,
			long clusterX, double alpha);

	/**
	  *  Determines if the set cliqueC \union {u} is a positive clique in the graph g.
	  *  By definition, cliqueC is already a positive clique.
	*/
	bool isPositiveClique(SignedGraph *g, std::list<long>& cliqueC,
			ClusterArray& cliqueCClusterArray, long u);

	/**
	  *  Determines if the set cliqueC \union {u, v} is a positive clique in the graph g.
	  *  By definition, cliqueC is already a positive clique.
	*/
	bool isPositiveClique(SignedGraph *g, ClusterArray& cliqueCClusterArray, long cliqueSize,
			long u, long v);

	/**
	  *  Calculates the positive and negative degrees of each vertex v in clusterX
	  *  *relative to clusterX only*.
	*/
	std::list<VertexDegree> calculateDegreesInsideCluster(SignedGraph *g,
			Clustering& globalClustering, long clusterX);

	long chooseRandomVertex(std::list<VertexDegree>& vertexList, int x);

};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_IMBALANCESUBGRAPHPARALLELILS_H_ */
