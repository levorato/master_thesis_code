/*
 * ImbalanceSubgraphParallelILS.h
 *
 *  Created on: Mar 10, 2017
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_ILS_SPLITGRAPH_PARALLELGRAPHDISTRIBUTEDILS_H_
#define SRC_RESOLUTION_ILS_SPLITGRAPH_PARALLELGRAPHDISTRIBUTEDILS_H_

#include "../../include/ParallelILS.h"
#include "ImbalanceSubgraphParallelILS.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <list>

#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "SplitgraphUtil.h"
#include "util/include/MPIMessage.h"

using namespace boost::numeric::ublas;

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Defines the parameters used in global VND neighborhoods
#define PERCENTAGE_OF_PROCESS_PAIRS 0.2
#define PERCENTAGE_OF_MOST_IMBALANCED_VERTICES_TO_BE_MOVED 0.25
#define PERCENTAGE_OF_MOST_IMBALANCED_CLUSTERS_TO_BE_MOVED 0.05
// Defines the number of times each neighborhood structure will rerun local ILS while global solution does not improve
#define MAX_ILS_RETRIES 2
// Defines the time limit (in seconds) of the local ILS procedure on each node
#define LOCAL_ILS_TIME_LIMIT 60
#define EPS 10e-6

namespace resolution {
namespace ils {

/**
 * This class is responsible for all logic related to the SplitGraph Distributed ILS algorithm,
 * based on Boost Parallel BGL.
 * Splits a large graph into smaller subgraphs, distributes them among several processes,
 * invokes a local ILS procedure for each one, and then merges all partial solutions into a global one.
 * After that, the algorithm iteratively tries to improve the global solution by moving clusters
 * between processes, using 4 different neighborhood structures.
 */
class ParallelGraphDistributedILS: public ImbalanceSubgraphParallelILS {
public:
	ParallelGraphDistributedILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
			const bool& split = true, const bool& cuda = true, const bool& parallelgraph = true);
	virtual ~ParallelGraphDistributedILS();

	/**
	 * Triggers the parallel execution of the ILS algorithm using Imbalance ratio / SplitGraph.
	 */
	virtual Clustering executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	ProcessClustering preProcessSplitgraphPartitioning(SignedGraph *g, ClusteringProblem& problem,
			bool partitionByVertex);

	/**
	 * Triggers the distributed ILS (split graph) resolution, invoking the local ILS
	 * resolution in each process (distributeSubgraphsBetweenProcessesAndRunILS)
	 * as many times as needed to improve the global solution.
	 */
	Clustering runDistributedILS(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
			ProcessClustering& splitgraphClustering, Clustering& currentSolution);

	/**
	 * Executes a local ILS in each subgraph / process, gathers the individual results and
	 * merges them into a global CC solution.
	 */
	Clustering distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
			ProcessClustering& splitgraphClustering);

	/**
	 * Executes 2 local ILS procedures for each cluster movement (2 subgraphs), gathers the individual results,
	 * merges them into a global CC solution and returns the best movement between all.
	 */
	Clustering distributeClusterMovementsAndRunILS(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
			ClusterArray& splitgraphClusterArray, ImbalanceMatrix& processClusterImbMatrix);

	/**
	 * Rebalances clusters between processes without running ILS (zero-cost moves).
	 */
	void rebalanceClustersBetweenProcessesWithZeroCost(SignedGraph* g, ClusteringProblem& problem,
			ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering,	const int& numberOfProcesses);

	/**
	 * 1-move-cluster: moves a cluster from one process to another.
	 */
	bool moveCluster1opt(SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * 	swap-cluster: in this neighborhood the priority is swapping a bigger cluster in an overloaded process
	 * 	with a smaller cluster from a less-loaded process.
	 */
	bool swapCluster1opt(SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * 2-move-cluster: try to move 2 clusters from an overloaded process to 2 less-loaded processes.
	 */
	bool twoMoveCluster(SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * quasi-clique-move or pseudo-clique-move: moves a quasi-clique from an existing cluster to
	 * a new cluster in another process.
	 */
	bool splitClusterMove(SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering,
			const int& numberOfProcesses,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info);

	/**
	 * Distributed VND: alternates between 4 types of neighborhood structures that move clusters
	 * between processes, in a distributed fashion.
	 */
	long variableNeighborhoodDescent(SignedGraph* g, ProcessClustering& bestSplitgraphClustering,
			Clustering& bestClustering, const int& numberOfProcesses,
			ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
			ClusteringProblem& problem, ExecutionInfo& info, const double& timeSpentSoFar, int invocationNumber);

protected:
	int machineProcessAllocationStrategy;
	unsigned int numberOfSearchSlaves;
	unsigned int numberOfSlaves;
	Clustering CCclustering;
	bool splitGraph;
	bool cudaEnabled;
	// counts the number of times the local ILS found solutions worse than the current solution (worse than zero-cost move).
	long numberOfFrustratedSolutions;

	// data structures containing the imbalance contribution of each vertex and between processes
	std::vector< pair<long, double> > vertexImbalance;

	// time spent by each process in each iteration
	std::vector< std::vector<double> > timeSpentAtIteration;

	SplitgraphUtil util;

	/**
	 * Executes ILS locally (execution performed by leader process, rank 0) on the subgraph of g, induced by vertexList.
	 */
	OutputMessage runILSLocallyOnSubgraph(ConstructClustering *construct,
			VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
			const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info, std::vector<long>& vertexList);

	void moveClusterToDestinationProcessZeroCost(SignedGraph *g, Clustering& bestClustering,
			ProcessClustering& bestSplitgraphClustering, long clusterToMove, unsigned int sourceProcess, unsigned int destinationProcess);
};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_PARALLELGRAPHDISTRIBUTEDILS_H_ */
