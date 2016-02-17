/*
 * SplitgraphUtil.h
 *
 *  Created on: Feb 17, 2016
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_SPLITGRAPHUTIL_H_
#define SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_SPLITGRAPHUTIL_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <list>
#include <vector>

#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "graph/include/Imbalance.h"

using namespace boost::numeric::ublas;
using namespace clusteringgraph;

/** Clique detection parameters */
#define POS_EDGE_PERC_RELAX 0.7  // minimum percentage of positive edges in pseudo-clique
#define NEG_EDGE_PERC_RELAX 0.3  // maximum percentage of negative edges in pseudo-clique
// Defines the parameters used in pseudo-clique (splitcluster) neighborhood
#define CLUSTERING_ALPHA 0.2  // pseudo-clique alpha constructive parameter

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

struct weighted_pos_neg_degree_ordering_desc
{
    inline bool operator() (const VertexDegree& struct1, const VertexDegree& struct2)
    {
    	// Sorts the vertex list according to a weighted formula of positive and negative degrees (descending order)
    	double p1 = (POS_EDGE_PERC_RELAX * struct1.positiveDegree) - (NEG_EDGE_PERC_RELAX * struct1.negativeDegree);
    	double p2 = (POS_EDGE_PERC_RELAX * struct2.positiveDegree) - (NEG_EDGE_PERC_RELAX * struct2.negativeDegree);
    	return (p1 > p2);
    }
};

struct ImbalanceMatrix {
	matrix<double> pos, neg;
	ImbalanceMatrix() : pos(), neg() { }
	ImbalanceMatrix(int nc) : pos(zero_matrix<double>(nc, nc)), neg(zero_matrix<double>(nc, nc)) { }
};

class SplitgraphUtil {
public:
	SplitgraphUtil();
	virtual ~SplitgraphUtil();

	ImbalanceMatrix calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster,
			std::vector< pair<long, double> >& vertexImbalance, const int& numberOfProcesses);

	void updateProcessToProcessImbalanceMatrix(SignedGraph& g, const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix, const int& numberOfProcesses);

	std::vector<Coordinate> obtainListOfClustersFromProcess(SignedGraph& g,
				Clustering& globalClustering, int processNumber);

	std::vector<long> getListOfVeticesInCluster(SignedGraph& g, Clustering& globalClustering,
			long clusterNumber);

	Imbalance calculateExternalImbalanceSumBetweenProcesses(ImbalanceMatrix& processClusterImbMatrix);

	Imbalance calculateInternalImbalanceSumOfAllProcesses(std::vector<Imbalance>& internalProcessImbalance);

	Imbalance calculateProcessInternalImbalance(SignedGraph *g, Clustering& c, unsigned int processNumber);

	/**
	  *  Calculates the positive and negative degrees of each vertex v in clusterX
	  *  *relative to clusterX only*.
	*/
	std::list<VertexDegree> calculateDegreesInsideCluster(SignedGraph *g,
			Clustering& globalClustering, long clusterX);

	long chooseRandomVertex(std::list<VertexDegree>& vertexList, long x);

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

	std::vector<Coordinate> obtainListOfOverloadedProcesses(SignedGraph& g,
					Clustering& splitGraphClustering, long maximumNumberOfVertices);

	/**
	  *  Returns a list containing the vertices that belong to the pseudo clique C+,
	  *  found inside global cluster X. This is a greedy heuristic.
	*/
	std::vector<long> findPseudoCliqueC(SignedGraph *g, Clustering& globalClustering,
			long clusterX);

	/**
	  *  Determines if the set cliqueC \union {u} is a pseudo clique in the graph g.
	  *  By definition, cliqueC is already a pseudo clique.
	*/
	bool isPseudoClique(SignedGraph *g, std::list<long>& cliqueC,
			ClusterArray& cliqueCClusterArray, long u);

	/**
	  *  Determines if the set cliqueC \union {u, v} is a pseudo clique in the graph g.
	  *  By definition, cliqueC is already a pseudo clique.
	*/
	bool isPseudoClique(SignedGraph *g, ClusterArray& cliqueCClusterArray, long cliqueSize,
			long u, long v);

};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_SPLITGRAPHUTIL_H_ */
