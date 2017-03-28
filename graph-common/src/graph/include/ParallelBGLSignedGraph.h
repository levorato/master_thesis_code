/*
 * ParallelBGLSignedGraph.h
 *
 *  Created on: 14 de mar de 2017
 *      Author: mlevorato
 */

#ifndef SRC_GRAPH_INCLUDE_PARALLELBGLSIGNEDGRAPH_H_
#define SRC_GRAPH_INCLUDE_PARALLELBGLSIGNEDGRAPH_H_

#include "Graph.h"

// includes from parallel BGL
#include <boost/mpi.hpp>
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/local_subgraph.hpp>

namespace clusteringgraph {

using namespace boost;
using boost::graph::distributed::mpi_process_group;

typedef boost::adjacency_list<vecS, distributedS<mpi_process_group, vecS>, boost::undirectedS,
		VertexProperty, EdgeProperty, no_property, vecS > ParallelGraph;

class ParallelBGLSignedGraph: public Graph {
public:
	ParallelBGLSignedGraph(const unsigned long &numberOfNodes, ParallelGraph *pgraph);
	~ParallelBGLSignedGraph();

	/**
	 * Returns the number of vertices of the graph.
	 */
	virtual unsigned long getN();

	/**
	 * Sets the number of vertices of the graph.
	 */
	void setN(long num_vertices);

	/**
	 * Returns the number of edges of the graph
	 */
	virtual unsigned long getM();

	/**
	 * Return the id of the graph.
	 */
	virtual unsigned int getId();

	virtual void setId(const unsigned int& i);

	/**
	 * Add an edge to the graph. Accepts only edges whose weight is
	 * equal to -1, 0 or 1.
	 */
	virtual void addEdge(unsigned long source, unsigned long target, Edge edge);

	/**
	 * Returns the degree of vertex a.
	 */
	virtual unsigned long getDegree(const unsigned long &a);

	/**
	 * Returns the out-degree of vertex a.
	 */
	virtual unsigned long getOutDegree(const unsigned long &a);

	/**
	 * Returns the negative degree of vertex a, that is, the sum of
	 * negative incoming edges.
	 */
	virtual unsigned long getNegativeDegree(const unsigned long &a);

	virtual unsigned long getPositiveDegree(const unsigned long &a);

	/**
	 * Returns the negative edge cardinality between a vertex ni and a vertex set Si.
	 * Counts incoming and outcoming edges.
	 */
	virtual double getNegativeEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray& cluster);

	/**
	 * Returns the positive edge cardinality between a vertex ni and a vertex set Si.
	 * Counts incoming and outcoming edges.
	 */
	virtual double getPositiveEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray& cluster);

	/**
	 * Returns the number of edges crossing a specific cluster and also internal to the same cluster.
	 */
	virtual long getNumberOfEdgesInClustering(const ClusterArray& cluster, const long& clusterNumber);

	virtual string getGraphFileLocation();

	virtual void setGraphFileLocation(string txt);

	ParallelGraph *graph;

	ParallelGraph& getGraph() {
		return *graph;
	}

private:


	/** the number of nodes of the graph */
	unsigned long n;

	/** the identifier of the graph */
	unsigned int id;

	/**
	 * The file location of the graph (for use on MPI messages).
	 */
	string graphFileLocation;
};

typedef ParallelBGLSignedGraph SignedGraph;
typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;
typedef local_subgraph<ParallelGraph> LocalSubgraph;

} /* namespace clusteringgraph */

#endif /* SRC_GRAPH_INCLUDE_PARALLELBGLSIGNEDGRAPH_H_ */
