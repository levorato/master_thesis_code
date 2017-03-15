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


namespace clusteringgraph {

typedef boost::adjacency_list<boost::setS, boost::distributedS<boost::graph::distributed::mpi_process_group, boost::vecS>, boost::undirectedS,
		VertexProperty, EdgeProperty, no_property, vecS > ParallelGraph;


class ParallelBGLSignedGraph: public SignedGraph {
public:
	ParallelBGLSignedGraph(const unsigned long &numberOfNodes);
	virtual ~ParallelBGLSignedGraph();

	/**
	 * Returns the numbers of vertices of the graph.
	 */
	virtual unsigned long getN();

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
	virtual void addEdge(unsigned long a, unsigned long b, Edge edge);

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

	ParallelGraph graph;

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

} /* namespace clusteringgraph */

#endif /* SRC_GRAPH_INCLUDE_PARALLELBGLSIGNEDGRAPH_H_ */
