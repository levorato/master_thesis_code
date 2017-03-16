/*
 * BGLSignedGraph.h
 *
 *  Created on: 14 de mar de 2017
 *      Author: mlevorato
 */

#ifndef SRC_GRAPH_INCLUDE_BGLSIGNEDGRAPH_H_
#define SRC_GRAPH_INCLUDE_BGLSIGNEDGRAPH_H_

#include "Graph.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>

namespace clusteringgraph {

typedef adjacency_list< setS, vecS, undirectedS, VertexProperty,
		EdgeProperty, no_property, vecS > UDGraph;
typedef subgraph< UDGraph > SubGraph;
typedef subgraph< UDGraph > UndirectedGraph;


class BGLSignedGraph: public Graph {
public:
	BGLSignedGraph(const unsigned long &numberOfNodes);
	/**
	 * Builds a subgraph based on the graph g provided as parameter, induced by the
	 * vertex node list subGraphNodeList.
	 */
	BGLSignedGraph(UndirectedGraph &g, std::vector<long> subGraphNodeList);
	virtual ~BGLSignedGraph();

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

	/**
	 *  Creates a string to be used as input format for the Graclus clustering program,
	 *  which clusters undirected graphs.
	 *  If edge weights (must be integer values) are different,
			  10
		  1 ----- 2
		  |	  |
		 9|	  |6
		  |   7   |
		  3 ----- 5
		   \     /
		  11\   /28
			 \ /
			  4

		 then the matrix representation becomes

		 5 6 1		<--- # of nodes and edges and format
		 2 10 3 9	<--- nodes adjacent to 1 and corresponding edge weight
		 1 10 5 6	.
		 1 9 4 11 5 7	.
		 3 11 5 28	.
		 2 6 3 7 4 28	<--- nodes adjacent to 5 and corresponding edge weight
	 *
	 */
	string convertToGraclusInputFormat();

	UndirectedGraph graph;

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

#endif /* SRC_GRAPH_INCLUDE_BGLSIGNEDGRAPH_H_ */
