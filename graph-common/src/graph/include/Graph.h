/*
 * Graph.h
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#include "../../util/serialization/dynamic_bitset.hpp"

// Maximum number of nodes in a graph
// #define MAX_NODES 200000

using namespace boost;
using namespace std;

struct Edge {
    double weight;
    std::size_t vertex_index_t;
    Edge() : weight(0) { }
    Edge(double w) : weight(w) { }
};
struct Vertex {
    int id;
    std::size_t edge_index_t;
    Vertex() : id(0) { }
    Vertex(int w) : id(w) { }
};


enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, properties);
	BOOST_INSTALL_PROPERTY(edge, properties);
}

namespace clusteringgraph {

// typedef property< edge_index_t, size_t, Edge > EdgeProp;

typedef property< vertex_properties_t, Vertex,
            property< vertex_index_t, std::size_t > > VertexProperty;

typedef property< edge_properties_t, Edge, property< edge_index_t, std::size_t > > EdgeProperty;
// typedef property< edge_index_t, std::size_t, Edge > EdgeProperty;

// typedef property<vertex_index_t, Vertex> vertex_prop;
// typedef property<edge_index_t, Edge> edge_prop;

typedef adjacency_list< vecS, vecS, bidirectionalS, VertexProperty,
		EdgeProperty, no_property, vecS > DGraph;
typedef subgraph< DGraph > SubGraph;
typedef subgraph< DGraph > DirectedGraph;

/**
 *  uses dynamic_bitset for bool array, a high performance and space saving structure
 *  based on real bits
 *  the following array is initially empty and needs to be dynamically intialized.
 *  DISABLED.
 */
// typedef dynamic_bitset<> BoolArray;

// Defines the cluster list
// the list is made of boolean arrays, indicating that node i is in the cluster
// typedef std::vector<BoolArray> ClusterList;

typedef std::vector<long> ClusterArray;

class SignedGraph {
public:
	SignedGraph(const unsigned long &numberOfNodes);

	/**
	 * Builds a subgraph based on the graph g provided as parameter, induced by the
	 * vertex node list subGraphNodeList.
	 */
	SignedGraph(DirectedGraph &g, std::vector<long> subGraphNodeList);
	virtual ~SignedGraph();

	/**
	 * Returns the numbers of vertices of the graph.
	 */
	unsigned long getN();

	/**
	 * Returns the number of edges of the graph
	 */
	unsigned long getM();

	/**
	 * Return the id of the graph.
	 */
	unsigned int getId();

	void setId(const unsigned int& i);

	/**
	 * Add an edge to the graph. Accepts only edges whose weight is
	 * equal to -1, 0 or 1.
	 */
	void addEdge(unsigned long a, unsigned long b, Edge edge);

	/**
	 * Returns the degree of vertex a.
	 */
	unsigned long getDegree(const unsigned long &a);

	/**
	 * Returns the negative degree of vertex a, that is, the sum of
	 * negative incoming edges.
	 */
	unsigned long getNegativeDegree(const unsigned long &a);

	unsigned long getPositiveDegree(const unsigned long &a);

	/**
	 * Returns the negative edge cardinality between a vertex ni and a vertex set Si.
	 * Counts incoming and outcoming edges.
	 */
	double getNegativeEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray& cluster);

	/**
	 * Returns the positive edge cardinality between a vertex ni and a vertex set Si.
	 * Counts incoming and outcoming edges.
	 */
	double getPositiveEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray& cluster);

	string getGraphFileLocation();

	void setGraphFileLocation(string txt);

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


	DirectedGraph graph;

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

typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
