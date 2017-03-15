/*
 * Graph.h
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <boost/config.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
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
	virtual SignedGraph(const unsigned long &numberOfNodes);

	/**
	 * Builds a subgraph based on the graph g provided as parameter, induced by the
	 * vertex node list subGraphNodeList.
	 */
	SignedGraph(UndirectedGraph &g, std::vector<long> subGraphNodeList);
	virtual ~SignedGraph();

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

};

typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
