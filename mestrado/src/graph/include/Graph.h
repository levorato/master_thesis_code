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
#include <boost/graph/graph_utility.hpp>
#include <boost/shared_ptr.hpp>

// Maximum number of nodes in a graph
// #define MAX_NODES 200000

using namespace boost;
using namespace std;

namespace clusteringgraph {

struct Edge {
    double weight;
    Edge() : weight(0) { }
    Edge(double w) : weight(w) { }
};
struct Vertex {
    int id;
    Vertex() : id(0) { }
    Vertex(int w) : id(w) { }
};
typedef adjacency_list<vecS, vecS, bidirectionalS, Vertex, Edge, no_property, vecS > DirectedGraph;

class SignedGraph {
public:
	SignedGraph(const unsigned long &numberOfNodes);
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

	string getGraphAsText();

	void setGraphAsText(string txt);

	DirectedGraph graph;
private:


	/** the number of nodes of the graph */
	unsigned long n;

	/** the identifier of the graph */
	unsigned int id;

	/**
	 * The text representation of the graph (for use on MPI messages).
	 */
	string graphAsText;
};

typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
