/*
 * Graph.h
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <boost/config.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/shared_ptr.hpp>

// Maximum number of nodes in a graph
#define MAX_NODES 1000

using namespace boost;
using namespace std;

namespace clusteringgraph {

struct Edge {
    double weight;
    Edge() : weight(0) { }
    Edge(double w) : weight(w) { }
};
typedef adjacency_matrix<directedS, no_property, Edge > DirectedGraph;

class SignedGraph {
public:
	SignedGraph(const int &numberOfNodes);
	virtual ~SignedGraph();

	/**
	 * Returns the numbers of vertices of the graph.
	 */
	int getN();

	/**
	 * Return the number of edges of the graph
	 */
	int getM();

	/**
	 * Add an edge to the graph. Accepts only edges whose weight is
	 * equal to -1, 0 or 1.
	 */
	void addEdge(int a, int b, Edge edge);
	/**
	 * Returns the value of the corresponding edge.
	 */
	double getEdge(const int &a, const int &b);

	bool isPositiveEdge(const int &a, const int &b);

	bool isNegativeEdge(const int &a, const int &b);

	/**
	 * Returns the degree of vertex a.
	 */
	int getDegree(const int &a);

	/**
	 * Returns the negative degree of vertex a, that is, the sum of
	 * negative incoming edges.
	 */
	int getNegativeDegree(const int &a);

	void printGraph();

	string getGraphAsText();

	void setGraphAsText(string txt);

private:
	DirectedGraph graph;

	/* the number of nodes of the graph */
	int n;

	/**
	 * The text representation of the graph (for use on MPI messages).
	 */
	string graphAsText;
};

typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
