/*
 * Graph.h
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <boost/config.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>

// Maximum number of nodes in a graph
#define MAX_NODES 1000

using namespace boost;

namespace clusteringgraph {

struct Edge {
    float weight;
    Edge() : weight(0) { }
    Edge(float w) : weight(w) { }
};
typedef adjacency_matrix<directedS, no_property, Edge > DirectedGraph;
typedef boost::scoped_ptr<DirectedGraph> DigraphPtr;

class SignedGraph {
public:
	SignedGraph(int numberOfNodes);
	virtual ~SignedGraph();

	/**
	 * Returns the numbers of verticves of the graph.
	 */
	int getN();

	/**
	 * Add an edge to the graph. Accepts only edges whose weight is
	 * equal to -1, 0 or 1.
	 */
	void addEdge(int a, int b, Edge edge);
	/**
	 * Return the value of the corresponding edge.
	 */
	float getEdge(int a, int b);
	void printGraph();

private:
	DigraphPtr digraphPtr;

};

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
