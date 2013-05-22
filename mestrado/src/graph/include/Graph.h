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
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

// Maximum number of nodes in a graph
#define MAX_NODES 1000

using namespace boost;

namespace clusteringgraph {

struct Edge {
    float weight;
    Edge() : weight(0) { }
    Edge(float w) : weight(w) { }
};
typedef adjacency_matrix<undirectedS, no_property, Edge > UndirectedGraph;
// the modularity matrix: a matrix of float
typedef multi_array<float, 2> ModularityMatrix;

class SignedGraph {
public:
	SignedGraph(int numberOfNodes);
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
	float getEdge(const int &a, const int &b);

	/**
	 * Returns the degree of vertex a.
	 */
	int getDegree(const int &a);

	/**
	 * Calculates the modularity matrix for this graph.
	 */
	void calculateModularityMatrix();

	ModularityMatrix& getModularityMatrix();

	void printGraph();

private:
	UndirectedGraph graph;
	/** the modularity matrix */
	ModularityMatrix modularityMatrix;
	bool modularityMatrixCalculated;
	/* the number of nodes of the graph */
	int n;
};

typedef boost::shared_ptr<SignedGraph> SignedGraphPtr;

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
