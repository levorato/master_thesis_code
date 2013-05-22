/*
 * Graph.cpp
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#include "include/Graph.h"
#include <boost/config.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <cassert>

/*
 * Defines the boost::graph_traits class template.
 */
#include <boost/graph/graph_traits.hpp>

/*
 * Defines the vertex and edge property tags.
 */
#include <boost/graph/properties.hpp>

using namespace std;
using namespace boost;

namespace clusteringgraph {

SignedGraph::SignedGraph(int numberOfNodes) : graph(numberOfNodes),
		modularityMatrix(boost::extents[numberOfNodes][numberOfNodes]),
		modularityMatrixCalculated(false), n(numberOfNodes) {

}

SignedGraph::~SignedGraph() {
	// TODO Auto-generated destructor stub
}


int SignedGraph::getN() {
	cout << "Graph n is " << n << endl;
	return n;
}

int SignedGraph::getM() {
	return graph.m_num_edges;
}

void SignedGraph::addEdge(int a, int b, Edge edge) {
	add_edge(a, b, edge, graph);
}

float SignedGraph::getEdge(const int &a, const int &b) {
	return ((Edge)graph.get_edge(a, b).second).weight;
}

int SignedGraph::getDegree(const int &a) {
	return in_degree(a, graph);
}

// TODO calculate the modularity matrix of weighed graphs
void SignedGraph::calculateModularityMatrix() {
	int m = this->getM();
	int numberOfNodes = this->getN();
	int degree[numberOfNodes];
	// Prestore the degrees for optimezed lookup
	for(int i = 0; i < numberOfNodes; i++) {
		degree[i] = this->getDegree(i);
	}

	for(int i = 0; i < numberOfNodes; i++) {
		for(int j = 0; j < numberOfNodes; j++) {
			int a = (this->getEdge(i, j) != 0) ? 1 : 0;
			modularityMatrix[i][j] = a - ( (degree[i] * degree[j]) / (2 * m) );
		}
	}
	modularityMatrixCalculated = true;
}

ModularityMatrix& SignedGraph::getModularityMatrix() {
	assert(modularityMatrixCalculated);
	return modularityMatrix;
}

void SignedGraph::printGraph() {
	const char* name = "ABCDEF";

	int numberOfInEdges = boost::in_degree(1, graph);
	std::cout << "numberOfInEdges: " << numberOfInEdges << std::endl;

	cout << "vertex set: ";
	print_vertices(graph, name);
	cout << std::endl;

	cout << "edge set: ";
	print_edges(graph, name);
	cout << std::endl;

	cout << "out-edges: " << std::endl;
	print_graph(graph, name);
	cout << std::endl;
}

} /* namespace graph */
