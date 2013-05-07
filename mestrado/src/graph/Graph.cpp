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

SignedGraph::SignedGraph(int numberOfNodes) : digraphPtr(new DirectedGraph(numberOfNodes)) {

}

SignedGraph::~SignedGraph() {
	// TODO Auto-generated destructor stub
}


int SignedGraph::getN() {
	return num_vertices(*digraphPtr);
}

void SignedGraph::addEdge(int a, int b, int value) {
	if(value == -1 || value == 0 || value == -1)
		add_edge(a, b, Edge(value), *digraphPtr);
}

void SignedGraph::printGraph() {
	const char* name = "ABCDEF";

	int numberOfInEdges = boost::in_degree(1,*digraphPtr);
	std::cout << "numberOfInEdges: " << numberOfInEdges << std::endl;

	cout << "vertex set: ";
	print_vertices(*digraphPtr, name);
	cout << std::endl;

	cout << "edge set: ";
	print_edges(*digraphPtr, name);
	cout << std::endl;

	cout << "out-edges: " << std::endl;
	print_graph(*digraphPtr, name);
	cout << std::endl;
}

} /* namespace graph */
