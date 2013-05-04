/*
 * Graph.cpp
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#include "Graph.h"
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

namespace graph {
typedef boost::graph_traits < DirectedGraph >::vertex_descriptor Vertex;

SignedGraph::SignedGraph(int numberOfNodes) : digraphPtr(new DirectedGraph(numberOfNodes)) {


}

SignedGraph::~SignedGraph() {
	// TODO Auto-generated destructor stub
}

void SignedGraph::addEdge(int a, int b, int value) {
	Vertex v0 = boost::add_vertex(a, *digraphPtr);
	Vertex v1 = boost::add_vertex(b, *digraphPtr);
	typedef typename DirectedGraph::edge_property_type Weight;
	add_edge(v0, v1, Weight(value), *digraphPtr);
}

void SignedGraph::printGraph() {
	const char* name = "ABCDEF";

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
