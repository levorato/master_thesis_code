/*
 * Graph.cpp
 *
 *  Created on: May 3, 2013
 *      Author: mario
 */

#include "include/Graph.h"
#include <boost/config.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <cassert>
#include <limits>

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

SignedGraph::SignedGraph(const unsigned long &numberOfNodes) :
		graph(numberOfNodes), n(numberOfNodes) {

}

SignedGraph::~SignedGraph() {
	// TODO Auto-generated destructor stub
}


unsigned long SignedGraph::getN() {
	return n;
}

unsigned long SignedGraph::getM() {
	return num_edges(graph);
}

unsigned int SignedGraph::getId() {
	return id;
}

void SignedGraph::setId(const unsigned int& i) {
	id = i;
}

void SignedGraph::addEdge(unsigned long a, unsigned long b, Edge edge) {
	add_edge(a, b, edge, graph);
}

unsigned long SignedGraph::getDegree(const unsigned long &a) {
	return out_degree(a, graph);
}

unsigned long SignedGraph::getNegativeDegree(const unsigned long &a) {
	unsigned long sum = 0;
	// O(n)
	DirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		if(((Edge*)f->get_property())->weight < 0) {
			++sum;
		}
	}
	return sum;
}

unsigned long SignedGraph::getPositiveDegree(const unsigned long &a) {
	unsigned long sum = 0;
	// O(n)
	DirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		if(((Edge*)f->get_property())->weight > 0) {
			++sum;
		}
	}
	return sum;
}

void SignedGraph::setGraphFileLocation(string txt) {
	graphFileLocation = txt;
}

string SignedGraph::getGraphFileLocation() {
	return graphFileLocation;
}

} /* namespace graph */
