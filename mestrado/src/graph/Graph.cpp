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

SignedGraph::SignedGraph(const int &numberOfNodes) : graph(numberOfNodes),
		n(numberOfNodes) {

}

SignedGraph::~SignedGraph() {
	// TODO Auto-generated destructor stub
}


int SignedGraph::getN() {
	return n;
}

int SignedGraph::getM() {
	return graph.m_num_edges;
}

void SignedGraph::addEdge(int a, int b, Edge edge) {
	add_edge(a, b, edge, graph);
}

double SignedGraph::getEdge(const int &a, const int &b) {
	return graph.get_edge(a, b).second.weight;
}

bool SignedGraph::isPositiveEdge(const int &a, const int &b) {
	double value = getEdge(a, b);
	return (value > 2 * std::numeric_limits<double>::epsilon());
}

bool SignedGraph::isNegativeEdge(const int &a, const int &b) {
	 double value = getEdge(a, b);
        return (value < (-2) * std::numeric_limits<double>::epsilon());
}

int SignedGraph::getDegree(const int &a) {
	return in_degree(a, graph);
}

int SignedGraph::getNegativeDegree(const int &a) {
	int sum = 0;
	// O(n)
	typename adjacency_matrix<directedS, no_property, Edge >::in_edge_iterator f, l;
	for (boost::tie(f, l) = in_edges(a, graph); f != l; ++f) {
		if(((Edge*)f->get_property())->weight < 0) {
			++sum;
		}
	}
	// cout << "Negative degree of vertex " << a << " is " << sum << endl;
	// cout << "Positive degree of vertex " << a << " is " << getDegree(a) << endl;
	return sum;
}

void SignedGraph::printGraph() {
	for(int i = 0; i < this->getN(); i++) {
		for(int j = 0; j < this->getN(); j++) {
			cout << this->getEdge(i, j) << "  ";
		}
		cout << endl;
	}
}

void SignedGraph::setGraphAsText(string txt) {
	graphAsText = txt;
}

string SignedGraph::getGraphAsText() {
	return graphAsText;
}

} /* namespace graph */
