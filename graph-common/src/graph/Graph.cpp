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
#include <sstream>
#include <boost/log/trivial.hpp>

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

SignedGraph::SignedGraph(DirectedGraph &g, std::vector<long> subGraphNodeList) :
		graph(),
		n(subGraphNodeList.size()) {

	graph = g.create_subgraph(subGraphNodeList.begin(), subGraphNodeList.end());
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
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	DirectedGraph::edge_descriptor e;
	// O(n)
	DirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		e = *f;
		if(ew[e].weight < 0) {
			++sum;
		}
	}
	return sum;
}

unsigned long SignedGraph::getPositiveDegree(const unsigned long &a) {
	unsigned long sum = 0;
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	DirectedGraph::edge_descriptor e;
	// O(n)
	DirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		e = *f;
		if(ew[e].weight > 0) {
			++sum;
		}
	}
	return sum;
}

double SignedGraph::getNegativeEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray &cluster) {
	double sum = double(0.0);

	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	DirectedGraph::edge_descriptor e;
	// O(n)
	DirectedGraph::in_edge_iterator f, l;
	for (boost::tie(f, l) = in_edges(ni, graph); f != l; ++f) {
		e = *f;
		Vertex src = source(e, this->graph), targ = target(e, this->graph);
		long j = src.id;
		if((ew[e].weight < 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
		}
	}
	DirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(ni, graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, this->graph);
		if((ew[e].weight < 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
		}
	}
	return sum;
}

double SignedGraph::getPositiveEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray &cluster) {
	double sum = double(0.0);

	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	DirectedGraph::edge_descriptor e;
	// O(n)
	DirectedGraph::in_edge_iterator f, l;
	for (boost::tie(f, l) = in_edges(ni, graph); f != l; ++f) {
		e = *f;
		Vertex src = source(e, this->graph), targ = target(e, this->graph);
		long j = src.id;
		if((ew[e].weight > 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
		}
	}
	DirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(ni, graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, this->graph);
		if((ew[e].weight > 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
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

string SignedGraph::convertToGraclusInputFormat() {
	stringstream ss;
	int n = this->getN();
	// # of nodes and edges and format (format number 1 for integer-weighted edges)
	ss << n << " " << this->getM() << " 1\n";
	// nodes adjacent to vertex i and corresponding edge weight
	// vertex numbers start at 1
	// TODO: check if input graph can be directed!
	long edgeCount = 0;
	DirectedGraph::edge_descriptor e;
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	for(long i = 0; i < n; i++) {
		DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, this->graph); f != l; ++f) {
			e = *f;
			double weight = ew[e].weight;
			long j = target(*f, this->graph);
			// edge weights must be integer values!
			ss << (j+1) << " " << int(weight) << " ";
			edgeCount++;
		}
		// iterates over in-edges of vertex i
		DirectedGraph::in_edge_iterator in_i, in_end;
		// std::cout << "in-edges of " << i << ": ";
		for (tie(in_i, in_end) = in_edges(i, this->graph); in_i != in_end; ++in_i) {
			e = *in_i;
			Vertex src = source(e, this->graph), targ = target(e, this->graph);
			double weight = ew[e].weight;
			long j = src.id;
			// edge weights must be integer values!
			ss << (j+1) << " " << int(weight) << " ";
			edgeCount++;
		}
		ss << "\n";
	}
	BOOST_LOG_TRIVIAL(info) << "EdgeCount is " << edgeCount;
	return ss.str();
}

} /* namespace graph */
