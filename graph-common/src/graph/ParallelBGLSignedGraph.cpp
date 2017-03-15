/*
 * ParallelBGLSignedGraph.cpp
 *
 *  Created on: 14 de mar de 2017
 *      Author: mlevorato
 */

#include "include/ParallelBGLSignedGraph.h"


namespace clusteringgraph {

ParallelBGLSignedGraph::~ParallelBGLSignedGraph() {
	// TODO Auto-generated destructor stub
}

// TODO implementar
ParallelBGLSignedGraph::ParallelBGLSignedGraph(const unsigned long &numberOfNodes) :
		graph(numberOfNodes), n(numberOfNodes), id(0) {

}

unsigned long ParallelBGLSignedGraph::getN() {
	return n;
}

unsigned long ParallelBGLSignedGraph::getM() {
	return num_edges(graph);
}

unsigned int ParallelBGLSignedGraph::getId() {
	return id;
}

void ParallelBGLSignedGraph::setId(const unsigned int& i) {
	id = i;
}

void ParallelBGLSignedGraph::addEdge(unsigned long a, unsigned long b, Edge edge) {
	add_edge(a, b, edge, graph);
}

// TODO REIMPLEMENTAR
unsigned long ParallelBGLSignedGraph::getDegree(const unsigned long &a) {
	return degree(a, graph);
}

// TODO REIMPLEMENTAR
unsigned long ParallelBGLSignedGraph::getOutDegree(const unsigned long &a) {
	return out_degree(a, graph);
}

// TODO REIMPLEMENTAR
unsigned long ParallelBGLSignedGraph::getNegativeDegree(const unsigned long &a) {
	unsigned long sum = 0;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	UndirectedGraph::edge_descriptor e;
	// O(n)
	UndirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		e = *f;
		if(ew[e].weight < 0) {
			++sum;
		}
	}
	return sum;
}

// TODO REIMPLEMENTAR
unsigned long ParallelBGLSignedGraph::getPositiveDegree(const unsigned long &a) {
	unsigned long sum = 0;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	UndirectedGraph::edge_descriptor e;
	// O(n)
	UndirectedGraph::out_edge_iterator f, l;
	for (boost::tie(f, l) = out_edges(a, graph); f != l; ++f) {
		e = *f;
		if(ew[e].weight > 0) {
			++sum;
		}
	}
	return sum;
}

// TODO REIMPLEMENTAR
double ParallelBGLSignedGraph::getNegativeEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray &cluster) {
	double sum = double(0.0);

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	UndirectedGraph::edge_descriptor e;
	// O(n)
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(ni, graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, this->graph);
		if((ew[e].weight < 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
		}
	}
	return sum;
}

// TODO REIMPLEMENTAR
double ParallelBGLSignedGraph::getPositiveEdgeSumBetweenVertexAndClustering(const unsigned long &ni, const ClusterArray &cluster) {
	double sum = double(0.0);

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
	UndirectedGraph::edge_descriptor e;
	// O(n)
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(ni, graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, this->graph);
		if((ew[e].weight > 0) and (cluster[j] >= 0)) {
			sum += ew[e].weight;
		}
	}
	return sum;
}

// TODO REIMPLEMENTAR
long ParallelBGLSignedGraph::getNumberOfEdgesInClustering(const ClusterArray& cluster, const long& clusterNumber) {
	long sum = 0;
	for(long ni = 0; ni < n; ni++) {  // O(n)
		boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, this->graph);
		UndirectedGraph::edge_descriptor e;
		// O(m)
		UndirectedGraph::out_edge_iterator f2, l2;
		for (boost::tie(f2, l2) = out_edges(ni, graph); f2 != l2; ++f2) {
			e = *f2;
			long j = target(*f2, this->graph);
			if((cluster[ni] >= 0) or (cluster[j] >= 0)) {
				sum++;
			}
		}
	}
	return sum;
}

void ParallelBGLSignedGraph::setGraphFileLocation(string txt) {
	graphFileLocation = txt;
}

string ParallelBGLSignedGraph::getGraphFileLocation() {
	return graphFileLocation;
}

} /* namespace clusteringgraph */
