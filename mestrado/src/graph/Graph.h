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

using namespace boost;

namespace clusteringgraph {

struct Edge {
    int weight;
    Edge() : weight(0) { }
    Edge(int w) : weight(w) { }
};
typedef adjacency_matrix<directedS, no_property, Edge > DirectedGraph;
typedef boost::scoped_ptr<DirectedGraph> DigraphPtr;

class SignedGraph {
public:
	SignedGraph(int numberOfNodes);
	virtual ~SignedGraph();

	void addEdge(int a, int b, int value);
	void printGraph();

private:
	DigraphPtr digraphPtr;

};

} /* namespace clusteringgraph */
#endif /* GRAPH_H_ */
