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

using namespace boost;

typedef adjacency_matrix<directedS, no_property, property<edge_weight_t, int > > DirectedGraph;
typedef boost::scoped_ptr<DirectedGraph> DigraphPtr;


namespace clusteringgraph {
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
