/*
 * NetworkAutomata.cpp
 *
 *  Created on: 21/05/2014
 *      Author: mlevorato
 */

#include "include/NetworkAutomata.h"
#include <boost/log/trivial.hpp>

namespace automata {

using namespace boost;

NetworkAutomata::NetworkAutomata(SignedGraph *g) : graph(g) {

}

NetworkAutomata::~NetworkAutomata() {
	// TODO Auto-generated destructor stub
}

void NetworkAutomata::execute(int numberOfGenerations) {

	for(int i = 0; i < numberOfGenerations; i++) {
		BOOST_LOG_TRIVIAL(info) << "Generation " << (i+1) << endl;
		int n = graph->getN();
		double edgeSum = 0.0;

		// traverses all edges of the (symmetric) graph
		for(int i = 0; i < n; i++) {  // For each vertex i
			DirectedGraph::out_edge_iterator f, l;
			// For each out edge of i => (i, j)
			for (boost::tie(f, l) = out_edges(i, graph->graph); f != l; ++f) {
				Edge *e = (Edge*)f->get_property();
				double weight = e->weight;
				int j = target(*f, graph->graph);
				edgeSum += fabs(weight);
				// apply the rule to the edge (i, j)
				double newWeight = applyRule(graph, i, j, weight);
				e->weight = newWeight;
			}
		}
		BOOST_LOG_TRIVIAL(info) << "Current edge sum is: " << edgeSum;
	}
}

void NetworkAutomata::print() {
	// TODO implement me
}

double NetworkAutomata::applyRule(SignedGraph *g, int i, int j, double w) {
	// TODO implementar logica segundo artigo 'Study of sign adjustment in weighted signed networks'
	// Local pressure adjustment rule
	//  (1) Select at random any three nodes from the network and, if they are each linked with the others,
	//        in a 3-cycle, go step 2. Otherwise, select another three nodes randomly, until a completed
	//        3-cycle is selected.
	//  (2) Select at random one edge from the selected 3-cycle to adjust.
	//  (3) The selected edge is adjusted according to the equation:
	//        w(i,j)[t+] = w(i,j)[t] + w(i, k)[t] x w(k, j)[t] x (1 - abs(w(i, j)[t]))

	return w*2.0;
}

} /* namespace automata */
