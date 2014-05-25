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
		double changes = 0.0;

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
				changes += applyRule(graph, i, j, weight);
				// e->weight = newWeight;
			}
		}
		BOOST_LOG_TRIVIAL(info) << "Current edge sum is: " << edgeSum;
		BOOST_LOG_TRIVIAL(info) << "Number of edge value changes: " << changes;
	}
}

void NetworkAutomata::print() {
	// TODO implement me
}

double NetworkAutomata::applyRule(SignedGraph *g, unsigned long i, unsigned long j, double w) {
	// TODO implementar logica segundo artigo 'Study of sign adjustment in weighted signed networks'
	// Local pressure adjustment rule
	//  (1) Select at random any three nodes from the network and, if they are each linked with the others,
	//        in a 3-cycle, go step 2. Otherwise, select another three nodes randomly, until a completed
	//        3-cycle is selected.
	//  (2) Select at random one edge from the selected 3-cycle to adjust.
	//  (3) The selected edge is adjusted according to the equation:
	//        w(i,j)[t+] = w(i,j)[t] + w(i, k)[t] x w(k, j)[t] x (1 - abs(w(i, j)[t]))

	// check if this edge belongs to a 3-cycle, that is, if there is another vertex k
	// such that (i, k) and (j, k) exists (supposing a symmetric graph)
	// BOOST_LOG_TRIVIAL(info) << "Applying rule to edge (" << i << ", " << j << ")" << endl;
	unsigned long n = g->getN();
	double numChanges = 0.0;
	for(unsigned long k = 0; k < n; k++) {
		DirectedGraph::out_edge_iterator f, l;
		bool hit_i = false;
		bool hit_j = false;
		// For each out edge of k => (k, x)
		for (boost::tie(f, l) = out_edges(k, graph->graph); f != l; ++f) {
			Edge *e = (Edge*)f->get_property();
			double weight = e->weight;
			unsigned long x = target(*f, graph->graph);
			if(x == i) {  hit_i = true;  }
			if(x == j) {  hit_j = true;  }
			if(hit_i and hit_j)  break;
		}
		if(hit_i and hit_j) {  // i, j and k form a 3-cycle
			// BOOST_LOG_TRIVIAL(info) << "Vertex " << k << " forms a 3-cycle with i and j.";
			// select at random one edge of the 3-cycle to adjust
			// TODO: this is going to cause problems in the parallel version
			// TODO: check what is the best option for directed graphs... test all combinations of i, j and k?
			std::vector< std::pair<unsigned long, unsigned long> > edgeList;
			edgeList.push_back(make_pair(i, j));
			edgeList.push_back(make_pair(k, i));
			edgeList.push_back(make_pair(k, j));
			std::random_shuffle(edgeList.begin(), edgeList.end());
			// TODO shuffle not working
			std::pair<unsigned long, unsigned long> edge = edgeList[0];
			// retrieves the the selected edge
			DirectedGraph::out_edge_iterator f2, l2;
			for (boost::tie(f2, l2) = out_edges(edge.first, graph->graph); f2 != l2 && target(*f2, graph->graph) != edge.second; ++f2);
			Edge *e2 = (Edge*)f2->get_property();
			// BOOST_LOG_TRIVIAL(info) << "Random edge from 3-cycle is (" << edge.first << ", " << edge.second << ") = " << e2->weight;
			// retrieves the weight of the other 2 edges of the cycle
			DirectedGraph::out_edge_iterator f3, l3;
			for (boost::tie(f3, l3) = out_edges(edgeList[1].first, graph->graph); f3 != l3 && target(*f3, graph->graph) != edgeList[1].second; ++f3);
			double w_ik = ((Edge*)f3->get_property())->weight;
			for (boost::tie(f3, l3) = out_edges(edgeList[2].first, graph->graph); f3 != l3 && target(*f3, graph->graph) != edgeList[2].second; ++f3);
			double w_kj = ((Edge*)f3->get_property())->weight;
			// BOOST_LOG_TRIVIAL(info) << "Operands: " << w_ik << " and " << w_kj;
			// update the weight of the selected edge
			double currentW = e2->weight;
			//        w(i,j)[t+] = w(i,j)[t] + w(i, k)[t] x w(k, j)[t] x (1 - abs(w(i, j)[t]))
			double newW = currentW + (1.0 - fabs(currentW)) * w_ik * w_kj;
			if(newW != currentW) { numChanges++; }
			// BOOST_LOG_TRIVIAL(info) << "New weight for edge (" << edge.first << ", " << edge.second << ") = " << newW;
			e2->weight = newW;
			break;
		}
	}

	return numChanges;
}

} /* namespace automata */
