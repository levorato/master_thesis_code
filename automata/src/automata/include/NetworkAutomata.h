/*
 * NetworkAutomata.h
 *
 *  Created on: 21/05/2014
 *      Author: mlevorato
 */

#ifndef NETWORKAUTOMATA_H_
#define NETWORKAUTOMATA_H_

#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"

namespace automata {

using namespace boost;

class NetworkAutomata {
public:
	NetworkAutomata(SignedGraph *g);
	virtual ~NetworkAutomata();

	void execute(int numberOfGenerations);
	void print();

private:
	SignedGraph *graph;
	/**
	 * Aplica a regra do automato a uma aresta do grafo.
	 */
	double applyRule(SignedGraph *g, int i, int j, double w);
};

} /* namespace automata */

#endif /* NETWORKAUTOMATA_H_ */
