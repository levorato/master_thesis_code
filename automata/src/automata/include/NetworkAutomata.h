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
	 * Aplica a regra do jogo da vida considerando a particao do tabuleiro
	 * reservada ao respectivo processo. Nao aplica a regra nas fronteiras da particao.
	 */
	void applyRule(SignedGraph *g, int startIndex, int endIndex);
	/**
	 * Processa as regioes de intersecao entre uma divisao e outra do tabuleiro.
	 * O processamento deve ser feito sobre a geracao anterior.
	 */
	void processIntersection(SignedGraph *gdest, int line);
};

} /* namespace automata */

#endif /* NETWORKAUTOMATA_H_ */
