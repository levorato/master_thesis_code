/*
 * NetworkAutomata.cpp
 *
 *  Created on: 21/05/2014
 *      Author: mlevorato
 */

#include "include/NetworkAutomata.h"

namespace automata {

NetworkAutomata::NetworkAutomata(SignedGraph *g) : graph(g) {

}

NetworkAutomata::~NetworkAutomata() {
	// TODO Auto-generated destructor stub
}

void NetworkAutomata::execute(int numberOfGenerations) {
	// TODO implement me
	for(int i = 0; i < numberOfGenerations; i++) {
		/*
		for(int j = 0; j < height; j++) {
			for(int k = 0; k < width; k++) {
				int n = numberOfNeighbors(j, k, height, width);
				if(n == 3) {
					(*newmatrixPtr)[j][k] = 1;
				} else if (n == 2) {
					(*newmatrixPtr)[j][k] = (*matrixPtr)[j][k];
				} else {
					(*newmatrixPtr)[j][k] = 0;
				}
			}
		}
		*/
		cout << "Generation " << (i+1) << endl;
		// print(matrixPtr.get(), height, width);
		// matrixPtr.reset(newmatrixPtr);
		// newmatrixPtr = new Matrix(boost::extents[height][width]);
	}
}

void NetworkAutomata::print() {
	// TODO implement me
}

void NetworkAutomata::applyRule(SignedGraph *g) {
	int h = 0; // dest->shape()[0];
	int w = 0; // dest->shape()[1];
	cout << "Aplicando a regra para as linhas " << startIndex << " e " << endIndex << endl;

	for(int i = startIndex; i < endIndex; i++) {
		for(int j = 0; j < w; j++) {
			/*
			int n = numberOfNeighbors(i, j, h, w);
			if(n == 3) {
				(*dest)[i][j] = 1;
			} else if (n == 2) {
				(*dest)[i][j] = (*matrixPtr)[i][j];
			} else {
				(*dest)[i][j] = 0;
			}
			*/
		}
	}
}

} /* namespace automata */
