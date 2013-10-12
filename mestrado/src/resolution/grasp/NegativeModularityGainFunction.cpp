/*
 * NegativeModularityGainFunction.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#include "include/NegativeModularityGainFunction.h"

namespace resolution {
namespace grasp {

NegativeModularityGainFunction::NegativeModularityGainFunction(SignedGraph* g) :
		ModularityGainFunction::ModularityGainFunction(g) {
	// TODO Auto-generated constructor stub

}

NegativeModularityGainFunction::~NegativeModularityGainFunction() {
	// TODO Auto-generated destructor stub
}

// TODO calculate the modularity matrix of weighed graphs
void NegativeModularityGainFunction::calculateModularityMatrix() {
	int m = graph->getM();
	int numberOfNodes = graph->getN();
	double degree[numberOfNodes];
	// Prestore the degrees for optimezed lookup
	for(int i = 0; i < numberOfNodes; i++) {
		degree[i] = graph->getNegativeDegree(i);
	}
	/*
			 * TODO Alterar a maneira como as arestas sao varridas para esse calculo
	for(int i = 0; i < numberOfNodes; i++) {
		for(int j = 0; j < numberOfNodes; j++) {
			double a = (graph->getEdge(i, j) != 0) ? 1.0 : 0.0;
			modularityMatrix[i][j] = a - ( (degree[i] * degree[j]) / (2.0 * m) );
		}
	}
	*/
	modularityMatrixCalculated = true;
}

int NegativeModularityGainFunction::getType() {
	return GainFunction::NEGATIVE_MODULARITY;
}

} /* namespace grasp */
} /* namespace resolution */
