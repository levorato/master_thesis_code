/*
 * PositiveNegativeModularityGainFunction.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#include "include/PositiveNegativeModularityGainFunction.h"

namespace resolution {
namespace construction {

PositiveNegativeModularityGainFunction::PositiveNegativeModularityGainFunction(SignedGraph* g) :
				ModularityGainFunction::ModularityGainFunction(g) {
	// TODO Auto-generated constructor stub

}

PositiveNegativeModularityGainFunction::~PositiveNegativeModularityGainFunction() {
	// TODO Auto-generated destructor stub
}

// TODO calculate the modularity matrix of weighed graphs
void PositiveNegativeModularityGainFunction::calculateModularityMatrix() {
	int m = graph->getM();
	int numberOfNodes = graph->getN();
	double pos_degree[numberOfNodes], neg_degree[numberOfNodes];
	// Prestore the degrees for optimezed lookup
	for(int i = 0; i < numberOfNodes; i++) {
		neg_degree[i] = graph->getNegativeDegree(i);
		pos_degree[i] = graph->getPositiveDegree(i);
	}

	/*
	 * TODO Alterar a maneira como as arestas sao varridas para esse calculo
	for(int i = 0; i < numberOfNodes; i++) {
		for(int j = 0; j < numberOfNodes; j++) {
			double a = (graph->getEdge(i, j) != 0) ? 1.0 : 0.0;
			modularityMatrix[i][j] = a -
					( ( (neg_degree[i] * neg_degree[j]) + (pos_degree[i] * pos_degree[j]) ) / (2.0 * m) );
		}
	}
	*/
	modularityMatrixCalculated = true;
}

int PositiveNegativeModularityGainFunction::getType() {
	return GainFunction::POSITIVE_NEGATIVE_MODULARITY;
}

} /* namespace grasp */
} /* namespace resolution */
