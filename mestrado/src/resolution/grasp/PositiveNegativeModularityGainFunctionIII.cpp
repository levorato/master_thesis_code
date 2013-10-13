/*
 * PositiveNegativeModularityGainFunctionIII.cpp
 *
 *  Created on: Jul 18, 2013
 *      Author: mario
 */

#include "include/PositiveNegativeModularityGainFunctionIII.h"

namespace resolution {
namespace grasp {

PositiveNegativeModularityGainFunctionIII::PositiveNegativeModularityGainFunctionIII(SignedGraph* g, const unsigned long& s) :
				ModularityGainFunction::ModularityGainFunction(g, s) {
	// TODO Auto-generated constructor stub

}

PositiveNegativeModularityGainFunctionIII::~PositiveNegativeModularityGainFunctionIII() {
	// TODO Auto-generated destructor stub
}

// TODO calculate the modularity matrix of weighed graphs
void PositiveNegativeModularityGainFunctionIII::calculateModularityMatrix() {
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
			double a = 0.0;
			if(graph->getEdge(i, j) > 0) {
				a = 1.0;
				modularityMatrix[i][j] = a + (neg_degree[i] * neg_degree[j]) / (2.0 * m);
			} else if(graph->getEdge(i, j) < 0) {
				a = -1.0;
				modularityMatrix[i][j] = a - (pos_degree[i] * pos_degree[j]) / (2.0 * m);
			} else {  // a(i,j) == zero
				modularityMatrix[i][j] = a +
									( ( (neg_degree[i] * neg_degree[j]) - (pos_degree[i] * pos_degree[j]) ) / (2.0 * m) );
			}
		}
	}
	*/
	modularityMatrixCalculated = true;
}

int PositiveNegativeModularityGainFunctionIII::getType() {
	return GainFunction::POSITIVE_NEGATIVE_MODULARITY_III;
}

} /* namespace grasp */
} /* namespace resolution */
