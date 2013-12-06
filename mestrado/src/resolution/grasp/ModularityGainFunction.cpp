/*
 * ModularityGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ModularityGainFunction.h"

namespace resolution {
namespace grasp {

ModularityGainFunction::ModularityGainFunction(SignedGraph* g) :
		GainFunction::GainFunction(g),
		modularityMatrixCalculated(false),
		modularityMatrix(boost::extents[g->getN()][g->getN()]) {
	// TODO Auto-generated constructor stub

}

ModularityGainFunction::~ModularityGainFunction() {
	// TODO Auto-generated destructor stub
}

// TODO calculate the modularity matrix of weighed graphs
void ModularityGainFunction::calculateModularityMatrix() {
	int m = graph->getM();
	int numberOfNodes = graph->getN();
	double degree[numberOfNodes];
	// Prestore the degrees for optimezed lookup
	for(int i = 0; i < numberOfNodes; i++) {
		degree[i] = graph->getDegree(i);
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

ModularityMatrix& ModularityGainFunction::getModularityMatrix() {
	if(not modularityMatrixCalculated) {
		calculateModularityMatrix();
	}
	return modularityMatrix;
}

void ModularityGainFunction::calculateGainList(ClusteringProblem &p, Clustering &c,
		GainFunctionVertexSet& nodeList) {
	gainMap.clear();
	ModularityMatrix& modularityMatrix = getModularityMatrix();
	int n = graph->getN();
	// cout << "Calculating gain list..." << endl;
	// For each vertex a
	list<int, allocator<int> >::const_iterator pos;
	// cout << "Calculating gain list..." << endl;
	unsigned int i = 0;
	for(i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		int a = *pos;
		GainCalculation gainCalculation;
		double max = modularityMatrix[a][a];
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

		// For each cluster k...
		int nc = c.getNumberOfClusters();
		for(unsigned long k = 0; k < nc; k++) {
				double sum = 0.0;
				// Cluster(k)
				BoolArray cluster = c.getCluster(k);
				// j in Cluster(k)
				for(int j = 0; j < n; j++) {
						if(cluster[j]) {
								sum += 2 * modularityMatrix[a][j];
						}
				}
				sum += modularityMatrix[a][a];
				if(sum > max) {
						max = sum;
						gainCalculation.clusterNumber = k;
				}
		}
		gainCalculation.value = max;
		gainMap[a] = gainCalculation;
	}
}

/**
 * TODO For a given vertex a, calculates the minimum value of imbalance (I(P))
 * of inserting 'a' into a new or an existing clustering k. Returns the minimum imbalance
 * and the cluster corresponding to it.
 */
GainCalculation& ModularityGainFunction::gain(const int &a) {
	return gainMap[a];
}

bool ModularityGainFunction::operator () ( const int& a, const int& b ) {
	return gain(a).value > gain(b).value;
}

int ModularityGainFunction::getType() {
	return GainFunction::MODULARITY;
}

GainFunction::GainFunctionComparison ModularityGainFunction::getComparator() {
	return GainFunctionComparison(this, false);
}

} /* namespace grasp */
} /* namespace resolution */
