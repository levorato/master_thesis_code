/*
 * ModularityGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ModularityGainFunction.h"

namespace resolution {
namespace grasp {

ModularityGainFunction::ModularityGainFunction() {
	// TODO Auto-generated constructor stub

}

ModularityGainFunction::~ModularityGainFunction() {
	// TODO Auto-generated destructor stub
}

void ModularityGainFunction::calculateGainList(SignedGraph &g, Clustering &c,
		list<int>& nodeList) {
	gainMap.clear();
	ModularityMatrix& modularityMatrix = g.getModularityMatrix();
	int n = g.getN();
	list<int, allocator<int> >::const_iterator pos;
	// cout << "Calculating gain list..." << endl;
	// For each vertex a
	unsigned int i = 0;
	for(i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		int a = *pos;
		GainCalculation gainCalculation;
		float max = modularityMatrix[a][a];
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

		// For each cluster k...
		int nc = c.getNumberOfClusters();
		for(int k = 0; k < nc; k++) {
				int sum = 0;
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

} /* namespace grasp */
} /* namespace resolution */
