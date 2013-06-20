/*
 * ImbalanceGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ImbalanceGainFunction.h"

namespace resolution {
namespace grasp {

ImbalanceGainFunction::ImbalanceGainFunction(SignedGraph* g) : GainFunction::GainFunction(g) {
	// TODO Auto-generated constructor stub

}

ImbalanceGainFunction::~ImbalanceGainFunction() {
	// TODO Auto-generated destructor stub
}

void ImbalanceGainFunction::calculateGainList(Clustering &c, list<int>& nodeList) {
	gainMap.clear();
	list<int, allocator<int> >::const_iterator pos;
	// cout << "Calculating gain list..." << endl;
	unsigned int i = 0;
	for(i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		int a = *pos;
		// cout << "Vertex " << a << endl;
		GainCalculation gainCalculation;
		double min = std::numeric_limits<double>::max();
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

		// For each cluster k...
		int nc = c.getNumberOfClusters();
		for(int k = 0; k < nc; k++) {
			// cout << "Cluster " << k << endl;
			Imbalance delta = c.calculateDeltaObjectiveFunction(*graph, c.getCluster(k), a);
			if(delta.getValue() < min) {
				min = delta.getValue();
				gainCalculation.clusterNumber = k;
			}
		}
		// For a new cluster k+1
		// cout << "New cluster" << endl;
		BoolArray newCluster(MAX_NODES);
		newCluster[a] = true;
		Imbalance delta = c.calculateDeltaObjectiveFunction(*graph, newCluster, a);
		if(delta.getValue() < min) {
			min = delta.getValue();
			gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
		}
		gainCalculation.value = min;
		// cout << "gain(a) = " << min << endl;
		gainMap[a] = gainCalculation;
	}
}

/**
 * TODO For a given vertex a, calculates the minimum value of imbalance (I(P))
 * of inserting 'a' into a new or an existing clustering k. Returns the minimum imbalance
 * and the cluster corresponding to it.
 */
GainCalculation& ImbalanceGainFunction::gain(const int &a) {
	return gainMap[a];
}

bool ImbalanceGainFunction::operator () ( const int& a, const int& b ) {
	return this->gain(a).value < this->gain(b).value;
}

int ImbalanceGainFunction::getType() {
	return GainFunction::IMBALANCE;
}

GainFunction::GainFunctionComparison ImbalanceGainFunction::getComparator() {
	return GainFunctionComparison(this, true);
}

} /* namespace grasp */
} /* namespace resolution */
