/*
 * ImbalanceGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ImbalanceGainFunction.h"
#include "../../problem/include/CCProblem.h"
#include "../../problem/include/RCCProblem.h"

namespace resolution {
namespace construction {

ImbalanceGainFunction::ImbalanceGainFunction(SignedGraph* g) : GainFunction::GainFunction(g) {
	// TODO Auto-generated constructor stub

}

ImbalanceGainFunction::~ImbalanceGainFunction() {
	// TODO Auto-generated destructor stub
}

void ImbalanceGainFunction::calculateGainList(ClusteringProblem &p, Clustering &c,
		GainFunctionVertexSet& nodeList) {
	gainMap.clear();
	// cout << "Calculating gain list..." << endl;
	list<int, allocator<int> >::const_iterator pos;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}
	unsigned int i = 0;
	Clustering cTemp = c;
	for(i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
	    unsigned long a = *pos;
		// cout << "Vertex " << a << endl;
		GainCalculation gainCalculation;
		double min = std::numeric_limits<double>::max();
		gainCalculation.clusterNumber = 0;

		// For each cluster k...
		int nc = c.getNumberOfClusters();
		for(unsigned long cl = 0; cl < nc; cl++) {
			// cout << "Cluster " << k << endl;
			cTemp = c;
			cTemp.addNodeToCluster(*graph, p, a, cl);
			Imbalance imb = cTemp.getImbalance();
			if(imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = cl;
			}
		}
		// For a new cluster k+1
		// cout << "New cluster" << endl;
		if((p.getType() == ClusteringProblem::CC_PROBLEM) or ((p.getType() == ClusteringProblem::RCC_PROBLEM)
				and (c.getNumberOfClusters() < k) ) ) {
			cTemp = c;
			cTemp.addCluster(*graph, p, a);
			Imbalance imb = cTemp.getImbalance();
			if(imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
			}
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
