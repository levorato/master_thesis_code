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

ImbalanceGainFunction::ImbalanceGainFunction(SignedGraph* g) :
		GainFunction::GainFunction(g) {
	// TODO Auto-generated constructor stub

}

ImbalanceGainFunction::~ImbalanceGainFunction() {
	// TODO Auto-generated destructor stub
}

GainCalculation ImbalanceGainFunction::calculateIndividualGain(
		ClusteringProblem& p, Clustering& c, int i) {

	Clustering cTemp = c;
	GainCalculation gainCalculation;
	double min = std::numeric_limits<double>::max();
	gainCalculation.vertex = i;
	gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
	long k = 0;
	if (p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}

	// For each cluster k...
	int nc = c.getNumberOfClusters();
	for (unsigned long k = 0; k < nc; k++) {
		// cout << "Cluster " << k << endl;
		cTemp = c;
		cTemp.addNodeToCluster(*graph, p, i, k);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = k;
		}
	}
	// For a new cluster k+1
	// cout << "New cluster" << endl;
	if ((p.getType() == ClusteringProblem::CC_PROBLEM)
					or ((p.getType() == ClusteringProblem::RCC_PROBLEM)
							and (c.getNumberOfClusters() < k))) {
		cTemp = c;
		cTemp.addCluster(*graph, p, i);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
		}
	}
	gainCalculation.gainValue = min;
	// cout << "gain(a) = " << min << endl;
	return gainCalculation;
}

void ImbalanceGainFunction::calculateGainList(ClusteringProblem &p,
		Clustering &c, list<GainCalculation>& nodeList) {
	// cout << "Calculating gain list..." << endl;
	list<GainCalculation, allocator<GainCalculation> >::iterator pos;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if (p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}
	unsigned int i = 0;
	Clustering cTemp = c;
	for (i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		GainCalculation& gainCalculation = *pos;
		unsigned long a = gainCalculation.vertex;
		// cout << "Vertex " << a << endl;
		double min = std::numeric_limits<double>::max();
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

		// For each cluster k...
		int nc = c.getNumberOfClusters();
		for (unsigned long cl = 0; cl < nc; cl++) {
			// cout << "Cluster " << k << endl;
			cTemp = c;
			cTemp.addNodeToCluster(*graph, p, a, cl);
			Imbalance imb = cTemp.getImbalance();
			if (imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = cl;
			}
		}
		// For a new cluster k+1
		// cout << "New cluster" << endl;
		if ((p.getType() == ClusteringProblem::CC_PROBLEM)
				or ((p.getType() == ClusteringProblem::RCC_PROBLEM)
						and (c.getNumberOfClusters() < k))) {
			cTemp = c;
			cTemp.addCluster(*graph, p, a);
			Imbalance imb = cTemp.getImbalance();
			if (imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
			}
		}
		gainCalculation.gainValue = min;
		// cout << "gain(a) = " << min << endl;
	}
}

int ImbalanceGainFunction::getType() {
	return GainFunction::IMBALANCE;
}

GainFunction::GainFunctionComparison ImbalanceGainFunction::getComparator() {
	return GainFunctionComparison(true);
}

} /* namespace grasp */
} /* namespace resolution */
