/*
 * CUDAImbalanceGainFunction.cpp
 *
 *  Created on: 09/02/2015
 *      Author: czt0
 */

#include "include/CUDAImbalanceGainFunction.h"
#include "../vnd/include/CUDASearch.h"

namespace resolution {
namespace construction {

CUDAImbalanceGainFunction::CUDAImbalanceGainFunction(SignedGraph *g) :
		ImbalanceGainFunction(g) {

}

CUDAImbalanceGainFunction::~CUDAImbalanceGainFunction() {

}

GainCalculation ImbalanceGainFunction::calculateIndividualGain(
		ClusteringProblem& p, Clustering& c, int i,
		thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
		thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset) {

	GainCalculation gainCalculationCUDA;
	gainCalculationCUDA.vertex = i;
	gainCalculationCUDA.clusterNumber = Clustering::NEW_CLUSTER;

	// TODO transform into class constant
	// number of threads per block
	unsigned short threadsCount = 256;  // limited by shared memory size

	// Pass raw array and its size to kernel
	ulong clusterNumber = -1;
	double gainValue = 0.0;
	// objective function value
	thrust::host_vector<float> h_functionValue(1);
	h_functionValue[0] = c.getImbalance().getValue();
	thrust::host_vector<unsigned long> h_mycluster(c.getClusterArray());
	runConstructKernel(h_weights, h_dest, h_numedges, h_offset, h_mycluster, graph->getN(),
			graph->getM(), c.getNumberOfClusters(), threadsCount, i, clusterNumber, gainValue);

	gainCalculationCUDA.vertex = i;
	gainCalculationCUDA.clusterNumber = clusterNumber;
	gainCalculationCUDA.gainValue = gainValue;
	// return gainCalculationCUDA;

	// VALIDANDO o calculo acima com o calculo tradicional na CPU
	double min = std::numeric_limits<double>::max();
	Clustering cTemp = c;
	GainCalculation gainCalculation;
	// For each cluster ci...
	int nc = c.getNumberOfClusters();
	for (unsigned long ci = 0; ci < nc; ci++) {
		cTemp = c;
		cTemp.addNodeToCluster(*graph, p, i, ci);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = ci;
		}
	}
	// For a new cluster ci+1
	if ((p.getType() == ClusteringProblem::CC_PROBLEM))) {
		cTemp = c;
		cTemp.addCluster(*graph, p, i);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
		}
	}
	gainCalculation.gainValue = min;

	assert(gainCalculationCUDA.gainValue == gainCalculation.gainValue);
	assert(gainCalculationCUDA.clusterNumber == gainCalculation.clusterNumber);
}

void ImbalanceGainFunction::calculateGainList(ClusteringProblem &p,
		Clustering &c, list<GainCalculation>& nodeList, thrust::host_vector<float>& h_weights,
		thrust::host_vector<int>& h_dest, thrust::host_vector<int>& h_numedges,
		thrust::host_vector<int>& h_offset) {
	// cout << "Calculating gain list..." << endl;
	list<GainCalculation, allocator<GainCalculation> >::iterator pos;
	unsigned int i = 0;
	Clustering cTemp = c;
	for (i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		GainCalculation& gainCalculation = *pos;
		unsigned long a = gainCalculation.vertex;
		GainCalculation gain = calculateIndividualGain(p, c, a, h_weights, h_dest, h_numedges, h_offset);
		gainCalculation.clusterNumber = gain.clusterNumber;
		gainCalculation.gainValue = gain.gainValue;
	}
}

} /* namespace construction */
} /* namespace resolution */
