/*
 * CUDAImbalanceGainFunction.h
 *
 *  Created on: 09/02/2015
 *      Author: czt0
 */

#ifndef CUDAIMBALANCEGAINFUNCTION_H_
#define CUDAIMBALANCEGAINFUNCTION_H_

#include "ImbalanceGainFunction.h"

#include <thrust/host_vector.h>

using namespace thrust;

namespace resolution {
namespace construction {

class CUDAImbalanceGainFunction : public ImbalanceGainFunction {
public:
	CUDAImbalanceGainFunction(SignedGraph *g);
	virtual ~CUDAImbalanceGainFunction();

	virtual GainCalculation calculateIndividualGain(ClusteringProblem& p,
			Clustering& c, int i, thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset);

	virtual void calculateGainList(ClusteringProblem &p, Clustering &c,
			list<GainCalculation>& nodeList, thrust::host_vector<float>& h_weights, thrust::host_vector<int>& h_dest,
			thrust::host_vector<int>& h_numedges, thrust::host_vector<int>& h_offset);
};

} /* namespace construction */
} /* namespace resolution */
#endif /* CUDAIMBALANCEGAINFUNCTION_H_ */
