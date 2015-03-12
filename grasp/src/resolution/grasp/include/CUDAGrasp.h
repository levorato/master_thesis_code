/*
 * CUDAGrasp.h
 *
 *  Created on: Feb 4, 2015
 *      Author: mario
 */

#ifndef CUDAGRASP_H_
#define CUDAGRASP_H_

#include "Grasp.h"

using namespace problem;
using namespace resolution::construction;

namespace resolution {
namespace grasp {

class CUDAGrasp : public resolution::grasp::Grasp {
public:
	CUDAGrasp();
	virtual ~CUDAGrasp();

	/**
	 * Triggers the execution of the GRASP algorithm using CUDA VND.
	 */
	Clustering executeGRASP(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
			SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info);

private:

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* CUDAGRASP_H_ */
