/*
 * ParallelGrasp.cpp
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#include "include/ParallelGrasp.h"

namespace resolution {
namespace grasp {

ParallelGrasp::ParallelGrasp() {
	// TODO Auto-generated constructor stub

}

ParallelGrasp::~ParallelGrasp() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr ParallelGrasp::executeGRASP(SignedGraph *g, int iter, float alpha, int l,
			ClusteringProblem& problem, std::ostream& os) {
	// TODO implement this method with MPI
	int np = 5;
	return Grasp::executeGRASP(g, iter / np, alpha, l, problem, os);
}

} /* namespace grasp */
} /* namespace resolution */
