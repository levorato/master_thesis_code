/*
 * ImbalanceGainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef IMBALANCEGAINFUNCTION_H_
#define IMBALANCEGAINFUNCTION_H_

#include "GainFunction.h"
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

namespace resolution {
namespace construction {

using namespace boost::numeric::ublas;

class ImbalanceGainFunction: public resolution::construction::GainFunction {
public:
	ImbalanceGainFunction(SignedGraph* g);
	virtual ~ImbalanceGainFunction();

	virtual GainCalculation calculateIndividualGain(ClusteringProblem& p,
			Clustering& c, int i);

	GainCalculation calculateIndividualGainCCProblem(
			ClusteringProblem& p, Clustering& c, int v);

	virtual void calculateGainList(ClusteringProblem &p, Clustering &c,
			list<GainCalculation>& nodeList);

	virtual int getType();

	virtual GainFunction::GainFunctionComparison getComparator();

private:
	matrix<double> h_VertexClusterPosSum;
	matrix<double> h_VertexClusterNegSum;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* IMBALANCEGAINFUNCTION_H_ */
