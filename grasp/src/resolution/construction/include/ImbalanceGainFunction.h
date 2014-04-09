/*
 * ImbalanceGainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef IMBALANCEGAINFUNCTION_H_
#define IMBALANCEGAINFUNCTION_H_

#include "GainFunction.h"

namespace resolution {
namespace construction {

class ImbalanceGainFunction: public resolution::construction::GainFunction {
public:
	ImbalanceGainFunction(SignedGraph* g);
	virtual ~ImbalanceGainFunction();

	virtual GainCalculation calculateIndividualGain(ClusteringProblem& p,
			Clustering& c, int i);

	virtual void calculateGainList(ClusteringProblem &p, Clustering &c,
			list<GainCalculation>& nodeList);

	virtual int getType();

	virtual GainFunction::GainFunctionComparison getComparator();

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* IMBALANCEGAINFUNCTION_H_ */
