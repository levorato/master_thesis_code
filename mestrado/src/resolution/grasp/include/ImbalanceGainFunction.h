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
namespace grasp {

class ImbalanceGainFunction: public resolution::grasp::GainFunction {
public:
	ImbalanceGainFunction(SignedGraph* g, const unsigned long& randomSeed);
	virtual ~ImbalanceGainFunction();

	virtual GainCalculation& gain(const int &a);

	virtual void calculateGainList(Clustering &c, GainFunctionVertexSet& nodeList);

	virtual bool operator () ( const int& a, const int& b );

	virtual int getType();

	virtual GainFunction::GainFunctionComparison getComparator();

	unsigned int chooseRandomNumber(int x);

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* IMBALANCEGAINFUNCTION_H_ */
