/*
 * NegativeModularityGainFunction.h
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#ifndef NEGATIVEMODULARITYGAINFUNCTION_H_
#define NEGATIVEMODULARITYGAINFUNCTION_H_

#include "GainFunction.h"
#include "ModularityGainFunction.h"

namespace resolution {
namespace grasp {

class NegativeModularityGainFunction: public resolution::grasp::ModularityGainFunction {
public:
	NegativeModularityGainFunction(SignedGraph* g, const unsigned long& randomSeed);
	virtual ~NegativeModularityGainFunction();

	virtual int getType();
	virtual void calculateModularityMatrix();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* NEGATIVEMODULARITYGAINFUNCTION_H_ */
