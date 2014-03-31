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
namespace construction {

class NegativeModularityGainFunction: public resolution::construction::ModularityGainFunction {
public:
	NegativeModularityGainFunction(SignedGraph* g);
	virtual ~NegativeModularityGainFunction();

	virtual int getType();
	virtual void calculateModularityMatrix();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* NEGATIVEMODULARITYGAINFUNCTION_H_ */
