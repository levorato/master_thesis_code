/*
 * PositiveNegativeModularityGainFunctionII.h
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#ifndef POSITIVENEGATIVEMODULARITYGAINFUNCTIONII_H_
#define POSITIVENEGATIVEMODULARITYGAINFUNCTIONII_H_

#include "GainFunction.h"
#include "ModularityGainFunction.h"


namespace resolution {
namespace grasp {

class PositiveNegativeModularityGainFunctionII: public resolution::grasp::ModularityGainFunction {
public:
	PositiveNegativeModularityGainFunctionII(SignedGraph* g);
	virtual ~PositiveNegativeModularityGainFunctionII();

	virtual int getType();
	virtual void calculateModularityMatrix();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* POSITIVENEGATIVEMODULARITYGAINFUNCTIONII_H_ */
