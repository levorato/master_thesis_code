/*
 * PositiveNegativeModularityGainFunctionIII.h
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#ifndef POSITIVENEGATIVEMODULARITYGAINFUNCTIONIII_H_
#define POSITIVENEGATIVEMODULARITYGAINFUNCTIONIII_H_

#include "GainFunction.h"
#include "ModularityGainFunction.h"


namespace resolution {
namespace grasp {

class PositiveNegativeModularityGainFunctionIII: public resolution::grasp::ModularityGainFunction {
public:
	PositiveNegativeModularityGainFunctionIII(SignedGraph* g);
	virtual ~PositiveNegativeModularityGainFunctionIII();

	virtual int getType();
	virtual void calculateModularityMatrix();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* POSITIVENEGATIVEMODULARITYGAINFUNCTIONIII_H_ */
