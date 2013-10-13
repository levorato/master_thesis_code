/*
 * PositiveNegativeModularityGainFunction.h
 *
 *  Created on: Jun 19, 2013
 *      Author: mario
 */

#ifndef POSITIVENEGATIVEMODULARITYGAINFUNCTION_H_
#define POSITIVENEGATIVEMODULARITYGAINFUNCTION_H_

#include "GainFunction.h"
#include "ModularityGainFunction.h"


namespace resolution {
namespace grasp {

class PositiveNegativeModularityGainFunction: public resolution::grasp::ModularityGainFunction {
public:
	PositiveNegativeModularityGainFunction(SignedGraph* g, const unsigned long& randomSeed);
	virtual ~PositiveNegativeModularityGainFunction();

	virtual int getType();
	virtual void calculateModularityMatrix();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* POSITIVENEGATIVEMODULARITYGAINFUNCTION_H_ */
