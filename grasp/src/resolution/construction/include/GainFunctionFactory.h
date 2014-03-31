/*
 * GainFunctionFactory.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef GAINFUNCTIONFACTORY_H_
#define GAINFUNCTIONFACTORY_H_

#include "ImbalanceGainFunction.h"
#include "ModularityGainFunction.h"
#include "NegativeModularityGainFunction.h"
#include "PositiveNegativeModularityGainFunction.h"
#include "PositiveNegativeModularityGainFunctionII.h"
#include "PositiveNegativeModularityGainFunctionIII.h"
#include "GainFunction.h"

namespace resolution {
namespace construction {

class GainFunctionFactory {
public:
	ImbalanceGainFunction imbalanceFunction;
	ModularityGainFunction modularityFunction;
	NegativeModularityGainFunction negativeModularityGainFunction;
	PositiveNegativeModularityGainFunction positiveNegativeModularityGainFunction;
	PositiveNegativeModularityGainFunctionII positiveNegativeModularityGainFunctionII;
	PositiveNegativeModularityGainFunctionIII positiveNegativeModularityGainFunctionIII;

	GainFunctionFactory(SignedGraph* graph);
	virtual ~GainFunctionFactory();

	GainFunction& build(int functionType) {
		if(functionType == GainFunction::MODULARITY) {
			return modularityFunction;
		} else if(functionType == GainFunction::IMBALANCE) {
			return imbalanceFunction;
		} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY) {
			return positiveNegativeModularityGainFunction;
		} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY_II) {
			return positiveNegativeModularityGainFunctionII;
		} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY_III) {
			return positiveNegativeModularityGainFunctionIII;
		} else {
			return negativeModularityGainFunction;
		}
	}
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTIONFACTORY_H_ */
