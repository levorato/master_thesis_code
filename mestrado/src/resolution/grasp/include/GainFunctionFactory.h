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
#include "GainFunction.h"

namespace resolution {
namespace grasp {

class GainFunctionFactory {
public:
	ImbalanceGainFunction imbalanceFunction;
	ModularityGainFunction modularityFunction;
	NegativeModularityGainFunction negativeModularityGainFunction;

	GainFunctionFactory(SignedGraph* graph);
	virtual ~GainFunctionFactory();

	GainFunction& build(int functionType) {
		if(functionType == GainFunction::MODULARITY) {
			return modularityFunction;
		} else if(functionType == GainFunction::IMBALANCE) {
			return imbalanceFunction;
		} else {
			return negativeModularityGainFunction;
		}
	}
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTIONFACTORY_H_ */
