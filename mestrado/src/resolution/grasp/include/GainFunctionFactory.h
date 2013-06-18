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
#include "GainFunction.h"

namespace resolution {
namespace grasp {

class GainFunctionFactory {
public:
	ImbalanceGainFunction imbalanceFunction;
	ModularityGainFunction modularityFunction;

	GainFunctionFactory();
	virtual ~GainFunctionFactory();

	GainFunction& build(int functionType) {
		if(functionType == GainFunction::MODULARITY) {
			return modularityFunction;
		} else {
			return imbalanceFunction;
		}
	}
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTIONFACTORY_H_ */
