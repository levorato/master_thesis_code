/*
 * GainFunctionFactory.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/GainFunctionFactory.h"

namespace resolution {
namespace grasp {

GainFunctionFactory::GainFunctionFactory(SignedGraph* g) :
		imbalanceFunction(g), modularityFunction(g), negativeModularityGainFunction(g),
		positiveNegativeModularityGainFunction(g), positiveNegativeModularityGainFunctionII(g),
		positiveNegativeModularityGainFunctionIII(g) {
	// TODO Auto-generated constructor stub

}

GainFunctionFactory::~GainFunctionFactory() {
	// TODO Auto-generated destructor stub
}

} /* namespace grasp */
} /* namespace resolution */
