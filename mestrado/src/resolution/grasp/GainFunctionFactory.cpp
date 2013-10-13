/*
 * GainFunctionFactory.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/GainFunctionFactory.h"

namespace resolution {
namespace grasp {

GainFunctionFactory::GainFunctionFactory(SignedGraph* g, const unsigned long& s) :
		imbalanceFunction(g, s), modularityFunction(g, s), negativeModularityGainFunction(g, s),
		positiveNegativeModularityGainFunction(g, s), positiveNegativeModularityGainFunctionII(g, s),
		positiveNegativeModularityGainFunctionIII(g, s) {
	// TODO Auto-generated constructor stub

}

GainFunctionFactory::~GainFunctionFactory() {
	// TODO Auto-generated destructor stub
}

} /* namespace grasp */
} /* namespace resolution */
