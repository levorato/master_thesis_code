/*
 * ModularityGainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef MODULARITYGAINFUNCTION_H_
#define MODULARITYGAINFUNCTION_H_

#include "GainFunction.h"

namespace resolution {
namespace grasp {

class ModularityGainFunction: public resolution::grasp::GainFunction {
public:
	ModularityGainFunction();
	virtual ~ModularityGainFunction();

	virtual GainCalculation& gain(const int &a);

	/**
	 * Calculates the vertex gain list based on the modularity matrix
	 * of the graph.
	 */
	virtual void calculateGainList(SignedGraph &g, Clustering &c,
				list<int>& nodeList);

	virtual bool operator () ( const int& a, const int& b );

	virtual int getType();

	virtual GainFunction::GainFunctionComparison getComparator();

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* MODULARITYGAINFUNCTION_H_ */
