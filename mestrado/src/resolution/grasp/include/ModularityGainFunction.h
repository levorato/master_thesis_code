/*
 * ModularityGainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef MODULARITYGAINFUNCTION_H_
#define MODULARITYGAINFUNCTION_H_

#include "GainFunction.h"

#include <boost/multi_array.hpp>

namespace resolution {
namespace grasp {

// the modularity matrix: a matrix of double
typedef multi_array<double, 2> ModularityMatrix;

class ModularityGainFunction: public resolution::grasp::GainFunction {
public:
	ModularityGainFunction(SignedGraph* g);
	virtual ~ModularityGainFunction();

	virtual GainCalculation& gain(const int &a);

	/**
	 * Calculates the vertex gain list based on the modularity matrix
	 * of the graph.
	 */
	virtual void calculateGainList(ClusteringProblem &p, Clustering &c,
			GainFunctionVertexSet& nodeList);

	virtual bool operator () ( const int& a, const int& b );

	virtual int getType();

	virtual GainFunction::GainFunctionComparison getComparator();

protected:
	bool modularityMatrixCalculated;
	ModularityMatrix modularityMatrix;

	virtual void calculateModularityMatrix();
	ModularityMatrix& getModularityMatrix();

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* MODULARITYGAINFUNCTION_H_ */
