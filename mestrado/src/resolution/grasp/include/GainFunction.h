/*
 * GainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef GAINFUNCTION_H_
#define GAINFUNCTION_H_

#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"

using namespace clusteringgraph;

namespace resolution {
namespace grasp {

typedef struct {
	double value;
	int clusterNumber;
} GainCalculation;

class GainFunction {
public:
	static const int MODULARITY = 0, IMBALANCE = 1;

	GainFunction();
	virtual ~GainFunction();

	class GainFunctionComparison
	{
	        GainFunction* function;
	public:
	  GainFunctionComparison(GainFunction *f) : function(f)
	  {   }

	    bool operator () ( const int& a, const int& b ) const
	    {
	      return function->gain(a).value < function->gain(b).value;
	    }
	};

	/**
	 * Returns the gain of a given vertex
	 * and the number of the cluster where the insertion of the vertex
	 * brings the best gain possible (return type is GainCalculation).
	 */
	virtual GainCalculation& gain(const int &a) = 0;

	virtual void calculateGainList(SignedGraph &g, Clustering &c,
			list<int>& nodeList) = 0;

	GainFunctionComparison getComparator() {
		return GainFunctionComparison(this);
	}

	virtual int getType() = 0;

protected:
	/** the map of nodes' gain value */
	map<int, GainCalculation> gainMap;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTION_H_ */
