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
	static const int IMBALANCE = 0, MODULARITY = 1, NEGATIVE_MODULARITY = 2;
	static const int POSITIVE_NEGATIVE_MODULARITY = 3;

	GainFunction(SignedGraph* g);
	virtual ~GainFunction();

	class GainFunctionComparison
	{
	private:
	        GainFunction* function;
	        bool ascendingOrder;
	public:
			GainFunctionComparison(GainFunction *f, bool ascending) :
			  function(f), ascendingOrder(ascending)
			{   }

			bool operator () ( const int& a, const int& b ) const {
				if(ascendingOrder) {
					return function->gain(a).value < function->gain(b).value;
				} else {
					return function->gain(a).value > function->gain(b).value;
				}
			}
	};

	/**
	 * Returns the gain of a given vertex
	 * and the number of the cluster where the insertion of the vertex
	 * brings the best gain possible (return type is GainCalculation).
	 */
	virtual GainCalculation& gain(const int &a) = 0;

	virtual void calculateGainList(Clustering &c, list<int>& nodeList) = 0;

	virtual GainFunctionComparison getComparator() {
		return GainFunctionComparison(this, true);
	}

	virtual int getType() = 0;

protected:
	/** the graph this gain function refers to */
	SignedGraph* graph;
	/** the map of nodes' gain value */
	map<int, GainCalculation> gainMap;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTION_H_ */
