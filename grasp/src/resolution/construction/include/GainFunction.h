/*
 * GainFunction.h
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#ifndef GAINFUNCTION_H_
#define GAINFUNCTION_H_

#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"

#include <list>
#include <vector>

using namespace clusteringgraph;
using namespace std;

namespace resolution {
namespace construction {

class GainCalculation {
public:
	GainCalculation() : vertex(0), gainValue(0.0), clusterNumber(0) {

	}

	GainCalculation(const GainCalculation &c) : vertex(c.vertex), gainValue(c.gainValue),
			clusterNumber(c.clusterNumber) {

	}

	int vertex;
	// stores the best solution value and cluster number for the vertex
	// to be inserted
	double gainValue;
	int clusterNumber;
};

class GainFunction {
public:
	static const int IMBALANCE = 0, MODULARITY = 1, NEGATIVE_MODULARITY = 2;
	static const int POSITIVE_NEGATIVE_MODULARITY = 3;
	static const int POSITIVE_NEGATIVE_MODULARITY_II = 4;
	static const int POSITIVE_NEGATIVE_MODULARITY_III = 5;

	GainFunction(SignedGraph* g);
	virtual ~GainFunction();

	class GainFunctionComparison {
	private:
		bool ascendingOrder;
	public:
		GainFunctionComparison(bool ascending) :
				ascendingOrder(ascending) {
		}

		bool operator ()(const GainCalculation& a, const GainCalculation& b) const {
			if (ascendingOrder) {
				return a.gainValue < b.gainValue;
			} else {
				return a.gainValue > b.gainValue;
			}
		}
	};

	virtual GainCalculation calculateIndividualGain(ClusteringProblem& p,
			Clustering& c, int i) = 0;

	virtual void calculateGainList(ClusteringProblem& p, Clustering &c,
			list<GainCalculation>& nodeList) = 0;

	virtual GainFunctionComparison getComparator() {
		return GainFunctionComparison(true);
	}

	virtual int getType() = 0;

protected:
	/** the graph this gain function refers to */
	SignedGraph* graph;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GAINFUNCTION_H_ */
