/*
 * VertexSet.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef VERTEXSET_H_
#define VERTEXSET_H_

#include <list>

#include "GainFunction.h"
#include "graph/include/Clustering.h"

namespace resolution {
namespace construction {

using namespace std;
using namespace boost;
using namespace clusteringgraph;

class VertexSet {
public:
	list<GainCalculation> gain;

	/**
	 * Creates an empty set.
	 */
	VertexSet(unsigned long randomSeed);

	/**
	 * Creates a vertex set containing nodes 0..n-1.
	 */
	VertexSet(unsigned long randomSeed, int n);

	virtual ~VertexSet();

	/**
	 * Returns the size of the vertex list.
	 */
	int size();
	/**
	 * Picks a random vertex among the first x elements and removes it.
	 */
	GainCalculation chooseRandomVertex(int x);
	/**
	 * Returns the first element from the list and removes it.
	 */
	GainCalculation chooseFirstElement(GainFunction *function);

	/**
	 * Sorts the list according to the gain function.
	 */
	void sort(GainFunction* function);

	list<GainCalculation>& getVertexList();

	/**
	 * Adds a vertex to the set.
	 */
	void addVertex(int i);

private:
	/**
	 * The random seed used by random vertex choose.
	 */
	unsigned long seed;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* VERTEXSET_H_ */
