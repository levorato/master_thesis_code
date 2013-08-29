/*
 * VertexSet.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef VERTEXSET_H_
#define VERTEXSET_H_

#include "GainFunction.h"
#include "../../../graph/include/Clustering.h"

namespace resolution {
namespace grasp {

using namespace std;
using namespace boost;
using namespace clusteringgraph;

class VertexSet {
public:
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
	 * Removes the specified vertex from the list.
	 */
	void removeVertex(int i);
	/**
	 * Picks a random vertex among the first x elements.
	 */
	int chooseRandomVertex(int x);

	/**
	 * Sorts the list according to the gain function.
	 */
	void sort(GainFunction* function);

	GainFunctionVertexSet& getVertexList();
private:
	/**
	 * The random seed used by random vertex choose.
	 */
	unsigned long seed;

	GainFunctionVertexSet vertexSet;

};

} /* namespace grasp */
} /* namespace resolution */
#endif /* VERTEXSET_H_ */
