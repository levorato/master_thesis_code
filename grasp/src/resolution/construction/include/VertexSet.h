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
#include "../../../graph/include/Clustering.h"

namespace resolution {
namespace construction {

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
	 * Picks a random vertex among the first x elements and removes it.
	 */
	int chooseRandomVertex(int x);
	/**
	 * Returns the first element from the list and removes it.
	 */
	int chooseFirstElement();

	/**
	 * Sorts the list according to the gain function.
	 */
	void sort(GainFunction& function);

	list<int>& getVertexList();
private:
	/**
	 * The random seed used by random vertex choose.
	 */
	unsigned long seed;

	list<int> vertexSet;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* VERTEXSET_H_ */
