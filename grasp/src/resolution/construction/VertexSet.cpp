/*
 * VertexSet.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include <algorithm>

#include "include/VertexSet.h"
#include "util/include/RandomUtil.h"

namespace resolution {
namespace construction {

using namespace util;

VertexSet::VertexSet(unsigned long randomSeed) : seed(randomSeed) {
	// Empty set
}

VertexSet::VertexSet(unsigned long randomSeed, int n) : seed(randomSeed),
		gain() {
	for(int i = 0; i < n; i++) {
		GainCalculation v;
		v.vertex = i;
		v.clusterNumber = 0;
		v.gainValue = 0.0;
		gain.push_back(v);
	}
}

VertexSet::~VertexSet() {
	// TODO Auto-generated destructor stub
}

int VertexSet::size() {
	return gain.size();
}

// TODO aceitar parametro seed para a geracao do numero aleatorio
GainCalculation VertexSet::chooseRandomVertex(int x) {

	// Generates a random number between 1 and x
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	unsigned int selectedVertexSetIndex = randomUtil.next(0, x - 1);
	GainCalculation selectedVertex;

	// Finds the Vertex
	list<GainCalculation, allocator<GainCalculation> >::iterator pos = gain.begin();
	std::advance(pos, selectedVertexSetIndex);
	selectedVertex = *pos;
	gain.erase(pos);
	return selectedVertex;
}

GainCalculation VertexSet::chooseFirstElement(GainFunction *function) {
	list<GainCalculation, allocator<GainCalculation> >::iterator pos =
			std::min_element(gain.begin(), gain.end(), function->getComparator());
	GainCalculation x = *pos;
	gain.erase(pos);
	return x;
}

void VertexSet::sort(GainFunction *function) {
	gain.sort(function->getComparator());
	// cout << "Vertex set sorting:" << gainPtr->front() << ", " << function->gain(gainPtr->front()).value << endl;
	// cout << "Vertex set sorting:" << gainPtr->back() << ", " << function->gain(gainPtr->back()).value << endl;
}

list<GainCalculation>& VertexSet::getVertexList() {
	return gain;
}

void VertexSet::addVertex(int i) {
	GainCalculation v;
	v.vertex = i;
	v.clusterNumber = 0;
	v.gainValue = 0.0;
	gain.push_back(v);
}

} /* namespace grasp */
} /* namespace resolution */
