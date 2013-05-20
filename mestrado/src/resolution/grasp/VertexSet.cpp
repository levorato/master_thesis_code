/*
 * VertexSet.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "include/VertexSet.h"

namespace resolution {
namespace grasp {

VertexSet::VertexSet(int n) : vertexSetPtr(new GainFunctionVertexSet) {
	for(int i = 0; i < n; i++) {
		vertexSetPtr->push_back(i);
	}
}

VertexSet::~VertexSet() {
	// TODO Auto-generated destructor stub
}

int VertexSet::size() {
	return vertexSetPtr->size();
}

void VertexSet::removeVertex(int i) {
	vertexSetPtr->remove(i);
}

// TODO aceitar parametro seed para a geracao do numero aleatorio
int VertexSet::chooseRandomVertex(int x) {
	// Generates a random number between 1 and x
	boost::random::mt19937 gen;
	// distribution that maps to 1..x
	boost::random::uniform_int_distribution<> dist(1,x);
	int selectedVertexSetIndex = dist(gen);
	int selectedVertex = 0;

	// Returns the Vertex
	list<int, allocator<int> >::const_iterator pos;
	list<int, allocator<int> > vertexSet = *vertexSetPtr.get();
	unsigned int i = 0;
	for(i = 0, pos = vertexSet.begin(); i < vertexSet.size(); ++pos, ++i) {
		if(i == selectedVertexSetIndex) {
			selectedVertex = *pos;
			break;
		}
	}
	return selectedVertex;
}

void VertexSet::sort(SignedGraph* g, Clustering* c) {
	vertexSetPtr->sort(GainFunctionComparison(g, c));
}

} /* namespace grasp */
} /* namespace resolution */
