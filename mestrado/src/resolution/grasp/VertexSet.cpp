/*
 * VertexSet.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/variate_generator.hpp>

#include "include/VertexSet.h"

namespace resolution {
namespace grasp {

VertexSet::VertexSet(unsigned long randomSeed, int n) : seed(randomSeed),
		vertexSetPtr(new GainFunctionVertexSet) {
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
	// boost::random::mt19937 generator;  TODO Adaptar para o modo debug
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	boost::uniform_int<> dist(0,x-1);
	boost::minstd_rand generator(seed);
	generator.seed(boost::random::random_device()());
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uni(generator, dist);
	unsigned int selectedVertexSetIndex = uni();
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

void VertexSet::sort(GainFunction *function) {
	vertexSetPtr->sort(function->getComparator());
	// cout << "Vertex set sorting:" << vertexSetPtr->front() << ", " << function->gain(vertexSetPtr->front()).value << endl;
	// cout << "Vertex set sorting:" << vertexSetPtr->back() << ", " << function->gain(vertexSetPtr->back()).value << endl;
}

list<int>& VertexSet::getVertexList() {
	return *(this->vertexSetPtr.get());
}

} /* namespace grasp */
} /* namespace resolution */
