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
#include <algorithm>

#include "include/VertexSet.h"

namespace resolution {
namespace construction {

VertexSet::VertexSet(unsigned long randomSeed, int n) : seed(randomSeed),
		vertexSet() {
	for(int i = 0; i < n; i++) {
		vertexSet.push_back(i);
	}
}

VertexSet::~VertexSet() {
	// TODO Auto-generated destructor stub
}

int VertexSet::size() {
	return vertexSet.size();
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

	// Finds the Vertex
	list<int, allocator<int> >::iterator pos = vertexSet.begin();
	std::advance(pos, selectedVertexSetIndex);
	selectedVertex = *pos;
	vertexSet.erase(pos);
	return selectedVertex;
}

int VertexSet::chooseFirstElement(GainFunction &function) {
	list<int, allocator<int> >::iterator pos = std::min_element(vertexSet.begin(), vertexSet.end(), function.getComparator());
	int x = *pos;
	vertexSet.erase(pos);
	return x;
}

void VertexSet::sort(GainFunction &function) {
	vertexSet.sort(function.getComparator());
	// cout << "Vertex set sorting:" << vertexSetPtr->front() << ", " << function->gain(vertexSetPtr->front()).value << endl;
	// cout << "Vertex set sorting:" << vertexSetPtr->back() << ", " << function->gain(vertexSetPtr->back()).value << endl;
}

list<int>& VertexSet::getVertexList() {
	return vertexSet;
}

} /* namespace grasp */
} /* namespace resolution */
