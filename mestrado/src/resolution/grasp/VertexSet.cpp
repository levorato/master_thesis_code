/*
 * VertexSet.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/nondet_random.hpp>

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
	/*
		* Caveat: std::time(0) is not a very good truly-random seed.  When
		* called in rapid succession, it could return the same values, and
		* thus the same random number sequences could ensue.
		* Instead, we are using boost::random_device
		* http://stackoverflow.com/questions/4329284/c-boost-random-numeric-generation-problem
	    */
	boost::minstd_rand generator(1234u);
	// consider using seed + (long long)getpid() << 32 with more than one process (MPI)
	generator.seed(boost::random::random_device()());

	// Generates a random number between 1 and x
	// boost::random::mt19937 generator;  TODO Adaptar para o modo deug
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	boost::random::uniform_int_distribution<> dist(0,x-1);
	unsigned int selectedVertexSetIndex = dist(generator);
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
