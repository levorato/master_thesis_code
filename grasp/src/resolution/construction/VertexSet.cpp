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
