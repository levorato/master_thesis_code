/*
 * ImbalanceGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ImbalanceGainFunction.h"

#include <vector>
#include <algorithm>
#include <boost/math/special_functions/round.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/variate_generator.hpp>

namespace resolution {
namespace grasp {

using namespace std;

ImbalanceGainFunction::ImbalanceGainFunction(SignedGraph* g, const unsigned long& s) :
		GainFunction::GainFunction(g, s) {
	// TODO Auto-generated constructor stub

}

ImbalanceGainFunction::~ImbalanceGainFunction() {
	// TODO Auto-generated destructor stub
}

void ImbalanceGainFunction::calculateGainList(Clustering &c, GainFunctionVertexSet& nodeList) {
	gainMap.clear();
	// cout << "Calculating gain list..." << endl;
	list<int, allocator<int> >::const_iterator pos;
	unsigned int i = 0;
	for(i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
	    int a = *pos;
	    // gainList contains the list of all possible gains for vertex a
	    vector<GainCalculation> gainList;
		// cout << "Vertex " << a << endl;

		// For each existing cluster k...
		int nc = c.getNumberOfClusters();
		for(int k = 0; k < nc; k++) {
			// cout << "Cluster " << k << endl;
			GainCalculation gainCalculation;
			Imbalance delta = c.calculateDeltaObjectiveFunction(*graph, c.getCluster(k), a);
			gainCalculation.value = delta.getValue();
			gainCalculation.clusterNumber = k;
			gainList.push_back(gainCalculation);
		}
		// For a new cluster (k+1)
		// cout << "New cluster" << endl;
		BoolArray newCluster(graph->getN());
		newCluster[a] = true;
		GainCalculation gainCalculation;
		Imbalance delta = c.calculateDeltaObjectiveFunction(*graph, newCluster, a);
		gainCalculation.value = delta.getValue();
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
		gainList.push_back(gainCalculation);

		// Sorts the gainList according to the gain value (ascending order)
		std::sort(gainList.begin(), gainList.end(), GainCalculationComparison(true));

		// Chooses randomly one in the first (alpha x |gainList.size()|) destination clusters in gainList
		int x = chooseRandomNumber(boost::math::iround(c.getAlpha() * gainList.size()));
		gainMap[a] = gainList[x];
	}
}

/**
 * TODO For a given vertex a, calculates the minimum value of imbalance (I(P))
 * of inserting 'a' into a new or an existing clustering k. Returns the minimum imbalance
 * and the cluster corresponding to it.
 */
GainCalculation& ImbalanceGainFunction::gain(const int &a) {
	return gainMap[a];
}

bool ImbalanceGainFunction::operator () ( const int& a, const int& b ) {
	return this->gain(a).value < this->gain(b).value;
}

int ImbalanceGainFunction::getType() {
	return GainFunction::IMBALANCE;
}

GainFunction::GainFunctionComparison ImbalanceGainFunction::getComparator() {
	return GainFunctionComparison(this, true);
}

unsigned int ImbalanceGainFunction::chooseRandomNumber(int x) {

	// Generates a random number between 1 and x
	// boost::random::mt19937 generator;  TODO Adaptar para o modo debug
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	boost::uniform_int<> dist(0,x-1);
	boost::minstd_rand generator(randomSeed);
	generator.seed(boost::random::random_device()());
	boost::variate_generator<minstd_rand&, boost::uniform_int<> > uni(generator, dist);
	unsigned int selectedIndex = uni();

	return selectedIndex;
}

} /* namespace grasp */
} /* namespace resolution */
