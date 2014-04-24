/*
 * RCCLocalSearch.cpp
 *
 *  Created on: 28/11/2013
 *      Author: Mario Levorato
 */

#include "include/VariableNeighborhoodDescent.h"
#include "../../graph/include/SequentialNeighborhoodSearch.h"
#include "../../graph/include/ParallelNeighborhoodSearch.h"
#include "../../problem/include/ClusteringProblem.h"
#include "../../graph/include/NeighborhoodSearchFactory.h"
#include "../../graph/include/Imbalance.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/round.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/log/trivial.hpp>


using namespace problem;
using namespace boost;
using namespace clusteringgraph;

namespace resolution {
namespace vnd {

VariableNeighborhoodDescent::VariableNeighborhoodDescent(unsigned long seed) : timeResults(),
		randomSeed(seed), timeSum(0.0), timeSpentInSearch(0.0),	numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

VariableNeighborhoodDescent::~VariableNeighborhoodDescent() {
	// TODO Auto-generated destructor stub
}

Clustering VariableNeighborhoodDescent::localSearch(SignedGraph *g, Clustering& Cc, const int &l,
		const int& graspIteration, const bool& firstImprovementOnOneNeig,
		ClusteringProblem& problem, NeighborhoodSearch &neig,
		const long& timeLimit, const long& timeSpentSoFar, const int &numberOfSlaves, const int& myRank,
		const int& numberOfSearchSlaves) {
	// k is the current neighborhood distance in the local search
	int r = 1, iteration = 0;
	Clustering CStar = Cc; // C* := Cc
	double timeSpentOnLocalSearch = 0.0;
	BOOST_LOG_TRIVIAL(debug)<< "GRASP local search...\n";
	BOOST_LOG_TRIVIAL(debug)<< "Current neighborhood is " << r << endl;

	while (r <= l && (timeSpentSoFar + timeSpentOnLocalSearch < timeLimit)) {
		// cout << "Local search iteration " << iteration << endl;
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// N := Nl(C*)
		// apply a local search in CStar using the k-neighborhood
		Clustering Cl = neig.searchNeighborhood(r, g, &CStar, problem,
				timeSpentSoFar + timeSpentOnLocalSearch, timeLimit,
				randomSeed, myRank, firstImprovementOnOneNeig);
		// sums the number of tested combinations on local search
		numberOfTestedCombinations += neig.getNumberOfTestedCombinations();
		if (Cl.getImbalance().getValue() < 0.0) {
			BOOST_LOG_TRIVIAL(error)<< myRank << ": Objective function below zero. Obj = " << Cl.getImbalance().getValue();
			break;
		}
		// cout << "Comparing local solution value." << endl;
		Imbalance il = Cl.getImbalance();
		Imbalance ic = CStar.getImbalance();
		if (il < ic) {
			// BOOST_LOG_TRIVIAL(trace) << myRank << ": New local solution found: " << setprecision(2) << il.getValue() << endl;
			// Cl->printClustering();
			CStar = Cl;
			r = 1;
		} else {  // no better result found in neighborhood
			r++;
			// BOOST_LOG_TRIVIAL(debug) << "Changed to neighborhood size l = " << k;
		}
		iteration++;

		// => Finally: Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timeSpentOnLocalSearch += (end_time.wall - start_time.wall)
				/ double(1000000000);

		// Registers the best result at time intervals
		notifyNewValue(CStar, timeSpentOnLocalSearch, graspIteration);
	}
	BOOST_LOG_TRIVIAL(debug)<< "GRASP local search done. Time spent: " << timeSpentOnLocalSearch << " s";
	return CStar;
}

void VariableNeighborhoodDescent::measureTimeResults(Clustering &CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	Imbalance imbalance = CStar.getImbalance();
	timeResults << fixed << setprecision(4) << (timeSpentInSearch + timeSpentOnLocalSearch) << "," << imbalance.getValue()
			<< "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << "," << CStar.getNumberOfClusters()
			<< "," << (iteration+1) << "\n";
}

void VariableNeighborhoodDescent::notifyNewValue(Clustering& CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	measureTimeResults(CStar, timeSpentOnLocalSearch, iteration);
}

unsigned long VariableNeighborhoodDescent::getNumberOfTestedCombinations() {
	return numberOfTestedCombinations;
}


} /* namespace vnd */
} /* namespace resolution */
