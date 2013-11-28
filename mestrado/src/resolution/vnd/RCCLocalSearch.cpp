/*
 * RCCLocalSearch.cpp
 *
 *  Created on: 28/11/2013
 *      Author: czt0
 */

#include "include/RCCLocalSearch.h"
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

RCCLocalSearch::RCCLocalSearch() : timeResults(), randomSeed(seed), timeSum(0.0), timeSpentInSearch(0.0)
		numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

RCCLocalSearch::~RCCLocalSearch() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr RCCLocalSearch::localSearch(SignedGraph *g, Clustering& Cc, const int &l,
		const int& graspIteration,
		const bool& firstImprovementOnOneNeig, const ClusteringProblem& problem,
		NeighborhoodSearch &neig, const long& timeLimit, const int &numberOfSlaves,
		const int& myRank, const int& numberOfSearchSlaves) {
	// k is the current neighborhood distance in the local search
	int k = 1, iteration = 0;
	ClusteringPtr CStar = make_shared<Clustering>(Cc); // C* := Cc

	double timeSpentOnLocalSearch = 0.0;
	BOOST_LOG_TRIVIAL(trace) << "GRASP local search...\n";
	BOOST_LOG_TRIVIAL(trace) << "Current neighborhood is " << k << endl;

	while(k <= l && (timeSpentInSearch + timeSpentOnLocalSearch < timeLimit)) {
		// cout << "Local search iteration " << iteration << endl;
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// N := Nl(C*)
		// apply a local search in CStar using the k-neighborhood
		ClusteringPtr Cl = neig.searchNeighborhood(k, g, CStar.get(), problem,
				timeSpentInSearch + timeSpentOnLocalSearch, timeLimit, randomSeed, myRank, firstImprovementOnOneNeig);
		// sums the number of tested combinations on local search
		numberOfTestedCombinations += neig.getNumberOfTestedCombinations();
		if(Cl->getImbalance().getValue() < 0.0) {
			BOOST_LOG_TRIVIAL(error) << myRank << ": Objective function below zero. Error.";
			break;
		}
		// cout << "Comparing local solution value." << endl;
		Imbalance il = Cl->getImbalance();
		Imbalance ic = CStar->getImbalance();
		if(il < ic) {
			BOOST_LOG_TRIVIAL(trace) << myRank << ": New local solution found: " << setprecision(2) << il.getValue() << endl;
			// Cl->printClustering();
			CStar.reset();
			CStar = Cl;
			k = 1;
		} else {  // no better result found in neighborhood
			k++;
			BOOST_LOG_TRIVIAL(trace) << "Changed to neighborhood size l = " << k;
		}
		iteration++;

		// => Finally: Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timeSpentOnLocalSearch += (end_time.wall - start_time.wall) / double(1000000000);

		// Registers the best result at time intervals
		notifyNewValue(CStar, timeSpentOnLocalSearch, graspIteration);
	}
	BOOST_LOG_TRIVIAL(trace) << "GRASP local search done.\n";
	return CStar;
}

void RCCLocalSearch::generateOutputFile(stringstream& fileContents, const string& rootFolder,
		const string& fileId, const string& executionId, const int &processNumber, const string& fileSuffix,
		const double& alpha, const int& l, const int& numberOfIterations) {
	namespace fs = boost::filesystem;
	// Creates the output file (with the results of the execution)
	if (!fs::exists(fs::path(rootFolder))) {
		fs::create_directories(rootFolder);
	}
	string outputFolder = rootFolder + "/";
	if (!fs::exists(fs::path(outputFolder + fileId))) {
		fs::create_directories(
				fs::path(outputFolder + fileId));
	}
	if (!fs::exists(fs::path(outputFolder + fileId + "/" + executionId))) {
		fs::create_directories(
				fs::path(outputFolder + fileId + "/" + executionId));
	}
	stringstream filename;
	filename << outputFolder << fileId << "/" << executionId << "/"
			<< "Node" << processNumber << "-l" << l << "a" << std::setprecision(2)
			<< alpha << "-" << fileSuffix << ".csv";
	fs::path newFile(filename.str());
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		BOOST_LOG_TRIVIAL(fatal) << "Can't open output file!" << endl;
		throw "Cannot open output file.";
	}
	// Writes the parameters to the output file
	// Format: alpha,l,numberOfIterations,numberOfCombinations
	os << std::setprecision(2) << alpha << "," << l << ","
			<< numberOfIterations << "\n";
	// Writes file contents to the output file
	os << fileContents.str();

	os.close();
}

void RCCLocalSearch::measureTimeResults(const double& timeSpentOnLocalSearch, const int& graspIteration) {
	Imbalance imbalance = CBest->getImbalance();
	timeResults << fixed << setprecision(4) << (timeSpentInSearch + timeSpentOnLocalSearch) << "," << imbalance.getValue()
			<< "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << "," << CBest->getNumberOfClusters()
			<< "," << (graspIteration+1) << "\n";
}

void RCCLocalSearch::notifyNewValue(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& graspIteration) {
	Imbalance i1 = CStar->getImbalance();
	Imbalance i2 = CBest->getImbalance();
	if(i1 < i2) {
		CBest = CStar;
		measureTimeResults(timeSpentOnLocalSearch, graspIteration);
	}
}

long RCCLocalSearch::getNumberOfTestedCombinations() {
	return numberOfTestedCombinations;
}


} /* namespace vnd */
} /* namespace resolution */
