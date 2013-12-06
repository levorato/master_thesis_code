/*
 * RCCLocalSearch.cpp
 *
 *  Created on: 28/11/2013
 *      Author: Mario Levorato
 */

#include "include/RCCVariableNeighborhoodSearch.h"
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

RCCVariableNeighborhoodSearch::RCCVariableNeighborhoodSearch(unsigned long seed) : timeResults(),
		randomSeed(seed), timeSum(0.0), timeSpentInSearch(0.0),	numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

RCCVariableNeighborhoodSearch::~RCCVariableNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr RCCVariableNeighborhoodSearch::executeSearch(SignedGraph *g, Clustering& Cc, const int &l,
		const long unsigned int& k, const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
		NeighborhoodSearch &neig, string& executionId, string& fileId, string& outputFolder,
		const long& timeLimit, const int &numberOfSlaves, const int& myRank, const int& numberOfSearchSlaves) {
	// neighborhoodSize is the current neighborhood distance in the local search
	int neighborhoodSize = 1, iteration = 0;
	ClusteringPtr CStar = make_shared<Clustering>(Cc); // C* := Cc
	// Calculates the initial relaxed imbalance of the clustering
	CStar->setImbalance(problem.objectiveFunction(*g, *CStar));

	double timeSpentOnLocalSearch = 0.0;
	BOOST_LOG_TRIVIAL(trace) << "RCC local search...\n";
	BOOST_LOG_TRIVIAL(trace) << "Current neighborhood size is " << neighborhoodSize << endl;

	while(neighborhoodSize <= l && (timeSpentInSearch + timeSpentOnLocalSearch < timeLimit)) {
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// N := Nl(C*)
		// apply a local search in CStar using the neighborhoodSize
		ClusteringPtr Cl = neig.searchNeighborhood(neighborhoodSize, g, CStar.get(), problem,
				timeSpentInSearch + timeSpentOnLocalSearch, timeLimit, randomSeed, myRank, firstImprovementOnOneNeig, k);
		// sums the number of tested combinations on local search
		numberOfTestedCombinations += neig.getNumberOfTestedCombinations();
		// sanity check for obj function value
		if(Cl->getImbalance().getValue() < 0.0) {
			BOOST_LOG_TRIVIAL(error) << myRank << ": Objective function below zero. Error.";
			break;
		}
		Imbalance il = Cl->getImbalance();
		Imbalance ic = CStar->getImbalance();
		if(il < ic) {
			BOOST_LOG_TRIVIAL(trace) << myRank << ": New RCC solution found: " << setprecision(2) << il.getValue() << endl;
			CStar.reset();
			CStar = Cl;
			neighborhoodSize = 1;
		} else {  // no better result found in neighborhood
			neighborhoodSize++;
			BOOST_LOG_TRIVIAL(trace) << "RCC Search: Changed to neighborhood size l = " << neighborhoodSize;
		}
		iteration++;

		// => Finally: Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		timeSpentOnLocalSearch += (end_time.wall - start_time.wall) / double(1000000000);

		// Registers the best result at time intervals
		notifyNewValue(CStar, timeSpentOnLocalSearch, iteration);
	}
	stringstream iterationResults;
	iterationResults << "Best value," << fixed << setprecision(4) << CStar->getImbalance().getValue()
			<< "," << CStar->getImbalance().getPositiveValue()
			<< "," << CStar->getImbalance().getNegativeValue()
			<< setprecision(0)
			<< "," << CStar->getNumberOfClusters()
			<< "," << fixed << setprecision(4) << timeSpentInSearch
			<< "," << (iteration+1)
			<< "," << numberOfTestedCombinations << endl;

	// CStar->printClustering();
	CStar->printClustering(iterationResults);
	generateOutputFile(iterationResults, outputFolder, fileId, executionId, myRank, string("iterations"), k, l);
	// saves the best result to output time file
	measureTimeResults(CStar, 0.0, iteration);
	generateOutputFile(timeResults, outputFolder, fileId, executionId, myRank, string("timeIntervals"), k, l);

	BOOST_LOG_TRIVIAL(trace) << "RCC VNS done.\n";
	return CStar;
}

void RCCVariableNeighborhoodSearch::generateOutputFile(stringstream& fileContents, const string& rootFolder,
		const string& fileId, const string& executionId, const int &processNumber, const string& fileSuffix,
		const unsigned int& k, const int& l) {
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
			<< "RCC-Node" << processNumber << "-l" << l << "-k"
			<< k << "-" << fileSuffix << ".csv";
	fs::path newFile(filename.str());
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		BOOST_LOG_TRIVIAL(fatal) << "Can't open output file!" << endl;
		throw "Cannot open output file.";
	}
	// Writes the parameters to the output file
	// Format: k,l,numberOfIterations,numberOfCombinations
	os << std::setprecision(2) << k << "," << l << "\n";
	// Writes file contents to the output file
	os << fileContents.str();

	os.close();
}

void RCCVariableNeighborhoodSearch::measureTimeResults(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	Imbalance imbalance = CStar->getImbalance();
	timeResults << fixed << setprecision(4) << (timeSpentInSearch + timeSpentOnLocalSearch) << "," << imbalance.getValue()
			<< "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << "," << CStar->getNumberOfClusters()
			<< "," << (iteration+1) << "\n";
}

void RCCVariableNeighborhoodSearch::notifyNewValue(ClusteringPtr CStar, const double& timeSpentOnLocalSearch, const int& iteration) {
	measureTimeResults(CStar, timeSpentOnLocalSearch, iteration);
}

unsigned long RCCVariableNeighborhoodSearch::getNumberOfTestedCombinations() {
	return numberOfTestedCombinations;
}


} /* namespace vnd */
} /* namespace resolution */
