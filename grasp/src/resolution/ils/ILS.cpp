/*
 * ILS.cpp
 *
 *  Created on: 31/03/2014
 *      Author: Mario Levorato
 */

#include "include/ILS.h"
#include "../construction/include/VertexSet.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/Imbalance.h"
#include "graph/include/Perturbation.h"

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
namespace ils {

ILS::ILS() : timeSpentInILS(0.0),
		timeResults(), timeSum(0.0), CBest(),
		numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

ILS::~ILS() {
	// TODO Auto-generated destructor stub
}

Clustering ILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iterMax, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem,	ExecutionInfo& info) {
	BOOST_LOG_TRIVIAL(info) << "Initializing ILS "<< problem.getName() << " procedure for alpha = "
			<< construct->getAlpha() << " and l = " << vnd->getNeighborhoodSize();

	double totalTimeSpentOnConstruction = 0.0;
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	// 1. Initial clustering (construct)
	Clustering CStar = construct->constructClustering(g, problem, info.processRank);
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timeSpentInILS += timeSpentInConstruction;
	totalTimeSpentOnConstruction += timeSpentInConstruction;
	timer.resume();
	start_time = timer.elapsed();

	Clustering previousCc = CStar;
	CBest = CStar;
	Clustering Cc = CStar;
	int iterationValue = 0;
	double timeSpentOnBestSolution = 0.0;
	double initialImbalanceSum = 0.0;
	double timeSpentOnLocalSearch = 0.0;
	stringstream iterationResults;
	stringstream constructivePhaseResults;
	numberOfTestedCombinations = 0;

	// Multi-start ILS
	for (int i = 0; i < iterMax || iterMax < 0 ; i++, previousCc = Cc) {
		BOOST_LOG_TRIVIAL(debug) << "ILS iteration " << i;
		// cout << "Best solution so far: I(P) = " << fixed << setprecision(0) << bestValue.getValue() << endl;

		//    Store initial solution value in corresponding results file
		constructivePhaseResults << (i+1) << "," << Cc.getImbalance().getValue() << ","
				<< Cc.getImbalance().getPositiveValue()
				<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
				<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
		initialImbalanceSum += Cc.getImbalance().getValue();
		// Registers the Cc result
		notifyNewValue(Cc, 0.0, i);

		CStar = Cc;
		Clustering Cl = Cc;
		int perturbationLevel = 1;
		for(int j = 1, total = 0; j <= iterMaxILS; total++) {  // internal ILS loop
			// 2. Execute local search algorithm
			Cl = vnd->localSearch(g, Cl, i, problem, timeSpentInILS, info.processRank);
			numberOfTestedCombinations += vnd->getNumberOfTestedCombinations();
			timeSpentOnLocalSearch += vnd->getTimeSpentOnLocalSearch();
			// 3. Select the best clustring so far
			// if Q(Cl) > Q(Cstar)
			Imbalance newValue = Cl.getImbalance();
			Imbalance bestValue = CStar.getImbalance();
			if(newValue < bestValue) {  // solution improved
				// cout << "A better solution was found." << endl;
				CStar = Cl;
				bestValue = newValue;
				iterationValue = total;
				// Registers the best result at time intervals
				notifyNewValue(CStar, timeSpentInILS, i);
				// restarts the internal ILS loop and the perturbation
				j = 1;
				perturbationLevel = 1;
				timeSpentOnBestSolution = timeSpentInILS;
				if(newValue.getValue() == 0)  break;
			} else {  // did not improve solution
				j++;
				if(j > iterMaxILS) {
					BOOST_LOG_TRIVIAL(debug)<< "Increasing perturbation level...";
					perturbationLevel++;
					j = 1;
					if(perturbationLevel > perturbationLevelMax) {
						break;
					}
				}
			}
			// 4. Generate perturbation over C*
			if((problem.getType() == ClusteringProblem::RCC_PROBLEM) and (CStar.getNumberOfClusters() == 1)) {
				// perturbation does not work if k = 1 and RCC Problem
				break;
			}
			Perturbation perturbation(vnd->getRandomSeed());
			Cl = perturbation.randomMove(g, CStar, problem, perturbationLevel);

			// 5. Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();

			// 6. Write the results into ostream os, using csv format
			// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
			timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
			iterationResults << (total+1) << "," << bestValue.getValue() << "," << bestValue.getPositiveValue()
					<< "," << bestValue.getNegativeValue() << "," << CStar.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInILS << "\n";
			timer.resume();
			start_time = timer.elapsed();
			// if elapsed time is bigger than timeLimit, break
			if(timeSpentInILS >= vnd->getTimeLimit()) {
				BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
				break;
			}
		}
		Imbalance newValue = CStar.getImbalance();
		Imbalance bestValue = CBest.getImbalance();
		if(newValue < bestValue) {
			// cout << "A better solution was found." << endl;
			CBest = CStar;
			bestValue = newValue;
			iterationValue = i;
			timeSpentOnBestSolution = timeSpentInILS;
			if(newValue.getValue() == 0)  break;
		}
		timer.stop();
		end_time = timer.elapsed();
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInILS >= vnd->getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}
		timer.resume();
		start_time = timer.elapsed();
		// guarantees at least one execution of the ILS when the number of iterations is smaller than one
		if(iterMax <= 0) {  break;  }

		// Avoids constructClustering if loop break condition is met
		if((i + 1) < iterMax) {
			// 0. Triggers local processing time calculation
			start_time = timer.elapsed();

			// 1. Construct the next clustering
			Cc = construct->constructClustering(g, problem, info.processRank);

			timer.stop();
			end_time = timer.elapsed();
			timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
			totalTimeSpentOnConstruction += timeSpentInConstruction;
			timeSpentInILS += timeSpentInConstruction;
			timer.resume();
			start_time = timer.elapsed();
		}
	}
	Imbalance bestValue = CBest.getImbalance();
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
			<< "," << bestValue.getPositiveValue()
			<< "," << bestValue.getNegativeValue()
			<< setprecision(0)
			<< "," << CBest.getNumberOfClusters()
			<< "," << (iterationValue+1)
			<< "," << fixed << setprecision(4) << timeSpentOnBestSolution
			<< "," << iterMax
			<< "," << numberOfTestedCombinations << endl;

	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << (initialImbalanceSum / iterMax)
			<< endl;

	BOOST_LOG_TRIVIAL(info) << "ILS procedure done. Obj = " << fixed << setprecision(2) << bestValue.getValue();
	BOOST_LOG_TRIVIAL(info) << "Time spent on construction phase: " << fixed << setprecision(2) << totalTimeSpentOnConstruction << "s, " << (100 * totalTimeSpentOnConstruction / timeSpentInILS) << "%.";
	BOOST_LOG_TRIVIAL(info) << "Time spent on local search: " << fixed << setprecision(2) << timeSpentOnLocalSearch << "s, " << (100 * timeSpentOnLocalSearch / timeSpentInILS) << "%.";
	// CStar.printClustering();
	CStar.printClustering(iterationResults, g->getN());
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);
	// saves the best result to output time file
	measureTimeResults(0.0, iterMax);
	generateOutputFile(problem, timeResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timeIntervals"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);

	return CBest;
}

void ILS::generateOutputFile(ClusteringProblem& problem, stringstream& fileContents, const string& rootFolder,
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
			<< problem.getName() << "-Node" << processNumber << "-l" << l << "a" << std::setprecision(2)
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

void ILS::measureTimeResults(const double& timeSpentOnLocalSearch, const int& ILSIteration) {
	Imbalance imbalance = CBest.getImbalance();
	timeResults << fixed << setprecision(4) << (timeSpentInILS + timeSpentOnLocalSearch) << "," << imbalance.getValue()
			<< "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << "," << CBest.getNumberOfClusters()
			<< "," << (ILSIteration+1) << "," << numberOfTestedCombinations << "\n";
}

unsigned long ILS::getNumberOfTestedCombinations() {
	return numberOfTestedCombinations;
}

} /* namespace ILS */
} /* namespace resolution */
