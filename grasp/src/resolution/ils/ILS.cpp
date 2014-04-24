/*
 * ILS.cpp
 *
 *  Created on: 31/03/2014
 *      Author: Mario Levorato
 */

#include "include/ILS.h"
#include "../construction/include/VertexSet.h"
#include "../../graph/include/SequentialNeighborhoodSearch.h"
#include "../../graph/include/ParallelNeighborhoodSearch.h"
#include "../../problem/include/ClusteringProblem.h"
#include "../../problem/include/CCProblem.h"
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
namespace ils {

ILS::ILS() : timeSpentInILS(0.0),
		timeResults(), timeSum(0.0), CBest(),
		numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

ILS::~ILS() {
	// TODO Auto-generated destructor stub
}

Clustering ILS::executeILS(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
		SignedGraph *g, const int& iter, ClusteringProblem& problem, string& executionId,
		string& fileId, string& outputFolder, const int& myRank) {
	BOOST_LOG_TRIVIAL(info) << "Initializing ILS "<< problem.getName() << " procedure for alpha = " << construct.getAlpha() << " and l = " << vnd.getNeighborhoodSize();

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	Clustering CStar = construct.constructClustering(g, problem, myRank);
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timeSpentInILS += timeSpentInConstruction;
	timer.resume();
	start_time = timer.elapsed();

	Clustering previousCc = CStar;
	CBest = CStar;
	Clustering Cc = CStar;
	Imbalance bestValue = CStar.getImbalance();
	int iterationValue = 0;
	double timeSpentOnBestSolution = 0.0;
	double initialImbalanceSum = 0.0;
	stringstream iterationResults;
	stringstream constructivePhaseResults;
	int i = 0, totalIter = 0;
	numberOfTestedCombinations = 0;

	for (i = 0, totalIter = 0; i <= iter || iter < 0 ; i++, totalIter++, previousCc = Cc) {
		BOOST_LOG_TRIVIAL(debug) << "ILS iteration " << i;
		// cout << "Best solution so far: I(P) = " << fixed << setprecision(0) << bestValue.getValue() << endl;

		//    Store initial solution value in corresponding results file
		constructivePhaseResults << (totalIter+1) << "," << Cc.getImbalance().getValue() << "," << Cc.getImbalance().getPositiveValue()
						<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
						<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
		initialImbalanceSum += Cc.getImbalance().getValue();
		// Registers the Cc result
		notifyNewValue(Cc, 0.0, totalIter);

		// 2. Execute local search algorithm
		Clustering Cl = vnd.localSearch(g, Cc, totalIter, problem, timeSpentInILS, myRank);
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		Imbalance newValue = Cl.getImbalance();

		// 4. Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();

		// 5. Write the results into ostream os, using csv format
		// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
		// TODO melhorar formatacao do tempo
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		iterationResults << (totalIter+1) << "," << newValue.getValue() << "," << newValue.getPositiveValue()
				<< "," << newValue.getNegativeValue() << "," << CStar.getNumberOfClusters()
				<< "," << fixed << setprecision(4) << timeSpentInILS << "\n";
		timer.resume();
		start_time = timer.elapsed();

		if(newValue < bestValue) {
			// cout << "A better solution was found." << endl;
			CStar = Cl;
			bestValue = newValue;
			iterationValue = totalIter;
			timeSpentOnBestSolution = timeSpentInILS;
			i = 0;
			// TODO validar se essa saida eh valida: nao ha valor de FO menor que zero
			if(newValue.getValue() == 0)  break;
		}
		timer.stop();
		end_time = timer.elapsed();
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInILS >= vnd.getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}

		// 0. Triggers local processing time calculation
		timer.resume();
		start_time = timer.elapsed();

		// 1. Construct the next clustering
		Cc = construct.constructClustering(g, problem, myRank);

		timer.stop();
		end_time = timer.elapsed();
		timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
		timer.resume();

		// guarantees at least one execution of the ILS when the number of iterations is smaller than one
		if(iter <= 0) {  break;  }
	}
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
			<< "," << bestValue.getPositiveValue()
			<< "," << bestValue.getNegativeValue()
			<< setprecision(0)
			<< "," << CStar.getNumberOfClusters()
			<< "," << (iterationValue+1)
			<< "," << fixed << setprecision(4) << timeSpentOnBestSolution
			<< "," << (totalIter+1) 
			<< "," << numberOfTestedCombinations << endl;

	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << (initialImbalanceSum / (totalIter+1))
			<< endl;

	BOOST_LOG_TRIVIAL(info) << "ILS procedure done. Obj = " << fixed << setprecision(2) << bestValue.getValue();
	// CStar.printClustering();
	CStar.printClustering(iterationResults);
	generateOutputFile(problem, iterationResults, outputFolder, fileId, executionId, myRank, string("iterations"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);
	// saves the best result to output time file
	measureTimeResults(0.0, totalIter);
	generateOutputFile(problem, timeResults, outputFolder, fileId, executionId, myRank, string("timeIntervals"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, outputFolder, fileId, executionId, myRank, string("initialSolutions"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);

	return CStar;
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
