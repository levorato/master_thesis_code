/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: Mario Levorato
 */

#include "include/Grasp.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include "graph/include/Imbalance.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
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
namespace grasp {

Grasp::Grasp() :
		timeSpentInGRASP(0.0), timeResults(), timeSum(
				0.0), CBest(), numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

Clustering Grasp::executeGRASP(ConstructClustering &construct, VariableNeighborhoodDescent &vnd,
		SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info) {
	BOOST_LOG_TRIVIAL(info)<< "Initializing " << " GRASP "<< problem.getName() <<
	" procedure for alpha = " << construct.getAlpha() << " and l = " << vnd.getNeighborhoodSize();

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	Clustering CStar = construct.constructClustering(g, problem, info.processRank);
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timeSpentInGRASP += timeSpentInConstruction;
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

	while(i <= iter || iter < 0) {
		BOOST_LOG_TRIVIAL(debug) << "GRASP iteration " << i;
		// cout << "Best solution so far: I(P) = " << fixed << setprecision(0) << bestValue.getValue() << endl;

		//    Store initial solution value in corresponding results file
		constructivePhaseResults << (totalIter+1) << "," << Cc.getImbalance().getValue() << "," << Cc.getImbalance().getPositiveValue()
		<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
		<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
		initialImbalanceSum += Cc.getImbalance().getValue();
		// Registers the Cc result
		notifyNewValue(Cc, 0.0, totalIter);

		// 2. Execute local search algorithm: RVND
		Clustering Cl = vnd.localSearch(g, Cc, totalIter, problem, timeSpentInGRASP, info.processRank);
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		Imbalance newValue = Cl.getImbalance();

		// 4. Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();

		// 5. Write the results into ostream os, using csv format
		// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
		// TODO melhorar formatacao do tempo
		timeSpentInGRASP += (end_time.wall - start_time.wall) / double(1000000000);
		iterationResults << (totalIter+1) << "," << newValue.getValue() << "," << newValue.getPositiveValue()
		<< "," << newValue.getNegativeValue() << "," << CStar.getNumberOfClusters()
		<< "," << fixed << setprecision(4) << timeSpentInGRASP << "\n";
		timer.resume();
		start_time = timer.elapsed();

		if(newValue < bestValue) {
			// cout << "A better solution was found." << endl;
			CStar = Cl;
			bestValue = newValue;
			iterationValue = totalIter;
			timeSpentOnBestSolution = timeSpentInGRASP;
			// Registers the best result at time intervals
			notifyNewValue(CStar, timeSpentInGRASP, totalIter);
			i = 0;
			// TODO validar se essa saida eh valida: nao ha valor de FO menor que zero
			if(newValue.getValue() == 0) break;
		}
		timer.stop();
		end_time = timer.elapsed();
		timeSpentInGRASP += (end_time.wall - start_time.wall) / double(1000000000);
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInGRASP >= vnd.getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}

		// Increment to next loop
		i++, totalIter++, previousCc = Cc;

		// guarantees at least one execution of the GRASP when the number of iterations is smaller than one
		// used in vote/boem tests
		if(iter <= 0) {break;}

		// Avoids constructClustering if loop break condition is met
		if(i <= iter) {
			// 0. Triggers local processing time calculation
			timer.resume();
			start_time = timer.elapsed();

			// 1. Construct the next clustering
			Cc = construct.constructClustering(g, problem, info.processRank);

			timer.stop();
			end_time = timer.elapsed();
			timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
			timeSpentInGRASP += timeSpentInConstruction;
			timer.resume();
			start_time = timer.elapsed();
		}
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

	BOOST_LOG_TRIVIAL(info) << "GRASP procedure done. Obj = " << fixed << setprecision(2) << bestValue.getValue();
	// CStar.printClustering();
	CStar.printClustering(iterationResults);
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);
	// saves the best result to output time file
	measureTimeResults(0.0, totalIter);
	generateOutputFile(problem, timeResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timeIntervals"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions"), construct.getAlpha(), vnd.getNeighborhoodSize(), iter);

	return CStar;
}

void Grasp::generateOutputFile(ClusteringProblem& problem,
		stringstream& fileContents, const string& rootFolder,
		const string& fileId, const string& executionId,
		const int &processNumber, const string& fileSuffix, const double& alpha,
		const int& l, const int& numberOfIterations) {
	namespace fs = boost::filesystem;
	// Creates the output file (with the results of the execution)
	if (!fs::exists(fs::path(rootFolder))) {
		fs::create_directories(rootFolder);
	}
	string outputFolder = rootFolder + "/";
	if (!fs::exists(fs::path(outputFolder + fileId))) {
		fs::create_directories(fs::path(outputFolder + fileId));
	}
	if (!fs::exists(fs::path(outputFolder + fileId + "/" + executionId))) {
		fs::create_directories(
				fs::path(outputFolder + fileId + "/" + executionId));
	}
	stringstream filename;
	filename << outputFolder << fileId << "/" << executionId << "/"
			<< problem.getName() << "-Node" << processNumber << "-l" << l << "a"
			<< std::setprecision(2) << alpha << "-" << fileSuffix << ".csv";
	fs::path newFile(filename.str());
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		BOOST_LOG_TRIVIAL(fatal)<< "Can't open output file!" << endl;
		throw "Cannot open output file.";
	}
	// Writes the parameters to the output file
	// Format: alpha,l,numberOfIterations,numberOfCombinations
	os << std::setprecision(2) << alpha << "," << l << "," << numberOfIterations
			<< "\n";
	// Writes file contents to the output file
	os << fileContents.str();

	os.close();
}

void Grasp::measureTimeResults(const double& timeSpentOnLocalSearch,
		const int& graspIteration) {
	Imbalance imbalance = CBest.getImbalance();
	timeResults << fixed << setprecision(4)
			<< (timeSpentInGRASP + timeSpentOnLocalSearch) << ","
			<< imbalance.getValue() << "," << imbalance.getPositiveValue()
			<< "," << imbalance.getNegativeValue() << ","
			<< CBest.getNumberOfClusters() << "," << (graspIteration + 1) << ","
			<< numberOfTestedCombinations << "\n";
}

unsigned long Grasp::getNumberOfTestedCombinations() {
	return numberOfTestedCombinations;
}

} /* namespace grasp */
} /* namespace resolution */
