/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: Mario Levorato
 */

#include "include/Grasp.h"
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
namespace grasp {

Grasp::Grasp(GainFunction* f, unsigned long seed) :
		timeSpentInGRASP(0.0), gainFunction(f), randomSeed(seed), timeResults(), timeSum(
				0.0), CBest(), numberOfTestedCombinations(0) {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

Clustering Grasp::executeGRASP(SignedGraph *g, const int& iter,
		const double& alpha, const int& l,
		const bool& firstImprovementOnOneNeig, ClusteringProblem& problem,
		string& executionId, string& fileId, string& outputFolder,
		const long& timeLimit, const int &numberOfSlaves, const int& myRank,
		const int& numberOfSearchSlaves) {
	BOOST_LOG_TRIVIAL(info)<< "Initializing " << problem.getName() << " GRASP "<< problem.getName() <<
	" procedure for alpha = " << alpha << " and l = " << l;
	BOOST_LOG_TRIVIAL(trace) << "Random seed is " << randomSeed << std::endl;

	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();
	Clustering CStar = constructClustering(g, problem, alpha, myRank);
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
	// Chooses between the sequential or parallel VND algorithm
	NeighborhoodSearch* neig;
	NeighborhoodSearchFactory nsFactory(numberOfSlaves, numberOfSearchSlaves);
	if(numberOfSearchSlaves > 0) {
		neig = nsFactory.build(NeighborhoodSearchFactory::PARALLEL);
	} else {
		neig = nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	}
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

		// 2. Execute local search algorithm
		Clustering Cl = localSearch(g, Cc, l, totalIter, firstImprovementOnOneNeig,
				problem, *neig, timeLimit, numberOfSlaves, myRank, numberOfSearchSlaves);
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
			i = 0;
			// TODO validar se essa saida eh valida: nao ha valor de FO menor que zero
			if(newValue.getValue() == 0) break;
		}
		timer.stop();
		end_time = timer.elapsed();
		timeSpentInGRASP += (end_time.wall - start_time.wall) / double(1000000000);
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInGRASP >= timeLimit) {
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
			Cc = constructClustering(g, problem, alpha, myRank);

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
	generateOutputFile(problem, iterationResults, outputFolder, fileId, executionId, myRank, string("iterations"), alpha, l, iter);
	// saves the best result to output time file
	measureTimeResults(0.0, totalIter);
	generateOutputFile(problem, timeResults, outputFolder, fileId, executionId, myRank, string("timeIntervals"), alpha, l, iter);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, outputFolder, fileId, executionId, myRank, string("initialSolutions"), alpha, l, iter);

	return CStar;
}

Clustering Grasp::constructClustering(SignedGraph *g,
		ClusteringProblem& problem, double alpha, int myRank) {
	Clustering Cc; // Cc = empty
	VertexSet lc(randomSeed, g->getN()); // L(Cc) = V(G)
	BOOST_LOG_TRIVIAL(debug)<< "GRASP construct clustering...\n";

	while (lc.size() > 0) { // lc != empty
		GainCalculation gainCalculation;
		if (alpha == 1.0) {
			// alpha = 1.0 (completely random): no need to calculate all gains (saves time)
			int i = lc.chooseRandomVertex(lc.size()).vertex;
			gainCalculation = gainFunction->calculateIndividualGain(problem, Cc,	i);
		} else {
			// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
			// according to the value of the gain function
			gainFunction->calculateGainList(problem, Cc, lc.getVertexList());

			// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
			// (alpha x |lc|) is a rounded number
			// IMPORTANT: if alpha < 0, the constructClustering method will always
			//             choose the first vertex in the gainFunction list, that is,
			//             the one that minimizes the objective (VOTE algorithm).
			if (alpha <= 0.0) {
				// chooses first element (min) without ordering (saves time)
				gainCalculation = lc.chooseFirstElement(gainFunction.get());
			} else {
				lc.sort(gainFunction.get());
				gainCalculation = lc.chooseRandomVertex(
						boost::math::iround(alpha * lc.size()));
			}
			// std::cout << "Random vertex between 0 and " << boost::math::iround(alpha * lc.size()) << " is " << i << std::endl;
		}
		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		// cout << "Selected vertex is " << i << endl;
		// TODO ensure no clusters of size > k are created if RCC Problem
		if (gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			Cc.addCluster(*g, problem, gainCalculation.vertex);
		} else {
			// inserts i into existing cluster
			Cc.addNodeToCluster(*g, problem, gainCalculation.vertex, gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the choosing vertex i automatically removes it from the list
		// Removal already done by the chooseVertex methods above
		// Cc->printClustering();
	}
	BOOST_LOG_TRIVIAL(debug)<< myRank << ": Initial clustering completed.\n";
	Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	// Cc.printClustering();
	return Cc;
}

Clustering Grasp::localSearch(SignedGraph *g, Clustering& Cc, const int &l,
		const int& graspIteration, const bool& firstImprovementOnOneNeig,
		ClusteringProblem& problem, NeighborhoodSearch &neig,
		const long& timeLimit, const int &numberOfSlaves, const int& myRank,
		const int& numberOfSearchSlaves) {
	// k is the current neighborhood distance in the local search
	int r = 1, iteration = 0;
	Clustering CStar = Cc; // C* := Cc
	CBefore = CStar; // previous best solution
	double timeSpentOnLocalSearch = 0.0;
	BOOST_LOG_TRIVIAL(debug)<< "GRASP local search...\n";
	BOOST_LOG_TRIVIAL(debug)<< "Current neighborhood is " << r << endl;

	while (r <= l && (timeSpentInGRASP + timeSpentOnLocalSearch < timeLimit)) {
		// cout << "Local search iteration " << iteration << endl;
		// 0. Triggers local processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// N := Nl(C*)
		// apply a local search in CStar using the r-neighborhood	
		Clustering Cl = neig.searchNeighborhood(r, g, &CStar, problem,
				timeSpentInGRASP + timeSpentOnLocalSearch, timeLimit,
				randomSeed, myRank, firstImprovementOnOneNeig);
		// sums the number of tested combinations on local search
		numberOfTestedCombinations += neig.getNumberOfTestedCombinations();
		if (Cl.getImbalance().getValue() < 0.0000001) {
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
