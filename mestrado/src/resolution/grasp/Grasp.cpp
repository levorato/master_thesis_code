/*
 * Grasp.cpp
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#include "include/Grasp.h"
#include "include/VertexSet.h"
#include "../../graph/include/SequentialNeighborhoodGen.h"
#include "../../problem/include/ClusteringProblem.h"
#include "../../problem/include/CCProblem.h"

#include <algorithm>
#include <iostream>
#include <boost/math/special_functions/round.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace problem;
using namespace boost;
using namespace clusteringgraph;

namespace resolution {
namespace grasp {

Grasp::Grasp() {
	// TODO Auto-generated constructor stub

}

Grasp::~Grasp() {
	// TODO Auto-generated destructor stub
}

ClusteringPtr Grasp::executeGRASP(SignedGraph *g, const int& iter, const float& alpha, const int& l,
		ClusteringProblem* problem, string& timestamp, string& fileId, const int& myRank) {
	std::cout << "Initializing GRASP procedure...\n";
	unsigned int ramdomSeed = 0;
	ClusteringPtr CStar = constructClustering(g, alpha, ramdomSeed);
	ClusteringPtr previousCc = CStar, Cc;
	float bestValue = problem->objectiveFunction(g, CStar.get());
	int iterationValue = 0;
	// TODO alterar o tipo de gerador de vizinhos, para quem sabe, a versao paralelizada
	SequentialNeighborhoodGenerator neig(g->getN());
	stringstream ss;

	// TODO: Parallelize here! Divide iterations by n processors with MPI.
	for (int i = 0; i < iter; i++, previousCc.reset(), previousCc = Cc) {
		cout << "GRASP iteration " << i << endl;
		cout << "Best solution so far: I(P) = " << bestValue << endl;

		// 0. Triggers processing time calculation
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		// 1. Construct the clustering
		Cc = constructClustering(g, alpha, ramdomSeed);

		// 2. Execute local search algorithm
		ClusteringPtr Cl = localSearch(g, *Cc, l, *problem, neig);
		// 3. Select the best clustring so far
		// if Q(Cl) > Q(Cstar)
		float newValue = problem->objectiveFunction(g, Cl.get());

		// 4. Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();

		// TODO consertar metodo equals
		bool same = (previousCc->equals(*Cc.get()));
		// 4. Write the results into ostream os, using csv format
		// Format: iterationNumber,objectiveFunctionValue,time(ms),boolean
		// TODO melhorar formatacao do tempo
		ss << (i+1) << "," << newValue << "," << (end_time.wall - start_time.wall) / 1000000 << "," << same << "\n";

		if(newValue < bestValue) {
			cout << "A better solution was found." << endl;
			CStar.reset();
			CStar = Cl;
			bestValue = newValue;
			iterationValue = i + 1;
			// TODO validar se essa saida eh valida: nao ha valor de FO menor que zero
			if(newValue == 0)  break;
		}
	}
	ss << "Best value: " << bestValue << ", Iteration: " << iterationValue << endl;
	cout << "GRASP procedure done." << endl;
	CStar->printClustering();
	CStar->printClustering(ss);
	generateOutputFile(ss, fileId, timestamp, myRank, alpha, l, iter);
	CStar->setObjectiveFunctionValue(bestValue);

	return CStar;
}

ClusteringPtr Grasp::constructClustering(SignedGraph *g, float alpha, unsigned int ramdomSeed) {
	ClusteringPtr Cc = make_shared<Clustering>(); // Cc = empty
	VertexSet lc(g->getN()); // L(Cc) = V(G)
	std::cout << "GRASP construct clustering...\n";

	while(lc.size() > 0) { // lc != empty
		// cout << "Vertex list size is " << lc.size() << endl;

		// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
		// It is important to calculate the modularity matrix first (used by vertex sorting)
		g->calculateModularityMatrix();
		lc.sort(g, Cc.get());

		// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		int i = lc.chooseRandomVertex(boost::math::iround(alpha * lc.size()));
		// std::cout << "Random vertex is " << i << std::endl;

		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		GainCalculation gainCalculation = Cc->gain(*g, i);
		if(gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			int vertexList[1] = {i};
			Cc->addCluster(vertexList, 1);
		} else {
			// inserts i into existing cluster
			Cc->addNodeToCluster(i, gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the removal of vertex i automatically recalculates the gain function
		lc.removeVertex(i);

		// Cc->printClustering();
	}
	std::cout << "\nInitial clustering completed.\n";
	Cc->printClustering();
	return Cc;
}

ClusteringPtr Grasp::localSearch(SignedGraph *g, Clustering& Cc, const int &l,
		const ClusteringProblem& problem, NeighborhoodListGenerator &neig) {
	// k is the current neighborhood distance in the local search
	int k = 1, iteration = 0;
	ClusteringPtr CStar = make_shared<Clustering>(Cc, g->getN()); // C* := Cc
	std::cout << "GRASP local search...\n";
	cout << "Current neighborhood is " << k << endl;

	while(k <= l) {
		// cout << "Local search iteration " << iteration << endl;
		// N := Nl(C*)
		// apply a local search in CStar using the k-neighborhood

		// TODO Parellelize here!
		ClusteringPtr Cl = neig.generateNeighborhood(k, g, CStar.get(), problem);
		// cout << "Comparing local solution value." << endl;
		if(Cl.get() != NULL && CStar != NULL) {
			if(problem.objectiveFunction(g, Cl.get()) < problem.objectiveFunction(g, CStar.get())) {
				// cout << "New local solution found." << endl;
				// Cl->printClustering();
				CStar.reset();
				CStar = Cl;
				k = 1;
			} else {
				k++;
				cout << "Changed to neighborhood size l = " << k << endl;
			}
		} else {  // no better result found in neighborhood
			k++;
			cout << "Changed to neighborhood size l = " << k << endl;
		}
		iteration++;
	}
	std::cout << "GRASP local search done.\n";
	return CStar;
}

void Grasp::generateOutputFile(stringstream& fileContents, string& fileId,
		string& timestamp, const int &processNumber,
		const float& alpha, const int& l, const int& numberOfIterations) {
	namespace fs = boost::filesystem;
	// Creates the output file (with the results of the execution)
	if (!fs::exists(fs::path("./output"))) {
		fs::create_directories(fs::path("./output"));
	}
	if (!fs::exists(fs::path("./output/" + fileId))) {
		fs::create_directories(
				fs::path("./output/" + fileId));
	}
	if (!fs::exists(fs::path("./output/" + fileId + "/" + timestamp))) {
		fs::create_directories(
				fs::path("./output/" + fileId + "/" + timestamp));
	}
	stringstream filename;
	filename << "./output/" << fileId << "/" << timestamp << "/"
			<< "Node" << processNumber << "-l" << l << "a" << std::setprecision(2) << alpha	<< ".csv";
	fs::path newFile(filename.str());
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		cerr << "Can't open output file!" << endl;
		throw "Cannot open output file.";
	}
	// Writes the parameters to the output file
	// Format: alpha,l,numberOfIterations
	os << std::setprecision(2) << fixed << alpha << "," << l << ","
			<< numberOfIterations << "\n";
	// Writes file contents to the output file
	os << fileContents.str();

	os.close();
}

} /* namespace grasp */
} /* namespace resolution */
