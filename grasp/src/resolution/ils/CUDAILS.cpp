/*
 * CUDAILS.cpp
 *
 *  Created on: Mar 2, 2015
 *      Author: mlevorato
 */

#include "include/CUDAILS.h"
#include "../construction/include/VertexSet.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/Imbalance.h"
#include "graph/include/Perturbation.h"
#include "../vnd/include/CUDASearch.h"

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
#include <thrust/host_vector.h>


using namespace problem;
using namespace boost;
using namespace clusteringgraph;
using namespace thrust;
using namespace util;

namespace resolution {
namespace ils {

CUDAILS::CUDAILS() {
	// TODO Auto-generated constructor stub

}

CUDAILS::~CUDAILS() {
	// TODO Auto-generated destructor stub
}

Clustering CUDAILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iterMax, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem,	ExecutionInfo& info) {
	BOOST_LOG_TRIVIAL(info) << "Initializing CUDA ILS "<< problem.getName() << " procedure for alpha = "
			<< construct->getAlpha() << " and l = " << vnd->getNeighborhoodSize();

	stringstream iterationResults;
	stringstream constructivePhaseResults;

	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// CUDA variables
	unsigned long n = g->getN();
	unsigned long m = g->getM();
	// CUDA graph structure (adapted adjacency list)
	thrust::host_vector<float> h_weights(2 * m);  // in/out edge weights
	thrust::host_vector<int> h_dest(2 * m);  // edge destination (vertex j)
	thrust::host_vector<int> h_numedges(n);  // number of edges of each vertex i
	thrust::host_vector<int> h_offset(n);  // initial edge number for vertex i
	// For each vertex, creates a list of in and out edges
	int i = 0, offset = 0;
	for(long edge = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		h_offset[i] = offset;
		for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {  // out edges of i
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, g->graph); f2 != l2; ++f2) {  // in edges of i
			double weight = ((Edge*)f2->get_property())->weight;
			int j = source(*f2, g->graph);
			h_dest[edge] = j;
			h_weights[edge] = weight;
			count++; edge++;
		}
		h_numedges[i] = count;
		offset += count;
	}
	// TODO transform into class constant
	// number of threads per block
	unsigned short threadsCount = 1024;  // limited by shared memory size
	// Pass raw array and its size to kernel
	int totalIter = 0;
	Clustering CStar;
	CStar.setImbalance(Imbalance(-1.0, -1.0));
	double totalTimeSpentOnConstruction = 0, timeSpentOnLocalSearch = 0;
	runILSKernel(problem, *construct, g, info.processRank, vnd->getTimeLimit(),
			iterMax, iterMaxILS, perturbationLevelMax,
			h_weights, h_dest, h_numedges, h_offset, n, m, threadsCount,
			vnd->isFirstImprovementOnOneNeig(), CStar, totalIter, totalTimeSpentOnConstruction, timeSpentInILS,
			constructivePhaseResults, iterationResults);
	timeSpentOnLocalSearch = timeSpentInILS - totalTimeSpentOnConstruction;

	// h_mycluster and h_functionValue
	Imbalance bestValue = CStar.getImbalance();
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentOnBestSolution = (end_time.wall - start_time.wall) / double(1000000000);

	int iterationValue = 0;
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
			<< "," << bestValue.getPositiveValue()
			<< "," << bestValue.getNegativeValue()
			<< setprecision(0)
			<< "," << CStar.getNumberOfClusters()
			<< "," << (iterationValue+1)
			<< "," << fixed << setprecision(4) << timeSpentOnBestSolution
			<< "," << iterMax
			<< "," << numberOfTestedCombinations << endl;

	if(CStar.getImbalance().getValue() < -1.0) {
		BOOST_LOG_TRIVIAL(error) << "An error occurred during CUDA ILS procedure. Is the graph too big for ILS auxiliary matrices?";
		BOOST_LOG_TRIVIAL(error) << "CUDA ILS requires at least 3 x (n x numberOfClusters) bytes of available GPU global memory to run.";
		std::cerr << "An error occurred during CUDA ILS procedure. Is the graph too big for ILS auxiliary matrices?\n";
		std:cerr << "CUDA ILS requires at least 3 x (n x numberOfClusters) bytes of available GPU global memory to run.\n";
	} else {
		BOOST_LOG_TRIVIAL(info) << "ILS procedure done. Obj = " << fixed << setprecision(2) << bestValue.getValue();
	}
	// CStar.printClustering();
	CStar.printClustering(iterationResults, g->getN());
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);
	// saves the best result to output time file
	measureTimeResults(0.0, iterMax);
	generateOutputFile(problem, timeResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timeIntervals"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions"), construct->getAlpha(), vnd->getNeighborhoodSize(), iterMax);

	return CStar;
}

} /* namespace ils */
} /* namespace resolution */
