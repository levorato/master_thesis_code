/*
 * CUDAGrasp.cpp
 *
 *  Created on: 4/2/2015
 *      Author: Mario Levorato
 */

#include "include/CUDAGrasp.h"
#include "problem/include/ClusteringProblem.h"
#include "graph/include/Imbalance.h"
#include "util/include/RandomUtil.h"
#include "graph/include/Graph.h"
#include "graph/include/SequentialNeighborhoodSearch.h"

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
#include "../vnd/include/CUDASearch.h"
#include "util/include/RandomUtil.h"
#include <limits>
#include <cstdio>
#include <thrust/host_vector.h>
#include <vector>

#include "include/Grasp.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include "graph/include/Imbalance.h"

using namespace problem;
using namespace boost;
using namespace clusteringgraph;
using namespace thrust;
using namespace util;

namespace resolution {
namespace grasp {

CUDAGrasp::CUDAGrasp() :
		Grasp() {
	// TODO Auto-generated constructor stub

}

CUDAGrasp::~CUDAGrasp() {
	// TODO Auto-generated destructor stub
}

Clustering CUDAGrasp::executeGRASP(ConstructClustering &construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, ClusteringProblem& problem, ExecutionInfo& info) {
	BOOST_LOG_TRIVIAL(info)<< "Initializing " << " CUDA GRASP "<< problem.getName() <<
	" procedure for alpha = " << construct.getAlpha() << " and l = " << vnd->getNeighborhoodSize();

	int iterationValue = 0;
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
	for(int edge = 0; i < n; i++) {  // For each vertex i
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
	unsigned short threadsCount = 256;  // limited by shared memory size
	// Pass raw array and its size to kernel
	int totalIter = 0;
	Clustering CStar;
	double totalTimeSpentOnConstruction = 0, timeSpentOnLocalSearch = 0;
	runGRASPKernel(problem, construct, g, info.processRank, vnd->getTimeLimit(), iter,
			h_weights, h_dest, h_numedges, h_offset, n, m, threadsCount,
			vnd->isFirstImprovementOnOneNeig(), CStar, totalIter, totalTimeSpentOnConstruction, timeSpentInGRASP);
	timeSpentOnLocalSearch = timeSpentInGRASP - totalTimeSpentOnConstruction;

	// h_mycluster and h_functionValue
	// TODO capture CUDA GRASP kernel results
	Imbalance bestValue = CStar.getImbalance();
	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentOnBestSolution = (end_time.wall - start_time.wall) / double(1000000000);

	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
	<< "," << bestValue.getPositiveValue()
	<< "," << bestValue.getNegativeValue()
	<< setprecision(0)
	<< "," << CStar.getNumberOfClusters()
	<< "," << (iterationValue+1)
	<< "," << fixed << setprecision(4) << timeSpentOnBestSolution
	<< "," << (totalIter+1)
	<< "," << numberOfTestedCombinations
	<< "," << fixed << setprecision(2) << 0.0 << endl;

	BOOST_LOG_TRIVIAL(info) << "GRASP procedure done. Obj = " << fixed << setprecision(2) << bestValue.getValue();
	BOOST_LOG_TRIVIAL(info) << "Time spent on construction phase: " << fixed << setprecision(2) << totalTimeSpentOnConstruction << "s, " << (100 * totalTimeSpentOnConstruction / timeSpentInGRASP) << "%.";
	BOOST_LOG_TRIVIAL(info) << "Time spent on local search: " << fixed << setprecision(2) << timeSpentOnLocalSearch << "s, " << (100 * timeSpentOnLocalSearch / timeSpentInGRASP) << "%.";
	// CStar.printClustering();
	CStar.printClustering(iterationResults, g->getN());
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations"), construct.getAlpha(), vnd->getNeighborhoodSize(), iter);
	// saves the best result to output time file
	measureTimeResults(0.0, totalIter);
	generateOutputFile(problem, timeResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timeIntervals"), construct.getAlpha(), vnd->getNeighborhoodSize(), iter);
	// saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions"), construct.getAlpha(), vnd->getNeighborhoodSize(), iter);

	return CStar;
}


} /* namespace grasp */
} /* namespace resolution */


