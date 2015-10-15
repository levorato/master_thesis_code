//============================================================================
// Name        : TestController.h (test suites)
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2015
// Description : Test controller class.
//============================================================================

#ifndef TestController_H_
#define TestController_H_

#include <string>
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>

#include "graph/include/Graph.h"

using namespace std;
namespace fs = boost::filesystem;

namespace controller {

typedef boost::error_info<struct tag_stack_str,std::string> stack_info;
// undirectedS or bidirectionalS
typedef adjacency_list< setS, vecS, undirectedS, property<edge_index_t, int>, property<edge_index_t, int>, vecS > UDGraph;
typedef subgraph< UDGraph > UndirectedSubGraph;
typedef subgraph< UDGraph > UndirectedGraph;


class TestController {
public:
	TestController();
	virtual ~TestController();

	enum StategyName {GRASP, ILS};
	enum SearchName {SEQUENTIAL_SEARCH, PARALLEL_SEARCH};

	string getTimeAndDateAsString();
	int runTestSuite();

	static void terminateMPIProcessesIfAny(int np, int machineProcessAllocationStrategy,
				int numberOfMasters, int numberOfSearchSlavesPerMaster);

private:
	void testSubgraphCreationUndirectedGraph(string executionId, unsigned long seed);
	void testSubgraphCreationPerformance(string executionId, unsigned long seed);
	void testSplitGraphParallelILS(string executionId, unsigned long seed);
	unsigned long mix(unsigned long a, unsigned long b, unsigned long c);

	static void handler();
	UndirectedGraph readGraphFromFile(const string& filepath);
	UndirectedGraph readGraphFromString(const string& graphContents);
	std::string get_file_contents(const char *filename);

};
} /* namespace controller */

#endif /* TestController_H_ */
