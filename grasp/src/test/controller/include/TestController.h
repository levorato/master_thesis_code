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

using namespace std;
namespace fs = boost::filesystem;

namespace controller {

typedef boost::error_info<struct tag_stack_str,std::string> stack_info;

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
	void testSubgraphCreationPerformance(string executionId, unsigned long seed);
	void testSplitGraphParallelILS(string executionId, unsigned long seed);
	unsigned long mix(unsigned long a, unsigned long b, unsigned long c);

	static void handler();

};
} /* namespace controller */

#endif /* TestController_H_ */
