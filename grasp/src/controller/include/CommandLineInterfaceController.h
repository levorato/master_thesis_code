//============================================================================
// Name        : CommandLineInterfaceController.h
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Command Line Interface (CLI) class. Processes program
//               arguments and invokes problem solving.
//============================================================================

#ifndef COMMANDLINEINTERFACECONTROLLER_H_
#define COMMANDLINEINTERFACECONTROLLER_H_

#include <string>
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>

using namespace std;
namespace fs = boost::filesystem;

namespace controller {

typedef boost::error_info<struct tag_stack_str,std::string> stack_info;

class CommandLineInterfaceController {
public:
	CommandLineInterfaceController();
	virtual ~CommandLineInterfaceController();

	enum StategyName {GRASP, ILS};
	string getTimeAndDateAsString();
	int processArgumentsAndExecute(int argc, char *argv[],
			const unsigned int &myRank, const int &np);

private:
	void processInputFile(fs::path filePath, string& outputFolder, string& timestamp,
			const bool& debug, const double& alpha, const int& l, const bool& firstImprovementOnOneNeig,
			const int& numberOfIterations, const long& timeLimit, const int& machineProcessAllocationStrategy,
			const int& numberOfSlaves, const int& numberOfSearchSlaves, const int& myRank,
			const int& functionType, const unsigned long& seed,	const bool& CCEnabled,
			const bool& RCCEnabled, long k, const StategyName& resolutionStrategy,
			const int& iterMaxILS, const int& perturbationLevelMax);

	void readPropertiesFile();

	void initLogging(string executionId, int myRank);

	static void handler();

	string logSeverity;
};
} /* namespace controller */

#endif /* COMMANDLINEINTERFACECONTROLLER_H_ */
