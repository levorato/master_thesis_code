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

	string getTimeAndDateAsString();
	int processArgumentsAndExecute(int argc, char *argv[]);

private:
	void processInputFile(fs::path filePath, string& outputFolder, string& timestamp,
			const bool& debug, const int& numberOfIterations, const long& timeLimit,
			const unsigned long& seed);

	void readPropertiesFile();

	void initLogging(string executionId);

	static void handler();

	string logSeverity;
};
} /* namespace controller */

#endif /* COMMANDLINEINTERFACECONTROLLER_H_ */
