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

	enum StategyName {GRASP, GRASP_PR};
	static string getTimeAndDateAsString();
	static int processArgumentsAndExecute(int argc, char *argv[],
			const int &myRank, const int &np);

private:
	static void processInputFile(fs::path filePath, bool debug, float alpha, int l,
			int numberOfIterations, int np, int myRank);

	static void handler();
};
} /* namespace controller */

#endif /* COMMANDLINEINTERFACECONTROLLER_H_ */
