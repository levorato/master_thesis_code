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

namespace controller {

class CommandLineInterfaceController {
public:
	CommandLineInterfaceController();
	virtual ~CommandLineInterfaceController();

	enum StategyName {GRASP, GRASP_PR};
	static int processArguments(int argc, char *argv[]);
};
} /* namespace controller */

#endif /* COMMANDLINEINTERFACECONTROLLER_H_ */
