//============================================================================
// Name        : main.cpp
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Program entry point. Includes main function and CLI.
//============================================================================

#include "controller/include/CommandLineInterfaceController.h"

int main(int ac, char* av[])
{
	return controller::CommandLineInterfaceController::processArgumentsAndExecute(ac, av);
}
