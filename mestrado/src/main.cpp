//============================================================================
// Name        : main.cpp
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Program entry point. Includes main function and CLI.
//============================================================================

#include "controller/include/CommandLineInterfaceController.h"

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;
using namespace controller;

int main(int ac, char* av[])
{
	// Inicializacao do MPI
	mpi::environment env(ac, av);
	mpi::communicator world;

	CommandLineInterfaceController controller;
	int return_value = controller.processArgumentsAndExecute(ac, av, world.rank(), world.size());

	return return_value;
}
