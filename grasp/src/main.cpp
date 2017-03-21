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

#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>

#include "graph/include/ParallelBGLSignedGraph.h"

using namespace controller;
using namespace boost;
using boost::graph::distributed::mpi_process_group;


int main(int ac, char* av[])
{
	// Inicializacao do MPI
	boost::mpi::environment env(ac, av);
	boost::mpi::communicator world;

	int rank = process_id(mpi_process_group());
	bool i_am_root = rank == 0;

	std::cout << "Creating distributed graph...\n";
	clusteringgraph::ParallelGraph pgraph;
	std::cout << "Build successfully.\n";

	// cout << "I am rank " << world.rank() << " of " << (world.size()-1) << " running on machine " << env.processor_name() << endl;

	CommandLineInterfaceController controller;
	int return_value = controller.processArgumentsAndExecute(ac, av, world.rank(), world.size(), &pgraph);

	return return_value;
}
