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


int main(int argc, char* argv[])
{
	// Inicializacao do MPI
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	// captures the number of vertices as command line argument
    std::string destination;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-v") || (arg == "--vertices")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                destination = argv[i+1]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--vertices option requires one argument." << std::endl;
                return 1;
            }
            break;
        }
    }
    istringstream ss(destination);
    int num_vertices = 0;
    if (!(ss >> num_vertices)) {
        cerr << "Invalid number: " << argv[1] << '\n';
        std::cerr << "--vertices option requires one argument." << std::endl;
        return 1;
    }

	std::cout << "Creating distributed graph with " << num_vertices << " vertices...\n";
	clusteringgraph::ParallelGraph pgraph(num_vertices);
	std::cout << "Graph built successfully.\n";

	// cout << "I am rank " << world.rank() << " of " << (world.size()-1) << " running on machine " << env.processor_name() << endl;

	CommandLineInterfaceController controller;
	int return_value = controller.processArgumentsAndExecute(argc, argv, world.rank(), world.size(), &pgraph);

	return return_value;
}
