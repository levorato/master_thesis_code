//============================================================================
// Name        : main.cpp
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Program entry point. Includes main function and CLI.
//============================================================================

#include "controller/include/CommandLineInterfaceController.h"

#include <mpi.h>

int main(int ac, char* av[])
{
	MPI_Status status;
	int myRank;

	// Inicializacao do MPI
	MPI_Init(&ac, &av);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	int return_value = controller::CommandLineInterfaceController::processArgumentsAndExecute(ac, av, myRank);

	// Finalizacao do MPI
	MPI_Finalize();

	return return_value;
}
