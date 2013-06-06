/*
 * CommandLineInterfaceController.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: mario
 */

#include "include/CommandLineInterfaceController.h"
#include "include/SimpleTextGraphFileReader.h"
#include "../graph/include/Graph.h"
#include "../resolution/grasp/include/Grasp.h"
#include "../resolution/grasp/include/ParallelGrasp.h"
#include "../graph/include/Clustering.h"
#include "../problem/include/ClusteringProblem.h"
#include "../problem/include/CCProblem.h"
#include "../util/include/TimeDateUtil.h"
#include "../util/include/MPIMessage.h"
#include "../problem/include/ClusteringProblemFactory.h"

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace boost::mpi;

#include <cstdlib>
#include <stdexcept>
#include <execinfo.h>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
using namespace std;
using namespace clusteringgraph;
using namespace resolution::grasp;
using namespace problem;
using namespace util;

namespace controller {

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

std::istream& operator>>(std::istream& in, CommandLineInterfaceController::StategyName& strategy)
{
    std::string token;
    in >> token;
    if (token == "GRASP")
        strategy = CommandLineInterfaceController::GRASP;
    else if (token == "GRASP_PR")
        strategy = CommandLineInterfaceController::GRASP_PR;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}


CommandLineInterfaceController::CommandLineInterfaceController() {
	// TODO Auto-generated constructor stub

}

CommandLineInterfaceController::~CommandLineInterfaceController() {
	// TODO Auto-generated destructor stub
}

void CommandLineInterfaceController::processInputFile(fs::path filePath, string& timestamp,
		bool debug, float alpha, int l, int numberOfIterations, int np, int myRank) {
	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g = reader.readGraphFromFile(filePath.string());
		if (debug) {
			g->printGraph();
		}

		// Triggers the execution of the GRASP algorithm
		CCProblem problem;
		ClusteringPtr c;
		string fileId = filePath.filename().string();

		if(np == 1) {	// sequential version of GRASP
			Grasp resolution;
			resolution.executeGRASP(g.get(), numberOfIterations,
					alpha, l, problem, timestamp, fileId, myRank);
		} else {  // parallel version
			// distributes GRASP processing among the processes and summarizes the result
			ParallelGrasp parallelResolution;
			parallelResolution.executeGRASP(g.get(), numberOfIterations,
					alpha, l, problem, timestamp, fileId, np, myRank);
		}
	} else {
		cerr << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
	}
}

int CommandLineInterfaceController::processArgumentsAndExecute(int argc, char *argv[],
		const int &myRank, const int &np) {

	// used for debugging purpose
	std::set_terminate( handler );
	// timestamp used for output files
	string timestamp = TimeDateUtil::getTimeAndDateAsString();

	// codigo do processor lider
	if(myRank == 0) {
		cout << "Correlation clustering problem solver" << endl << endl;

		try {
			float alpha = 0.5F;
			int numberOfIterations = 500, k = -1, l = 1;
			bool debug = false, RCC = false, profile = false;
			string inputFileDir;
			CommandLineInterfaceController::StategyName strategy = CommandLineInterfaceController::GRASP;

			po::options_description desc("Available options:");
			desc.add_options()
				("help", "show program options")
				("alpha,a", po::value<float>(&alpha)->default_value(0.5F),
					  "alpha - randomness factor")
				("iterations,it", po::value<int>(&numberOfIterations)->default_value(500),
					  "number of iterations")
				("neighborhood_size,l", po::value<int>(&l)->default_value(1),
					  "neighborhood size")
				("k,k", po::value<int>(&k)->default_value(-1), "k parameter (RCC problem)")
				("input-file", po::value< vector<string> >(), "input file")
				("debug", po::value<bool>(&debug)->default_value(false), "enable debug mode")
				("profile", po::value<bool>(&profile)->default_value(false), "enable profile mode")
				("input-file-dir", po::value<string>(&inputFileDir), "input file directory (processes all files inside)")
				/* TODO Resolver problema com o parametro da descricao
				("strategy",
							 po::typed_value<Resolution::StategyName, char *>(&strategy).default_value(strategy, "GRASP"),
							 "Resolution strategy to be used. Accepted values: GRASP, GRASP_PR.") */
			;

			po::positional_options_description p;
			p.add("input-file", -1);

			po::variables_map vm;
			po::store(po::command_line_parser(argc, argv).
					  options(desc).positional(p).run(), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << "Usage: MestradoMario [options]\n";
				cout << desc;
				return 0;
			}

			if (k != -1) {
				cout << "k value is " << k << ". RCC is enabled." << endl;
				RCC = true;
			}

			cout << "Resolution strategy is " << strategy << endl;
			cout << "Neighborhood size (l) is " << l << "\n";
			vector<fs::path> fileList;

			if (vm.count("input-file")) {
				cout << "Input files are: "
					 << vm["input-file"].as< vector<string> >() << "\n";
				fs::path filePath (vm["input-file"].as< vector<string> >().at(0));
				fileList.push_back(filePath);
			} else if(vm.count("input-file-dir")) {
				cout << "Input file dir is: " << inputFileDir << endl;
				fs::path inputDir(inputFileDir);
				fs::directory_iterator end_iter;
				if (!fs::exists(inputDir) || !fs::is_directory(inputDir)) {
					cout << "Input file directory not found. Exiting." << endl;
					return 1;
				}
				for( fs::directory_iterator dir_iter(inputDir) ; dir_iter != end_iter ; ++dir_iter) {
					string ext = dir_iter->path().extension().string();
					boost::algorithm::to_lower(ext);
					if (fs::is_regular_file(dir_iter->status()) &&
							(ext == ".g" || ext == ".net")) {
						fs::path filePath = *dir_iter;
						fileList.push_back(filePath);
					}
				}
			} else {
				cout << "Please specify at least one input file.";
				return 1;
			}

			for(unsigned int i = 0; i < fileList.size(); i++) {
				fs::path filePath = fileList.at(i);

				if(not profile) { // default mode: use specified parameters
					cout << "Alpha value is " << std::setprecision(2) << fixed << alpha << "\n";
					cout << "Number of iterations is " << numberOfIterations << "\n";
					processInputFile(filePath, timestamp, debug, alpha, l, numberOfIterations, np, myRank);
				} else {
					cout << "Profile mode on." << endl;
					for(float alpha2 = 0.0F; alpha2 < 1.1F; alpha2 += 0.1F) {
						cout << "Processing GRASP with alpha = " << std::setprecision(2) << alpha2 << endl;
						processInputFile(filePath, timestamp, debug, alpha2, l, numberOfIterations, np, myRank);
					}
				}
			}
		}
		catch(std::exception& e)
		{
			cout << "Abnormal program termination. Stracktrace: " << endl;
			cout << e.what() << "\n";
			if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
				std::cout << stack << endl;
			}
			std::cerr << diagnostic_information(e);
			return 1;
		}
	} else { // processos seguidores
		try {
			// Receives a message with GRASP parameters and triggers local execution
			InputMessage imsg;
			mpi::communicator world;
			mpi::status stat = world.recv(mpi::any_source, ParallelGrasp::INPUT_MSG_TAG, imsg);

			cout << "Process " << myRank << ": Received message from leader." << endl;
			// reconstructs the graph from its text representation
			SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
			SignedGraphPtr g = reader.readGraphFromString(imsg.graphInputFileContents);

			// trigggers the local GRASP routine
			ClusteringProblemFactory factory;
			Grasp resolution;
			ClusteringPtr bestClustering = resolution.executeGRASP(g.get(), imsg.iter, imsg.alpha,
						imsg.l, factory.build(imsg.problemType), timestamp, imsg.fileId, myRank);

			// Sends the result back to the leader process
			OutputMessage omsg(*bestClustering);
			world.send(ParallelGrasp::LEADER_ID, ParallelGrasp::OUTPUT_MSG_TAG, omsg);
			cout << "Process " << myRank << ": Message sent to leader." << endl;
		}
		catch(std::exception& e)
		{
			cout << "Abnormal program termination. Stracktrace: " << endl;
			cout << e.what() << "\n";
			if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
				std::cout << stack << endl;
			}
			std::cerr << diagnostic_information(e);
			return 1;
		}
	}
	return 0;
}

void CommandLineInterfaceController::handler()
{
    void *trace_elems[20];
    int trace_elem_count(backtrace( trace_elems, 20 ));
    char **stack_syms(backtrace_symbols( trace_elems, trace_elem_count ));
    for ( int i = 0 ; i < trace_elem_count ; ++i )
    {
        std::cout << stack_syms[i] << "\n";
    }
    free( stack_syms );
}

} // namespace controller
