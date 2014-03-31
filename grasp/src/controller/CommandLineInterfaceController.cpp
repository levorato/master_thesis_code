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
#include "../graph/include/Imbalance.h"
#include "../problem/include/ClusteringProblem.h"
#include "../problem/include/CCProblem.h"
#include "../util/include/TimeDateUtil.h"
#include "../util/include/EnumUtil.h"
#include "../util/include/MPIMessage.h"
#include "../util/parallel/include/MPIUtil.h"
#include "../problem/include/ClusteringProblemFactory.h"
#include "../resolution/grasp/include/GainFunctionFactory.h"
#include "../resolution/grasp/include/GainFunction.h"
#include "../graph/include/NeighborhoodSearchFactory.h"
#include "../graph/include/ParallelNeighborhoodSearch.h"
#include "../graph/include/SequentialNeighborhoodSearch.h"
#include "../resolution/vnd/include/RCCVariableNeighborhoodSearch.h"

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/support/date_time.hpp>

using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace boost::mpi;

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <execinfo.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <vector>
using namespace std;
using namespace clusteringgraph;
using namespace resolution::grasp;
using namespace resolution::vnd;
using namespace problem;
using namespace util;
using namespace util::parallel;

namespace controller {

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const std::vector<T>& v)
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


CommandLineInterfaceController::CommandLineInterfaceController() : logSeverity("warning") {
	// TODO Auto-generated constructor stub

}

CommandLineInterfaceController::~CommandLineInterfaceController() {
	// TODO Auto-generated destructor stub
}

void CommandLineInterfaceController::processInputFile(fs::path filePath, string& outputFolder,
		string& executionId, const bool& debug, const double& alpha, const int& l, const bool& firstImprovementOnOneNeig,
		const int& numberOfIterations, const long& timeLimit, const int& numberOfSlaves, const int& numberOfSearchSlaves,
		const int& myRank, const int& problemType, const int& functionType, const unsigned long& seed,
		const bool& RCCEnabled) {
	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g = reader.readGraphFromFile(filePath.string());

		Clustering c;
		string fileId = filePath.filename().string();
		ClusteringProblemFactory problemFactory;

		// -------------------  C C    G R A S P     P R O C E S S I N G -------------------------
		GainFunctionFactory functionFactory(g.get());
		// medicao de tempo do CC
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		if(numberOfSlaves == 0) {	// sequential version of GRASP
			Grasp resolution(functionFactory.build(functionType), seed);
			c = resolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
					problemFactory.build(problemType), executionId, fileId, outputFolder,
					timeLimit, numberOfSlaves, myRank, numberOfSearchSlaves);
		} else {  // parallel version
			// distributes GRASP processing among numberOfSlaves processes and summarizes the result
			ParallelGrasp parallelResolution(functionFactory.build(functionType), seed);
			c = parallelResolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
					problemFactory.build(problemType), executionId, fileId, outputFolder, timeLimit,
					numberOfSlaves, myRank, numberOfSearchSlaves);
		}

		// Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		double timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
		// Saves elapsed time and best solution to output file
		string filename = outputFolder + "/" + fileId + "/" + executionId + "/cc-result.txt";
		ofstream out(filename.c_str(), ios::out | ios::trunc); 
		if(!out) { 
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << filename;
		    // TODO tratar excecao 
  		} 

		BOOST_LOG_TRIVIAL(info) << "Global time spent: " << timeSpent << " s";
		out << "Global time spent: " << timeSpent << endl;
		Imbalance imb = c.getImbalance();
		out << "I(P) = " << imb.getValue() << endl;
		stringstream ss;
		c.printClustering(ss);
		out << ss.str();
 		out.close();

 		// -------------------  R C C    G R A S P     P R O C E S S I N G -------------------------
 		// 		TODO: Uses GRASP best result (cluster c) - number of clusters as input (k)
 		if(RCCEnabled) {
 			BOOST_LOG_TRIVIAL(info) << "Starting RCC GRASP..." << endl;
			// medicao de tempo do RCC
			timer.start();
			start_time = timer.elapsed();
			Clustering RCCCluster;

			GainFunctionFactory functionFactory(g.get());
			if(numberOfSlaves == 0) {	// sequential version of GRASP
				Grasp resolution(functionFactory.build(functionType), seed);
				RCCCluster = resolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
						problemFactory.build(ClusteringProblem::RCC_PROBLEM), executionId, fileId, outputFolder,
						timeLimit, numberOfSlaves, myRank, numberOfSearchSlaves);
			} else {  // parallel version
				// distributes GRASP processing among numberOfSlaves processes and summarizes the result
				ParallelGrasp parallelResolution(functionFactory.build(functionType), seed);
				RCCCluster = parallelResolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
						problemFactory.build(ClusteringProblem::RCC_PROBLEM), executionId, fileId, outputFolder, timeLimit,
						numberOfSlaves, myRank, numberOfSearchSlaves);
			}
			/*
			Clustering RCCCluster = vns.executeSearch(g.get(), c, l, c.getNumberOfClusters(), false,
									problemFactory.build(ClusteringProblem::RCC_PROBLEM), *neig, executionId, fileId, outputFolder,
									timeLimit, numberOfSlaves, myRank, numberOfSearchSlaves);
			*/
			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			// Saves elapsed time and best solution to output file
			filename = outputFolder + "/" + fileId + "/" + executionId + "/rcc-result.txt";
			ofstream out(filename.c_str(), ios::out | ios::trunc);
			if(!out) {
				BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << filename;
				// TODO tratar excecao
			}

			BOOST_LOG_TRIVIAL(info) << "Global time spent: " << timeSpent << " s";
			out << "RCC Global time spent: " << timeSpent << endl;
			Imbalance imb = RCCCluster.getImbalance();
			out << "SRI(P) = " << imb.getValue() << endl;
			stringstream ss;
			RCCCluster.printClustering(ss);
			out << ss.str();
			out.close();
			BOOST_LOG_TRIVIAL(info) << "RCC Search done. Obj = " << imb.getValue();
 		}

	} else {
		BOOST_LOG_TRIVIAL(fatal) << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
	}
}

// http://www.concentric.net/~Ttwang/tech/inthash.htm
unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

int CommandLineInterfaceController::processArgumentsAndExecute(int argc, char *argv[],
		const unsigned int &myRank, const int &np) {

	// used for debugging purpose
	std::set_terminate( handler );

	mpi::communicator world;
	cout << "Correlation clustering problem solver" << endl;

	// random seed used in grasp
	/*
	* Caveat: std::time(0) is not a very good truly-random seed.  When
	* called in rapid succession, it could return the same values, and
	* thus the same random number sequences could ensue.
	* Instead, we are using boost::random_device
	* http://stackoverflow.com/questions/4329284/c-boost-random-numeric-generation-problem
	* http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
	*/
	unsigned long seed = mix(clock(), time(NULL), getpid());
	// boost::minstd_rand generator(seed);

	// reads the system properties from ini file
	this->readPropertiesFile();

	string s_alpha;
	int numberOfIterations = 500, l = 1;
	bool debug = false, profile = false, firstImprovementOnOneNeig = true, RCCEnabled = true;
	string inputFileDir, outputFolder;
	int timeLimit = 1800;
	int problemType = ClusteringProblem::CC_PROBLEM, functionType = GainFunction::IMBALANCE;
	int numberOfSlaves = np - 1;
	string jobid;
	CommandLineInterfaceController::StategyName strategy = CommandLineInterfaceController::GRASP;

	po::options_description desc("Available options:");
	desc.add_options()
		("help", "show program options")
		("alpha,a", po::value<string>(&s_alpha),
			  "alpha - randomness factor")
		("iterations,it", po::value<int>(&numberOfIterations)->default_value(500),
			  "number of iterations")
		("neighborhood_size,l", po::value<int>(&l)->default_value(1),
			  "neighborhood size")
		("rcc", po::value<bool>(&RCCEnabled)->default_value(true), "Enable RCC Problem resolution after GRASP CC")
		("time-limit", po::value<int>(&timeLimit)->default_value(1800), "maximum execution time")
		("input-file", po::value< std::vector<string> >(), "input file")
		("debug", po::value<bool>(&debug)->default_value(false), "enable debug mode")
		("profile", po::value<bool>(&profile)->default_value(false), "enable profile mode")
		("input-file-dir", po::value<string>(&inputFileDir), "input file directory (processes all files inside)")
		("output-folder", po::value<string>(&outputFolder), "output folder for results files")
		("gain-function-type", po::value<int>(&functionType),
				"0 for min imbalance, 1 for max modularity gain function, 2 for max negative modularity gain function, "
				"3 for max positive-negative modularity gain function, 4 for max pos-neg mod gain function II, "
				"5 for max pos-neg mod gain function III")
		("slaves,s", po::value<int>(&numberOfSlaves)->default_value(np - 1), "number of GRASP processes in parallel")
		("firstImprovementOnOneNeig", po::value<bool>(&firstImprovementOnOneNeig)->default_value(true), "first improvement in 1-opt neighborhood (** sequential VNS only **)")
		("jobid", po::value<string>(&jobid), "job identifier (string)")
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

	// job id is obtained through command line parameter from PBS Scheduler
	cout << "Job id is " << jobid << "\n";
	// initializes the logging subsystem
	if(jobid.length() == 0) {
		jobid = TimeDateUtil::generateRandomId();
	} else {
		jobid += "-" + TimeDateUtil::getDateAndTime();
	}
	
	try {
        this->initLogging(jobid, myRank);
	}
	catch(std::exception& e)
	{
	     cerr << "Fatal application error in log init.\n";
		 cerr << "Abnormal program termination. Stracktrace: " << endl;
		 cerr << e.what() << "\n";
		 if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
			 cerr << stack << endl;
		 }
		 cerr << diagnostic_information(e);
		 return 1;
	 }

	// codigo do processor lider
	if(myRank == 0) {
		cout << "Correlation clustering problem solver" << endl;
		// id used for output folders
		string executionId = jobid;

		try {
			if (vm.count("help")) {
				cout << "Usage: MestradoMario [options]\n";
				cout << desc;
				return 0;
			}

			if(not vm.count("output-folder")) {
				cout << "Please specify the output folder." << endl;
				return 1;
			}

			float alpha = 0.0;
			if(vm.count("alpha")) {
				// std::istringstream i(s_alpha);
				// if (!(i >> alpha))
				//     throw BadConversion("convertToDouble(\"" + s + "\")");
				sscanf(s_alpha.c_str(), "%f", &alpha);				
			}

			BOOST_LOG_TRIVIAL(info) << "Resolution strategy is " << strategy << endl;
			BOOST_LOG_TRIVIAL(info) << "Neighborhood size (l) is " << l << "\n";
			BOOST_LOG_TRIVIAL(info) << "Gain function type is ";
			if(functionType == GainFunction::MODULARITY) {
				BOOST_LOG_TRIVIAL(info) << "max modularity\n";
			} else if(functionType == GainFunction::NEGATIVE_MODULARITY) {
				BOOST_LOG_TRIVIAL(info) << "max negative modularity\n";
			} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY) {
				BOOST_LOG_TRIVIAL(info) << "max positive-negative modularity\n";
			} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY_II) {
				BOOST_LOG_TRIVIAL(info) << "max positive-negative modularity II\n";
			} else if(functionType == GainFunction::POSITIVE_NEGATIVE_MODULARITY_III) {
				BOOST_LOG_TRIVIAL(info) << "max positive-negative modularity III\n";
			} else {
				BOOST_LOG_TRIVIAL(info) << "min imbalance\n";
			}
			BOOST_LOG_TRIVIAL(info) << "First improvement on 1-opt neighborhood is enabled? " << firstImprovementOnOneNeig;
			if (RCCEnabled) {
				BOOST_LOG_TRIVIAL(info) << "RCC is enabled. Will use GRASP CC results." << endl;
			} else {
				BOOST_LOG_TRIVIAL(info) << "RCC is disabled. Only GRASP CC will be executed." << endl;
			}
			if(alpha <= 0.0) {
				BOOST_LOG_TRIVIAL(info) << "VOTE is enabled. Will always choose best-gain vertex on constructionPhase." << endl;
			}

			BOOST_LOG_TRIVIAL(info) << "Total number of processes is " << np << endl;
            cout << "Total number of processes is " << np << endl;
			BOOST_LOG_TRIVIAL(info) << "Number of GRASP processes (masters) in parallel is " << (numberOfSlaves + 1);
			cout << "Number of GRASP processes (masters) in parallel is " << (numberOfSlaves + 1) << endl;
			BOOST_LOG_TRIVIAL(info) << "Number of VNS processes in parallel is " << (np - numberOfSlaves - 1) << endl;
            cout << "Number of VNS processes in parallel is " << (np - numberOfSlaves - 1) << endl;
			std::vector<fs::path> fileList;

			if (vm.count("input-file")) {
				fs::path filePath (vm["input-file"].as< std::vector<string> >().at(0));
				BOOST_LOG_TRIVIAL(info) << "Input file is: "
									 << filePath.string() << "\n";
				fileList.push_back(filePath);
			} else if(vm.count("input-file-dir")) {
				BOOST_LOG_TRIVIAL(info) << "Input file dir is: " << inputFileDir << endl;
				fs::path inputDir(inputFileDir);
				fs::directory_iterator end_iter;
				if (!fs::exists(inputDir) || !fs::is_directory(inputDir)) {
					BOOST_LOG_TRIVIAL(fatal) << "Input file directory not found. Exiting." << endl;
					return 1;
				}
				for( fs::directory_iterator dir_iter(inputDir) ; dir_iter != end_iter ; ++dir_iter) {
					string ext = dir_iter->path().extension().string();
					boost::algorithm::to_lower(ext);
					if ((fs::is_regular_file(dir_iter->status())) &&
							((ext == ".g") || (ext == ".net") || (ext == ".dat"))) {
						fs::path filePath = *dir_iter;
						fileList.push_back(filePath);
					}
				}
			} else {
				BOOST_LOG_TRIVIAL(fatal) << "Please specify at least one input file.";
				return 1;
			}

			unsigned int numberOfSearchSlaves = MPIUtil::calculateNumberOfSearchSlaves(np, numberOfSlaves);
			if(np > 1) {
				// broadcasts the numberOfSlaves to all processes
				// mpi::broadcast(world, numberOfSlaves, 0);
				for(int i = 1; i < np; i++) {
					world.send(i, MPIMessage::INPUT_MSG_NUM_SLAVES_TAG, numberOfSlaves);
				}
			}

			// -------------------  G R A P H     F I L E S     P R O C E S S I N G -------------------------
			for(unsigned int i = 0; i < fileList.size(); i++) {
				fs::path filePath = fileList.at(i);

				if(not profile) { // default mode: use specified parameters
					BOOST_LOG_TRIVIAL(info) << "Alpha value is " << std::setprecision(2) << fixed << alpha << "\n";
					BOOST_LOG_TRIVIAL(info) << "Number of iterations is " << numberOfIterations << "\n";
					processInputFile(filePath, outputFolder, executionId, debug, alpha, l, firstImprovementOnOneNeig,
							numberOfIterations, timeLimit, numberOfSlaves, numberOfSearchSlaves,
							myRank, problemType, functionType, seed, RCCEnabled);
				} else {
					BOOST_LOG_TRIVIAL(info) << "Profile mode on." << endl;
					for(double alpha2 = 0.0F; alpha2 < 1.1F; alpha2 += 0.1F) {
						BOOST_LOG_TRIVIAL(info) << "Processing GRASP with alpha = " << std::setprecision(2) << alpha2 << endl;
						processInputFile(filePath, outputFolder, executionId, debug, alpha2, l, firstImprovementOnOneNeig,
								numberOfIterations, timeLimit, numberOfSlaves, numberOfSearchSlaves,
								myRank, problemType, functionType, seed, RCCEnabled);
					}
				}
			}

			// ------------------ M P I    T E R M I N A T I O N ---------------------
			BOOST_LOG_TRIVIAL(info) << "Terminating MPI slave processes...";
			if(np > 1) {
				InputMessageParallelGrasp imsgpgrasp;
                                InputMessageParallelVNS imsgpvns;
				for(int i = 1; i < np; i++) {
					if(MPIUtil::isGRASPSlave(i, numberOfSlaves, numberOfSearchSlaves)) {  // GRASP slave
						world.send(i, MPIMessage::TERMINATE_MSG_TAG, imsgpgrasp);
					} else {  // VNS slave
						world.send(i, MPIMessage::TERMINATE_MSG_TAG, imsgpvns);
					}
					BOOST_LOG_TRIVIAL(debug) << "Terminate message sent to process " << i << endl;
				}
			}
			BOOST_LOG_TRIVIAL(info) << "Done.";
		}
		catch(std::exception& e)
		{
			cerr << "Fatal application error.\n";
			BOOST_LOG_TRIVIAL(fatal) << "Abnormal program termination. Stracktrace: " << endl;
			BOOST_LOG_TRIVIAL(fatal) << e.what() << "\n";
			if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
				BOOST_LOG_TRIVIAL(fatal) << stack << endl;
			}
			BOOST_LOG_TRIVIAL(fatal) << diagnostic_information(e);
			return 1;
		}
	} else { // slave processes
		// First message received from leader contains the number of GRASP slaves
		// and is used for a process to discover if it is a GRASP slave or a VNS slave
		// Grasp slave: rank % (numberOfSearchSlaves+1) == 0; VNS slave: otherwise
		unsigned int numberOfSlaves = 0;
		world.recv(MPIMessage::LEADER_ID, mpi::any_tag, numberOfSlaves);
		BOOST_LOG_TRIVIAL(info) << "number of slaves received: " << numberOfSlaves << "\n";
		cout << "number of slaves received: " << numberOfSlaves << "\n";
		cout << "Process " << myRank << " ready [slave process]." << endl;
		unsigned int numberOfSearchSlaves = MPIUtil::calculateNumberOfSearchSlaves(np, numberOfSlaves);

		// Controls the number of messages received
		unsigned int messageCount = 0;
		// common slave variables
		ClusteringProblemFactory problemFactory;
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g;
		unsigned int previousId = 0;

		try {
			if(MPIUtil::isGRASPSlave(myRank, numberOfSlaves, numberOfSearchSlaves)) {  // GRASP slave
				BOOST_LOG_TRIVIAL(info) << "Process " << myRank << " ready [GRASP slave process].";
				while(true) {
					// Receives a message with GRASP parameters and triggers local GRASP execution
					InputMessageParallelGrasp imsgpg;
					// receives a message of type ParallelGrasp::MPIMessage::INPUT_MSG_PARALLEL_GRASP_TAG or a terminate msg
					mpi::status stat = world.recv(MPIMessage::LEADER_ID, mpi::any_tag, imsgpg);
					if(stat.tag() == MPIMessage::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": terminate message received.";
						return 0;
					}
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel GRASP]: Received message from leader.";
					messageCount++;

					if(imsgpg.l == 0) {
						BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty GRASP message received. Terminating.\n";
						return 1;
					}
					// builds a new graph object only if it has changed (saves time)
					if(previousId != imsgpg.id) {
						g.reset();
						// Reads the graph from the specified text file
						g = reader.readGraphFromFile(imsgpg.graphInputFilePath);
						previousId = imsgpg.id;
					}

					// trigggers the local GRASP routine
					GainFunctionFactory functionFactory(g.get());
					Grasp resolution(functionFactory.build(imsgpg.gainFunctionType), seed);
					Clustering bestClustering = resolution.executeGRASP(g.get(), imsgpg.iter, imsgpg.alpha,
							imsgpg.l, imsgpg.firstImprovementOnOneNeig, problemFactory.build(imsgpg.problemType),
							imsgpg.executionId, imsgpg.fileId, imsgpg.outputFolder, imsgpg.timeLimit,
							imsgpg.numberOfSlaves, myRank, imsgpg.numberOfSearchSlaves);

					// Sends the result back to the leader process
					OutputMessage omsg(bestClustering, resolution.getNumberOfTestedCombinations());
					world.send(MPIMessage::LEADER_ID, MPIMessage::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": GRASP Output Message sent to leader.";

				}
			} else {  // Parallel VNS slave
				BOOST_LOG_TRIVIAL(info) << "Process " << myRank << " ready [VNS slave process].";
				while(true) {
					// Receives a parallel VNS message with parameters and triggers local VNS execution
					InputMessageParallelVNS imsgvns;
					// receives a message of type ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG or a terminate msg
					mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, imsgvns);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel VNS]: Received message from leader." << endl;
					messageCount++;
					if(stat.tag() == MPIMessage::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": terminate msg received.\n";
						return 0;
					} else if(stat.tag() == MPIMessage::INTERRUPT_MSG_PARALLEL_VNS_TAG) {
						// ignores interrupt message from previous VNS processing
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": ignoring VNS interrupt msg received.\n";
						continue;
					}

					if(imsgvns.l == 0) {
						BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty VNS message received. Terminating.\n";
						return 1;
					}
					// builds a new graph object only if it has changed (saves time)
					if(previousId != imsgvns.id) {
						g.reset();
						g = reader.readGraphFromFile(imsgvns.graphInputFilePath);
						previousId = imsgvns.id;
					}

					// triggers the local partial VNS search to be done between initial and final cluster indices
					ParallelNeighborhoodSearch pnSearch(numberOfSlaves, numberOfSearchSlaves);
					Clustering bestClustering = pnSearch.searchNeighborhood(imsgvns.l, g.get(), &imsgvns.clustering,
							problemFactory.build(imsgvns.problemType), imsgvns.timeSpentSoFar, imsgvns.timeLimit, seed, myRank,
							imsgvns.initialClusterIndex, imsgvns.finalClusterIndex, false, imsgvns.k);

					// Sends the result back to the leader process
					int leader = stat.source();
					OutputMessage omsg(bestClustering, pnSearch.getNumberOfTestedCombinations());
					world.send(leader, MPIMessage::OUTPUT_MSG_PARALLEL_VNS_TAG, omsg);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": VNS Output Message sent to grasp leader number "
							<< leader << ". I(P) = " << bestClustering.getImbalance().getValue();
				}
			}
		} catch(std::exception& e) {
			BOOST_LOG_TRIVIAL(fatal) << "Abnormal program termination. Stracktrace: " << endl;
			BOOST_LOG_TRIVIAL(fatal) << e.what() << "\n";
			if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
				BOOST_LOG_TRIVIAL(fatal) << stack << endl;
			}
			BOOST_LOG_TRIVIAL(fatal) << diagnostic_information(e);
			return 1;
		}
	}
	return 0;
}

void CommandLineInterfaceController::readPropertiesFile() {
	namespace pt = boost::property_tree;
	pt::ptree propTree;

	if ( !boost::filesystem::exists( "config.ini" ) ) {
	  cout << "Can't find config.ini file! Assuming default properties." << endl;
	} else {
		read_ini("config.ini", propTree);

		boost::optional<string> severity = propTree.get_optional<std::string>("logging.severity");
		if(severity) {
			logSeverity = *severity;
		} else {
			cout << "WARNING: warn log level not specified, assuming warn level." << endl;
		}
	}
}

void CommandLineInterfaceController::initLogging(string executionId, int myRank) {
	using namespace boost::log;
	namespace keywords = boost::log::keywords;
	namespace sinks = boost::log::sinks;
	namespace expr = boost::log::expressions;
	namespace attrs = boost::log::attributes;
	string filename = string("logs/") + executionId + string("-Node") + lexical_cast<string>(myRank) + string(".log");
	cout << "Init logging..." << endl;

	LogSeverityEnumParser parser;
	logging::trivial::severity_level mySeverity = parser.ParseSomeEnum(logSeverity);

	boost::shared_ptr<log::core> logger = log::core::get();
	logger->set_logging_enabled( true );
	logger->add_global_attribute("TimeStamp", attrs::local_clock());
	logger->add_global_attribute("ProcessID", attrs::current_process_id());
	logger->add_global_attribute("LineID", attrs::counter< unsigned int >(1));

	// @see boost log bug at: https://github.com/azat/boostcache/commit/839d14fbb7285ba3d702ac5c5d1b4c5adce81706
	logging::add_file_log
	    (
	        keywords::file_name = filename.c_str(),
	        keywords::rotation_size = 10 * 1024 * 1024,
	        keywords::time_based_rotation = sinks::file::rotation_at_time_point(0, 0, 0),
	        keywords::auto_flush = true,
	        keywords::open_mode = (std::ios::out | std::ios::app),
	        keywords::format = (
				expressions::stream
					<< "[" << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%d.%m.%Y %H:%M:%S.%f") << "] "
					<< "[" << expr::attr< attrs::current_process_id::value_type >("ProcessID") << "] "
					<< "(" << expr::attr< unsigned int >("LineID") << ")"
					<< "[" << logging::trivial::severity << "]\t"
					<< expr::smessage
			)
	    );
	// add_common_attributes();

	// This makes the sink to write log records that look like this:
	// 1: <normal> A normal severity message
	// 2: <error> An error severity message
	/*
	sink->set_formatter
	(
		expr::format("%1%: <%2%> %3%")
			% expr::attr< unsigned int >("TimeStamp")
			% logging::trivial::severity
			% expr::smessage
	);*/

	logger->set_filter
    (
        logging::trivial::severity >= mySeverity
    );
}

void CommandLineInterfaceController::handler()
{
    void *trace_elems[20];
    int trace_elem_count(backtrace( trace_elems, 20 ));
    char **stack_syms(backtrace_symbols( trace_elems, trace_elem_count ));
    for ( int i = 0 ; i < trace_elem_count ; ++i )
    {
    	BOOST_LOG_TRIVIAL(fatal) << stack_syms[i] << "\n";
    }
    free( stack_syms );
}

} // namespace controller
