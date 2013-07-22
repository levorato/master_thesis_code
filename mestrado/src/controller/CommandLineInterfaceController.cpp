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
#include "../problem/include/ClusteringProblemFactory.h"
#include "../resolution/grasp/include/GainFunctionFactory.h"
#include "../resolution/grasp/include/GainFunction.h"
#include "../graph/include/NeighborhoodSearchFactory.h"
#include "../graph/include/ParallelNeighborhoodSearch.h"
#include "../graph/include/SequentialNeighborhoodSearch.h"

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


CommandLineInterfaceController::CommandLineInterfaceController() : logSeverity("warning") {
	// TODO Auto-generated constructor stub

}

CommandLineInterfaceController::~CommandLineInterfaceController() {
	// TODO Auto-generated destructor stub
}

void CommandLineInterfaceController::processInputFile(fs::path filePath, string& outputFolder,
		string& executionId, const bool& debug, const double& alpha, const int& l, const bool& firstImprovementOnOneNeig,
		const int& numberOfIterations, const long& timeLimit, const int& numberOfSlaves, const int& numberOfSearchSlaves,
		const int& myRank, const int& problemType, const int& functionType, const unsigned long& seed) {
	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g = reader.readGraphFromFile(filePath.string());
		if (debug) {
			g->printGraph();
		}

		// Triggers the execution of the GRASP algorithm
		ClusteringPtr c;
		string fileId = filePath.filename().string();
		ClusteringProblemFactory problemFactory;
		GainFunctionFactory functionFactory(g.get());
		// medicao de tempo
		boost::timer::cpu_timer timer;
		timer.start();
		boost::timer::cpu_times start_time = timer.elapsed();

		if(numberOfSlaves == 0) {	// sequential version of GRASP
			Grasp resolution(&functionFactory.build(functionType), seed);
			c = resolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
					problemFactory.build(problemType), executionId, fileId, outputFolder,
					timeLimit, numberOfSlaves, myRank, numberOfSearchSlaves);
		} else {  // parallel version
			// distributes GRASP processing among numberOfSlaves processes and summarizes the result
			ParallelGrasp parallelResolution(&functionFactory.build(functionType), seed);
			c = parallelResolution.executeGRASP(g.get(), numberOfIterations, alpha, l, firstImprovementOnOneNeig,
					problemFactory.build(problemType), executionId, fileId, outputFolder, timeLimit,
					numberOfSlaves, myRank, numberOfSearchSlaves);
		}

		// Stops the timer and stores the elapsed time
		timer.stop();
		boost::timer::cpu_times end_time = timer.elapsed();
		double timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
		// Saves elapsed time and best solution to output file
		string filename = outputFolder + fileId + "/" + executionId + "/result.txt";
		ofstream out(filename.c_str(), ios::out | ios::trunc); 
		if(!out) { 
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open output summary file.\n";
		    // TODO tratar excecao 
  		} 

		out << "Global time spent: " << timeSpent << endl;
		Imbalance imb = c->getImbalance();
		out << "I(P) = " << imb.getValue() << endl;
		stringstream ss;
		c->printClustering(ss);
		out << ss.str();
 		out.close();

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
	// initializes the logging subsystem
	this->initLogging(myRank);

	// codigo do processor lider
	if(myRank == 0) {
		cout << "Correlation clustering problem solver" << endl;
		// id used for output folders
		string executionId = TimeDateUtil::generateRandomId();

		try {
			string s_alpha;
			int numberOfIterations = 500, k = -1, l = 1;
			bool debug = false, profile = false, firstImprovementOnOneNeig = true;
			string inputFileDir, outputFolder;
			int timeLimit = 1800;
			int problemType = ClusteringProblem::CC_PROBLEM, functionType = GainFunction::IMBALANCE;
			int numberOfSlaves = np - 1;
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
				("k,k", po::value<int>(&k)->default_value(-1), "k parameter (RCC problem)")
				("time-limit", po::value<int>(&timeLimit)->default_value(1800), "maximum execution time")
				("input-file", po::value< vector<string> >(), "input file")
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

			if(not vm.count("output-folder")) {
				cout << "Please specify the output folder." << endl;
				return 1;
			}

			if (k != -1) {
				BOOST_LOG_TRIVIAL(debug) << "k value is " << k << ". RCC is enabled." << endl;
				problemType = ClusteringProblem::RCC_PROBLEM;
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
			BOOST_LOG_TRIVIAL(info) << "Number of GRASP processes in parallel is " << numberOfSlaves << endl;
			vector<fs::path> fileList;

			if (vm.count("input-file")) {
				fs::path filePath (vm["input-file"].as< vector<string> >().at(0));
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

			int numberOfSearchSlaves = 0;
			if(numberOfSlaves > 0) {
				numberOfSearchSlaves = calculateNumberOfSearchSlaves(np, numberOfSlaves);
			}
			if(np > 1) {
				mpi::communicator world;
				// broadcasts the numberOfSlaves to all processes
				// mpi::broadcast(world, numberOfSlaves, 0);
				for(int i = 1; i < np; i++) {
					world.send(i, ParallelGrasp::INPUT_MSG_NUM_SLAVES_TAG, numberOfSlaves);
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
							myRank, problemType, functionType, seed);
				} else {
					BOOST_LOG_TRIVIAL(info) << "Profile mode on." << endl;
					for(double alpha2 = 0.0F; alpha2 < 1.1F; alpha2 += 0.1F) {
						BOOST_LOG_TRIVIAL(info) << "Processing GRASP with alpha = " << std::setprecision(2) << alpha2 << endl;
						processInputFile(filePath, outputFolder, executionId, debug, alpha2, l, firstImprovementOnOneNeig,
								numberOfIterations, timeLimit, numberOfSlaves, numberOfSearchSlaves,
								myRank, problemType, functionType, seed);
					}
				}
			}
			// ------------------ M P I    T E R M I N A T I O N ---------------------
			if(np > 1) {
				mpi::communicator world;
				int i = 1;
				for(i = 1; i <= numberOfSlaves; i++) {
					InputMessageParallelGrasp imsg;
					world.send(i, ParallelGrasp::TERMINATE_MSG_TAG, imsg);
					BOOST_LOG_TRIVIAL(debug) << "Terminate message sent to process " << i << endl;
				}
				for(; i < np; i++) {
					InputMessageParallelVNS imsg;
					world.send(i, ParallelGrasp::TERMINATE_MSG_TAG, imsg);
					BOOST_LOG_TRIVIAL(debug) << "Terminate message sent to process " << i << endl;
				}
			}
			BOOST_LOG_TRIVIAL(info) << "Done.";
		}
		catch(std::exception& e)
		{
			BOOST_LOG_TRIVIAL(fatal) << "Abnormal program termination. Stracktrace: " << endl;
			BOOST_LOG_TRIVIAL(fatal) << e.what() << "\n";
			if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
				BOOST_LOG_TRIVIAL(fatal) << stack << endl;
			}
			BOOST_LOG_TRIVIAL(fatal) << diagnostic_information(e);
			return 1;
		}
	} else { // slave processes
		mpi::communicator world;
		// First message received from leader contains the number of GRASP slaves
		// and is used for a process to discover if it is a GRASP slave or a VNS slave
		// Grasp slave: 1 < rank < numberOfSlaves; VNS slave: otherwise
		unsigned int numberOfSlaves = 0;
		world.recv(ParallelGrasp::LEADER_ID, mpi::any_tag, numberOfSlaves);
		BOOST_LOG_TRIVIAL(trace) << "number of slaves received\n";
		int numberOfSearchSlaves = 0;
		if(numberOfSlaves > 0) {
			numberOfSearchSlaves = calculateNumberOfSearchSlaves(np, numberOfSlaves);
		}
		// Controls the number of messages received
		unsigned int messageCount = 0;
		// common slave variables
		ClusteringProblemFactory problemFactory;
		ParallelNeighborhoodSearch pnSearch(numberOfSlaves, numberOfSearchSlaves);
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g;
		unsigned int previousId = 0;

		while(true) {
			try {
				if(myRank <= numberOfSlaves) {  // GRASP slave
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " ready [GRASP slave process]." << endl;
					// Receives a message with GRASP parameters and triggers local GRASP execution
					InputMessageParallelGrasp imsgpg;
					// receives a message of type ParallelGrasp::INPUT_MSG_PARALLEL_GRASP_TAG or a terminate msg
					mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, imsgpg);
					if(stat.tag() == ParallelGrasp::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": terminate msg received.\n";
						return 0;
					}
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel GRASP]: Received message from leader." << endl;
					messageCount++;

					if(imsgpg.l == 0) {
						BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty GRASP message received. Terminating.\n";
						return 1;
					}
					// builds a new graph object only if it has changed (saves time)
					if(previousId != imsgpg.id) {
						// reconstructs the graph from its text representation
						g.reset();
						g = reader.readGraphFromString(imsgpg.graphInputFileContents);
						previousId = imsgpg.id;
					}

					// trigggers the local GRASP routine
					GainFunctionFactory functionFactory(g.get());
					Grasp resolution(&functionFactory.build(imsgpg.gainFunctionType), seed);
					ClusteringPtr bestClustering = resolution.executeGRASP(g.get(), imsgpg.iter, imsgpg.alpha,
							imsgpg.l, imsgpg.firstImprovementOnOneNeig, problemFactory.build(imsgpg.problemType),
							imsgpg.executionId, imsgpg.fileId, imsgpg.outputFolder, imsgpg.timeLimit,
							imsgpg.numberOfSlaves, myRank, imsgpg.numberOfSearchSlaves);

					// Sends the result back to the leader process
					OutputMessage omsg(*bestClustering);
					world.send(ParallelGrasp::LEADER_ID, ParallelGrasp::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": GRASP Output Message sent to leader." << endl;

				} else {  // VNS slave
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " ready [VNS slave process]." << endl;
					// Receives a parallel VNS message with parameters and triggers local VNS execution
					InputMessageParallelVNS imsgvns;
					// receives a message of type ParallelGrasp::INPUT_MSG_PARALLEL_VNS_TAG or a terminate msg
					mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, imsgvns);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel VNS]: Received message from leader." << endl;
					messageCount++;
					if(stat.tag() == ParallelGrasp::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": terminate msg received.\n";
						return 0;
					}

					if(imsgvns.l == 0) {
						BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty VNS message received. Terminating.\n";
						return 1;
					}
					// builds a new graph object only if it has changed (saves time)
					if(previousId != imsgvns.id) {
						// reconstructs the graph from its text representation
						g.reset();
						g = reader.readGraphFromString(imsgvns.graphInputFileContents);
						previousId = imsgvns.id;
					}

					// triggers the local partial VNS search to be done between initial and final cluster indices
					ClusteringPtr bestClustering = pnSearch.searchNeighborhood(imsgvns.l, g.get(), &imsgvns.clustering,
							problemFactory.build(imsgvns.problemType), imsgvns.timeSpentSoFar, imsgvns.timeLimit, seed, myRank,
							imsgvns.initialClusterIndex, imsgvns.finalClusterIndex, false);

					// Sends the result back to the leader process
					int leader = stat.source();
					OutputMessage omsg(*bestClustering);
					world.send(leader, ParallelGrasp::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
					BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": VNS Output Message sent to grasp leader number "
							<< leader << ". I(P) = " << bestClustering->getImbalance().getValue();
				}
			}
			catch(std::exception& e)
			{
				BOOST_LOG_TRIVIAL(fatal) << "Abnormal program termination. Stracktrace: " << endl;
				BOOST_LOG_TRIVIAL(fatal) << e.what() << "\n";
				if ( std::string const *stack = boost::get_error_info<stack_info>(e) ) {
					BOOST_LOG_TRIVIAL(fatal) << stack << endl;
				}
				BOOST_LOG_TRIVIAL(fatal) << diagnostic_information(e);
				return 1;
			}
		}
	}
	return 0;
}

unsigned int CommandLineInterfaceController::calculateNumberOfSearchSlaves(const unsigned int& np, const unsigned int& numberOfSlaves) {
	// Number of remaining slaves after using 1 (leader) + 'numberOfSlaves' processes for parallel GRASP execution
	unsigned int remainingSlaves = np - numberOfSlaves - 1;
	// Divides the remaining slaves that will be used in parallel VNS processing
	unsigned int numberOfSearchSlaves = remainingSlaves / (numberOfSlaves + 1);
	BOOST_LOG_TRIVIAL(info) << "The number of VNS search slaves per master is " << numberOfSearchSlaves << endl;

	return numberOfSearchSlaves;
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

void CommandLineInterfaceController::initLogging(int myRank) {
	using namespace boost::log;
	namespace keywords = boost::log::keywords;
	namespace sinks = boost::log::sinks;
	namespace expr = boost::log::expressions;
	namespace attrs = boost::log::attributes;
	string filename = string("logs/Node") + lexical_cast<string>(myRank) + string(".log");

	LogSeverityEnumParser parser;
	logging::trivial::severity_level severity = parser.ParseSomeEnum(logSeverity);

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
					<< "[" << trivial::severity << "]\t"
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
        logging::trivial::severity >= severity
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
