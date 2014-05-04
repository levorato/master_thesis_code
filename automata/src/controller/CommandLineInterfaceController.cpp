/*
 * CommandLineInterfaceController.cpp
 *
 *  Created on: May 3, 2014
 *      Author: Mario Costa Levorato Junior
 */

#include "include/CommandLineInterfaceController.h"
#include "../../../grasp/src/controller/include/SimpleTextGraphFileReader.h"
#include "../../../grasp/src/graph/include/Graph.h"
#include "../../../grasp/src/util/include/TimeDateUtil.h"
#include "../../../grasp/src/util/include/EnumUtil.h"

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
using namespace util;

namespace controller {

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

CommandLineInterfaceController::CommandLineInterfaceController() : logSeverity("warning") {
	// TODO Auto-generated constructor stub

}

CommandLineInterfaceController::~CommandLineInterfaceController() {
	// TODO Auto-generated destructor stub
}

void CommandLineInterfaceController::processInputFile(fs::path filePath, string& outputFolder,
		string& executionId, const bool& debug, const int& numberOfIterations, const long& timeLimit,
		const unsigned long& seed) {
	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g = reader.readGraphFromFile(filePath.string());
		string fileId = filePath.filename().string();
		boost::timer::cpu_timer timer;
		boost::timer::cpu_times start_time, end_time;
		double timeSpent = 0.0;
		string filename;
		// Execution additional info
		// ExecutionInfo info(executionId, fileId, outputFolder, myRank);
		// medicao de tempo do algoritmo
		timer.start();
		start_time = timer.elapsed();

		// Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();
		timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
		// Saves elapsed time and best solution to output file
		filename = outputFolder + "/" + fileId + "/" + executionId + "/result.txt";
		ofstream out(filename.c_str(), ios::out | ios::trunc);
		if(!out) {
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << filename;
			// TODO tratar excecao
		}
		BOOST_LOG_TRIVIAL(info) << "Global time spent: " << timeSpent << " s";
		out << "Global time spent: " << timeSpent << endl;
		// Closes the file
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

int CommandLineInterfaceController::processArgumentsAndExecute(int argc, char *argv[]) {

	// used for debugging purpose
	std::set_terminate( handler );

	cout << "Network Automata Simulation Algorithm" << endl;

	// random seed used in algorithm
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

	int numberOfIterations = 500;
	bool debug = false;
	string inputFileDir, outputFolder;
	int timeLimit = 1800;
	string jobid;

	po::options_description desc("Available options:");
	desc.add_options()
		("help", "show program options")
		("iterations,iter", po::value<int>(&numberOfIterations)->default_value(400),
			  "number of iterations of multi-start heuristic")
		("time-limit", po::value<int>(&timeLimit)->default_value(1800), "maximum execution time (seconds)")
		("input-file", po::value< std::vector<string> >(), "graph input file")
		("debug", po::value<bool>(&debug)->default_value(false), "enable debug mode")
		("input-file-dir", po::value<string>(&inputFileDir), "input graph file directory (processes all files inside)")
		("output-folder", po::value<string>(&outputFolder), "output folder for result files")
		("jobid", po::value<string>(&jobid), "job identifier (string)")
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
        this->initLogging(jobid);
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

	cout << "Network Automata Simulation Algorithm" << endl;
	// id used for output folders
	string executionId = jobid;

	try {
		if (vm.count("help")) {
			cout << "Usage: automata [options]\n";
			cout << desc;
			return 0;
		}

		if(not vm.count("output-folder")) {
			cout << "Please specify the output folder." << endl;
			return 1;
		}

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

		// -------------------  G R A P H     F I L E S     P R O C E S S I N G -------------------------
		for(unsigned int i = 0; i < fileList.size(); i++) {
			fs::path filePath = fileList.at(i);

			BOOST_LOG_TRIVIAL(info) << "Number of iterations is " << numberOfIterations << "\n";
			processInputFile(filePath, outputFolder, executionId, debug, numberOfIterations, timeLimit,
					seed);
		}

		// ------------------  T E R M I N A T I O N ---------------------
		BOOST_LOG_TRIVIAL(info) << "Terminating process...";
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

void CommandLineInterfaceController::initLogging(string executionId) {
	using namespace boost::log;
	namespace keywords = boost::log::keywords;
	namespace sinks = boost::log::sinks;
	namespace expr = boost::log::expressions;
	namespace attrs = boost::log::attributes;
	string filename = string("logs/") + executionId + string(".log");
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
