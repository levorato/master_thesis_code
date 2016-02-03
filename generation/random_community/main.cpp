#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
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
#include <boost/make_shared.hpp>
#include <boost/exception/all.hpp>
#include <exception>
#include <boost/exception/info.hpp>

using namespace boost::program_options;

#include <iostream>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;
using namespace std;

#include "RandomSGCommunity.h"
#include "TimeDateUtil.h"
using namespace generation;

void initLogging(int myRank) {
	using namespace boost::log;
	namespace logging = boost::log;
	namespace keywords = boost::log::keywords;
	namespace sinks = boost::log::sinks;
	namespace expr = boost::log::expressions;
	namespace attrs = boost::log::attributes;
	string executionId = util::TimeDateUtil::getDateAndTime();
	string filename = string("logs/") + executionId + string("-Node") + boost::lexical_cast<string>(myRank) + string(".log");
	cout << "Init logging..." << endl;

	logging::trivial::severity_level mySeverity = logging::trivial::info;

	boost::shared_ptr<boost::log::core> logger = boost::log::core::get();
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

int main(int argc, char* argv[])
{
    try {
    	// Inicializacao do MPI
		mpi::environment env(argc, argv);
		mpi::communicator world;

        long c = 0, n = 0, k = 0;
		double p_in = -1.0, p_minus = -1.0, p_plus = -1.0;
        
        options_description desc("Generate random unweighted signed graphs with a predefined community structure.");
        desc.add_options()
        // First parameter describes option name/short name
        // The second is parameter to option
        // The third is description
        ("help,h", "print usage message")
        ("c", value<long>(&c), "<number of clusters>")
		("n", value<long>(&n), "<number of vertices in each cluster>")
		("k", value<long>(&k), "<degree of each vertex>")
		("pin", value<double>(&p_in), "<probability of each node connecting other nodes in the same community>")
        ("p_minus", value<double>(&p_minus), "<probability of negative links appearing within communities>")
		("p_plus", value<double>(&p_plus), "<probability of of positive links appearing between communities>")
        ;
		
		variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        try {
			initLogging(world.rank());
		}
		catch(std::exception& e)
		{
			 cerr << "Fatal application error in log init.\n";
			 cerr << "Abnormal program termination. Stracktrace: " << endl;
			 cerr << e.what() << "\n";
			 return 1;
		}

        if (vm.count("help")) {  
            cout << desc << "\n";
            return 0;
        }
        if(not ((p_in >= 0.0 and p_in <= 1.0) or (p_minus >= 0.0 and p_minus <= 1.0) or (p_plus >= 0.0 and p_plus <= 1.0))) {
        	BOOST_LOG_TRIVIAL(fatal) << "Please specify the parameters pin, p_minus and p_plus in the range [0, 1]" << endl;
			cout << "Please specify the parameters pin, p_minus and p_plus in the range [0, 1]." << endl;
			return 1;
		}

		cout << "Generating random signed graph with community structure with parameters c = " << c << ", n = " 
			<< n << ", k = " << k << ", pin = " << pin << ", p_minus = " << p_minus << ", p_plus = " << p_plus << endl;
		BOOST_LOG_TRIVIAL(info) << "Generating random signed graph with community structure with parameters c = " << c << ", n = " 
			<< n << ", k = " << k << ", pin = " << pin << ", p_minus = " << p_minus << ", p_plus = " << p_plus;
        
		RandomSGCommunity generator;
		generator.generateRandomSG(c, n, k, pin, p_minus, p_plus, world.rank(), world.size());

		// ------------------ M P I    T E R M I N A T I O N ---------------------
		if(world.size() > 1) {
			if(world.rank() == 0) {
				BOOST_LOG_TRIVIAL(info) << "Terminating MPI worker processes...";
				InputMessage imsg;
				for(int i = 1; i < world.size(); i++) {
					world.send(i, InputMessage::TERMINATE_MSG_TAG, imsg);
					BOOST_LOG_TRIVIAL(info) << "Terminate message sent to process " << i;
				}
			} else {
				// wait for terminate message
				InputMessage imsg;
				mpi::status msg = world.recv(0, InputMessage::TERMINATE_MSG_TAG, imsg);
			}
		}
    }
    catch(exception& e) {
        cerr << e.what() << "\n";
        BOOST_LOG_TRIVIAL(fatal) << e.what() << "\n";
    }
}



