/*
 * CommandLineInterfaceController.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: mario
 */

#include "include/CommandLineInterfaceController.h"
#include "graph/include/SimpleTextGraphFileReader.h"
#include "graph/include/Graph.h"
#include "../resolution/grasp/include/Grasp.h"
#include "../resolution/grasp/include/ParallelGrasp.h"
#include "../resolution/ils/include/ILS.h"
#include "../resolution/ils/include/ParallelILS.h"
#include "graph/include/Clustering.h"
#include "graph/include/Imbalance.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "util/include/TimeDateUtil.h"
#include "util/include/EnumUtil.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "problem/include/ClusteringProblemFactory.h"
#include "../resolution/construction/include/GainFunctionFactory.h"
#include "../resolution/construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "../resolution/construction/include/ConstructClustering.h"
#include "util/include/RandomUtil.h"
#include "../resolution/vnd/include/CUDAVariableNeighborhoodDescent.h"
#include "../resolution/vnd/include/CUDANeighborhoodSearch.h"
#include "../resolution/grasp/include/CUDAGrasp.h"
#include "../resolution/construction/include/CUDAConstructClustering.h"
#include "../resolution/construction/include/CUDAImbalanceGainFunction.h"
#include "../resolution/ils/include/CUDAILS.h"
#include "../test/controller/include/TestController.h"

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
using namespace resolution::construction;
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
    else if (token == "ILS")
    	strategy = CommandLineInterfaceController::ILS;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

std::istream& operator>>(std::istream& in, CommandLineInterfaceController::SearchName& search)
{
    std::string token;
    in >> token;
    if (token == "SEQUENTIAL")
        search = CommandLineInterfaceController::SEQUENTIAL_SEARCH;
    else if (token == "PARALLEL")
        search = CommandLineInterfaceController::PARALLEL_SEARCH;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

CommandLineInterfaceController::CommandLineInterfaceController() : logSeverity("info") {
	// TODO Auto-generated constructor stub

}

CommandLineInterfaceController::~CommandLineInterfaceController() {
	// TODO Auto-generated destructor stub
}

void CommandLineInterfaceController::processInputFile(fs::path filePath, string& outputFolder,
		string& executionId, const bool& debug, const double& alpha, const int& l, const bool& firstImprovementOnOneNeig,
		const int& numberOfIterations, const long& timeLimit, const int& machineProcessAllocationStrategy,
		const int& numberOfMasters, const int& numberOfSearchSlaves,
		const int& myRank, const int& functionType, const unsigned long& seed, const bool& CCEnabled, const bool& RCCEnabled,
		long k, const StategyName& resolutionStrategy, const SearchName& searchType, const int& iterMaxILS, const int& perturbationLevelMax,
		const bool& splitGraph, const bool& cuda) {
	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g = reader.readGraphFromFile(filePath.string());
		Clustering c;
		string fileId = filePath.filename().string();
		ClusteringProblemFactory problemFactory;
		GainFunctionFactory functionFactory(g.get());
		boost::timer::cpu_timer timer;
		boost::timer::cpu_times start_time, end_time;
		double timeSpent = 0.0;
		string filename;
		// Construct clustering module
		ConstructClustering construct(functionFactory.build(functionType), seed, alpha);
		// Chooses between the sequential or parallel search algorithm
		NeighborhoodSearch* neigborhoodSearch;
		NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves);
		if(searchType == CommandLineInterfaceController::PARALLEL_SEARCH) {
			neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::PARALLEL);
		} else {
			neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
		}
		// VND - local search module
		CUDANeighborhoodSearch neighborhoodSearchCUDA;
		SequentialNeighborhoodSearch neighborhoodSearchSeq;
		VariableNeighborhoodDescent vnd(neighborhoodSearchSeq, seed, l, firstImprovementOnOneNeig, timeLimit);
		CUDAVariableNeighborhoodDescent cudavnd(neighborhoodSearchSeq, seed, l, firstImprovementOnOneNeig, timeLimit);
		// Execution additional info
		ExecutionInfo info(executionId, fileId, outputFolder, myRank);

		// -------------------  C C    P R O C E S S I N G -------------------------
		if(CCEnabled) {
			// medicao de tempo do CC
			timer.start();
			start_time = timer.elapsed();

			if(resolutionStrategy == GRASP) {
				//   G R A S P
				if(numberOfMasters == 0) {	// sequential version of GRASP
					Grasp resolution;
					CUDAGrasp CUDAgrasp;
					if(cuda) {
						c = CUDAgrasp.executeGRASP(&construct, &cudavnd, g.get(), numberOfIterations,
												problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
					} else {
						c = resolution.executeGRASP(&construct, &vnd, g.get(), numberOfIterations,
												problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
					}
				} else {  // parallel version
					// distributes GRASP processing among numberOfSlaves processes and summarizes the result
					ParallelGrasp parallelResolution(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves, splitGraph, cuda);
					c = parallelResolution.executeGRASP(&construct, &vnd, g.get(), numberOfIterations,
												problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
				}
			} else if(resolutionStrategy == ILS) {
				//   I L S
				if(numberOfMasters == 0) {	// sequential version of ILS
					resolution::ils::ILS resolution;
					resolution::ils::CUDAILS CUDAils;
					if(cuda) {
						c = CUDAils.executeILS(&construct, &cudavnd, g.get(), numberOfIterations, iterMaxILS,
												perturbationLevelMax, problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
					} else {
						c = resolution.executeILS(&construct, &vnd, g.get(), numberOfIterations, iterMaxILS,
												perturbationLevelMax, problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
					}
				} else {  // parallel version
					// distributes ILS processing among numberOfMasters processes and summarizes the result
					resolution::ils::ParallelILS parallelResolution(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves, splitGraph, cuda);
					c = parallelResolution.executeILS(&construct, &vnd, g.get(), numberOfIterations, iterMaxILS,
												perturbationLevelMax, problemFactory.build(ClusteringProblem::CC_PROBLEM), info);
				}
			}
			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			// Saves elapsed time and best solution to output file
			filename = outputFolder + "/" + fileId + "/" + executionId + "/cc-result.txt";
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
			c.printClustering(ss, g->getN());
			out << ss.str();
			// Outputs additional graph analysis data
			CCProblem& ccp = static_cast<CCProblem&>(problemFactory.build(ClusteringProblem::CC_PROBLEM));
			string analysis = ccp.analyzeImbalance(*g, c);
			out << analysis << endl;
			// Closes the file
			out.close();
		}

 		// -------------------  R C C     P R O C E S S I N G -------------------------
 		if(RCCEnabled) {
 			ConstructClustering* RCCConstruct = &construct;
 			ConstructClustering NoConstruct(functionFactory.build(functionType), seed, alpha, &c);
 			// 	If k = 0 (not specified), uses CC best result (cluster c) number of clusters as input (k)
 			if(k == 0) {
 				k = c.getNumberOfClusters();
 				BOOST_LOG_TRIVIAL(info) << "RCC Problem: Using CC solution number of clusters as RCC k value: " << k << ".";
 			} else if(k < 0) {  // if k < 0, skips the construct clustering phase in the metaheuristic and uses the best solution from CC alg.
                RCCConstruct = &NoConstruct;
 			}

 			// medicao de tempo do RCC
			timer.start();
			start_time = timer.elapsed();
			Clustering RCCCluster;
			Imbalance CCimb = c.getImbalance();

			// if the CC result is already zero, no need to execute RCC at all
			// RCC problem solving does not have a CUDA-enabled version yet
			if(CCimb.getValue() > 0.0) {
				if(resolutionStrategy == GRASP) {
					//   G R A S P
					if(numberOfMasters == 0) {	// sequential version of GRASP
						Grasp resolution;
						RCCCluster = resolution.executeGRASP(RCCConstruct, &vnd, g.get(), numberOfIterations,
								problemFactory.build(ClusteringProblem::RCC_PROBLEM, k), info);
					} else {  // parallel version
						// distributes GRASP processing among numberOfMasters processes and summarizes the result
						ParallelGrasp parallelResolution(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves, splitGraph, cuda);
						RCCCluster = parallelResolution.executeGRASP(&construct, &vnd, g.get(), numberOfIterations,
								problemFactory.build(ClusteringProblem::RCC_PROBLEM, k), info);
					}
				} else if(resolutionStrategy == ILS) {
					//   I L S
					if(numberOfMasters == 0) {	// sequential version of ILS
						resolution::ils::ILS resolution;
						RCCCluster = resolution.executeILS(RCCConstruct, &vnd, g.get(), numberOfIterations, iterMaxILS,
								perturbationLevelMax, problemFactory.build(ClusteringProblem::RCC_PROBLEM, k), info);
					} else {  // parallel version
						// distributes ILS processing among numberOfMasters processes and summarizes the result
						resolution::ils::ParallelILS parallelResolution(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves, splitGraph, cuda);
						RCCCluster = parallelResolution.executeILS(&construct, &vnd, g.get(), numberOfIterations, iterMaxILS,
								perturbationLevelMax, problemFactory.build(ClusteringProblem::RCC_PROBLEM, k), info);
					}
				}
			} else {  // RCC equals CC solution ( I(P) = RI(P) = 0 )
				RCCCluster = c;
			}
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
			// If the RCC solution value is WORSE than the CC Solution, discard the RCC solution and
			// present the CC solution
			Imbalance RCCimb = RCCCluster.getImbalance();
			if(CCimb.getValue() < RCCimb.getValue()) {
				RCCCluster = c;
				RCCimb = CCimb;
			}
			out << "SRI(P) = " << RCCimb.getValue() << endl;
			stringstream ss;
			RCCCluster.printClustering(ss, g->getN());
			out << ss.str();
			// Outputs additional graph analysis data
			RCCProblem& rccp = static_cast<RCCProblem&>(problemFactory.build(ClusteringProblem::RCC_PROBLEM));
			string analysis = rccp.analyzeImbalance(*g, RCCCluster);
			out << analysis << endl;
			// Closes the file
			out.close();
			BOOST_LOG_TRIVIAL(info) << "RCC Solve done. Obj = " << RCCimb.getValue();
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
	//cout << "Correlation clustering problem solver" << endl;

	// random seed used in the algorithms
	/*
	* Caveat: std::time(0) is not a very good truly-random seed.  When
	* called in rapid succession, it could return the same values, and
	* thus the same random number sequences could ensue.
	* Instead, we are using boost::random_device
	* http://stackoverflow.com/questions/4329284/c-boost-random-numeric-generation-problem
	* http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
	*/
	unsigned long seed = mix(clock(), time(NULL), getpid());
	RandomUtil randomUtil;
	randomUtil.setSeed(seed);

	// reads the system properties from ini file
	this->readPropertiesFile();

	string s_alpha;
	int numberOfIterations = 500, l = 1;
	long k = 0;
	bool debug = false, profile = false, firstImprovementOnOneNeig = true, CCEnabled = true, RCCEnabled = true;
	string inputFileDir, outputFolder;
	int timeLimit = 1800;
	int functionType = GainFunction::IMBALANCE;
	int numberOfMasters = np - 1;
	int totalNumberOfVNDSlaves = -1;
	string jobid;
	CommandLineInterfaceController::StategyName strategy = CommandLineInterfaceController::GRASP;
	CommandLineInterfaceController::SearchName searchType = CommandLineInterfaceController::SEQUENTIAL_SEARCH;
	int iterMaxILS = 5, perturbationLevelMax = 30;  // for Slashdot and completely Random instances
	// int iterMaxILS = 5, perturbationLevelMax = 3; // for UNGA instances
	bool splitGraph = false;
	bool cuda = true;
	bool test = false;
	bool exportDOT = false;

	po::options_description desc("Available options:");
	desc.add_options()
		("help", "show program options")
		("split", po::value<bool>(&splitGraph)->default_value(false), "Enable graph partitioning when solving in parallel (requires slaves > 0)")
		("cuda", po::value<bool>(&cuda)->default_value(true), "Enable CUDA local search (requires NVIDIA GPU)")
		("test", po::value<bool>(&test)->default_value(false), "Enable test suite execution (for development purpose only)")
		("alpha,a", po::value<string>(&s_alpha),
			  "alpha - randomness factor of constructive phase")
		("iterations,iter", po::value<int>(&numberOfIterations)->default_value(400),
			  "number of iterations of multi-start heuristic")
		("neighborhood_size,l", po::value<int>(&l)->default_value(1),
			  "neighborhood size of local search")
		("cc", po::value<bool>(&CCEnabled)->default_value(true), "Enable CC Problem resolution")
		("rcc", po::value<bool>(&RCCEnabled)->default_value(true), "Enable RCC Problem resolution")
		("k", po::value<long>(&k)->default_value(0), "RCC Problem k parameter (max number of clusters) - optional")
		("perturbationLevelMax", po::value<int>(&perturbationLevelMax)->default_value(30), "ILS max perturbation level")
		("time-limit", po::value<int>(&timeLimit)->default_value(1800), "maximum execution time (seconds)")
		("input-file", po::value< std::vector<string> >(), "graph input file")
		("debug", po::value<bool>(&debug)->default_value(false), "enable debug mode")
		("profile", po::value<bool>(&profile)->default_value(false), "enable profile mode (varies alpha values)")
		("input-file-dir", po::value<string>(&inputFileDir), "input graph file directory (processes all files inside)")
		("output-folder", po::value<string>(&outputFolder), "output folder for result files")
		("gain-function-type", po::value<int>(&functionType),
				"0 for min imbalance, 1 for max modularity gain function, 2 for max negative modularity gain function, "
				"3 for max positive-negative modularity gain function, 4 for max pos-neg mod gain function II, "
				"5 for max pos-neg mod gain function III  (Gain Function is used in construction phase)")
		("slaves,s", po::value<int>(&numberOfMasters)->default_value(np - 1), "number of masters (main multi-start heuristic slaves) in parallel")
		("vndslaves", po::value<int>(&totalNumberOfVNDSlaves)->default_value(-1), "number of VND slaves in parallel")
		("firstImprovementOnOneNeig", po::value<bool>(&firstImprovementOnOneNeig)->default_value(true), "first improvement in 1-opt neighborhood (** sequential VND only **)")
		("jobid", po::value<string>(&jobid), "job identifier (string)")
		("strategy", po::value<StategyName>(&strategy),
					 "Resolution strategy to be used. Accepted values: GRASP (default), ILS.")
		("search", po::value<CommandLineInterfaceController::SearchName>(&searchType),
                                         "Local search to be used. Accepted values: SEQUENTIAL (default), PARALLEL.")
		("exportDOT", po::value<bool>(&exportDOT)->default_value(false), "Export graph file to GraphViz DOT format.")
		//("ils,i", po::value<int>(&iterMaxILS)->default_value(5), "number of iterations of internal ILS loop")
		//("perturb,p", po::value<int>(&perturbationLevelMax)->default_value(30), "maximum perturbation level in ILS")
	;
	po::positional_options_description p;
	p.add("input-file", -1);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), vm);
	po::notify(vm);

	// job id is obtained through command line parameter from PBS Scheduler
	//cout << "Job id is " << jobid << "\n";
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

	unsigned int numberOfSearchSlavesPerMaster = 0;
	int machineProcessAllocationStrategy = 0;
	// Leader process code (rank 0)
	if(myRank == 0) {
		//cout << "Correlation clustering problem solver" << endl;
		// id used for output folders
		string executionId = jobid;

		if(exportDOT) {
			cout << "Exporting the graph to DOT format." << endl;
			SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
			fs::path filePath (vm["input-file"].as< std::vector<string> >().at(0));
			SignedGraphPtr g = reader.readGraphFromFile(filePath.string());
			string loc = g->getGraphFileLocation();
			string filename(loc.substr(loc.find_last_of("/\\")));
			reader.exportGraphToGraphVizFile(*g, outputFolder, filename);

			return 0;
		}

		try {
			if (vm.count("help")) {
				cout << "Usage: graspcc [options]\n";
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
			if (splitGraph) {
				if(numberOfMasters > 0) {
					BOOST_LOG_TRIVIAL(info) << "Graph split is enabled. Partitioning graph between master and slaves.";
				} else {
					BOOST_LOG_TRIVIAL(error) << "WARNING: Graph split is disabled, since slaves = 0. Please set the slaves parameter accordingly.";
					splitGraph = false;
				}
			} else {
				BOOST_LOG_TRIVIAL(info) << "Graph split is disabled.";
			}

			if(cuda) {
				BOOST_LOG_TRIVIAL(info) << "CUDA local search is enabled.";
			} else {
				BOOST_LOG_TRIVIAL(info) << "CUDA local search is disabled.";
			}

			if(strategy == GRASP) {
				BOOST_LOG_TRIVIAL(info) << "Resolution strategy is GRASP.";
			} else if(strategy == ILS) {
				BOOST_LOG_TRIVIAL(info) << "Resolution strategy is ILS.";
				BOOST_LOG_TRIVIAL(info) << "Maximum perturbation level is " << perturbationLevelMax;
				BOOST_LOG_TRIVIAL(info) << "The number of ILS iterations is " << iterMaxILS;
			}
			if(searchType == CommandLineInterfaceController::SEQUENTIAL_SEARCH) {
					BOOST_LOG_TRIVIAL(info) << "Local search type is SEQUENTIAL.";
			} else if(searchType == CommandLineInterfaceController::PARALLEL_SEARCH) {
					BOOST_LOG_TRIVIAL(info) << "Local search type is PARALLEL.";
			}
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
			if (CCEnabled) {
				BOOST_LOG_TRIVIAL(info) << "CC is enabled.";
			} else {
				BOOST_LOG_TRIVIAL(info) << "CC is disabled.";
			}
			if (RCCEnabled) {
				BOOST_LOG_TRIVIAL(info) << "RCC is enabled.";
				// if k parameter was not informed, CC problem should be enabled
				// if CC is not enabled, issue an error
				if((k <= 0) and (not CCEnabled)) {
					BOOST_LOG_TRIVIAL(fatal) << "Please specify RCC k value or enable CC Problem resolution.";
					return 1;
				} else if(k != 0) {
					BOOST_LOG_TRIVIAL(info) << "RCC k value is " << k << ".";
				}
			} else {
				BOOST_LOG_TRIVIAL(info) << "RCC is disabled.";
			}
			if(alpha <= 0.0) {
				BOOST_LOG_TRIVIAL(info) << "VOTE is enabled. Will always choose best-gain vertex on constructionPhase." << endl;
			}
			if(totalNumberOfVNDSlaves < 0) {  // if numberOfVNDSlaves is not specified, defaults to numproc - (masters - 1)
				totalNumberOfVNDSlaves = np - numberOfMasters - 1;
			}

			BOOST_LOG_TRIVIAL(info) << "Total number of processes is " << np << endl;
            //cout << "Total number of processes is " << np << endl;
			BOOST_LOG_TRIVIAL(info) << "Number of master processes in parallel is " << (numberOfMasters + 1);
			//cout << "Number of master processes in parallel is " << (numberOfMasters + 1) << endl;
			BOOST_LOG_TRIVIAL(info) << "Number of VND processes in parallel is " << totalNumberOfVNDSlaves << endl;
            //cout << "Number of VND processes in parallel is " << totalNumberOfVNDSlaves << endl;
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

			MPIInitParams mpiparams(numberOfMasters, numberOfSearchSlavesPerMaster, MPIUtil::ALL_MASTERS_FIRST, searchType);
			if(totalNumberOfVNDSlaves > 0) {
				numberOfSearchSlavesPerMaster = MPIUtil::calculateNumberOfSearchSlaves(np, numberOfMasters);
				mpiparams.numberOfSearchSlavesPerMaster = numberOfSearchSlavesPerMaster;
				mpiparams.machineProcessAllocationStrategy = MPIUtil::MASTER_AND_VND_SLAVES_TOGETHER;
			} else if(np / (numberOfMasters+1) == 4) {
				mpiparams.machineProcessAllocationStrategy = MPIUtil::MASTER_AND_VND_SLAVES_TOGETHER;
			}
			machineProcessAllocationStrategy = mpiparams.machineProcessAllocationStrategy;
			if(np > 1) {
				// broadcasts the mpi params to all processes
				// mpi::broadcast(world, numberOfSlaves, 0);
				for(int i = 1; i < np; i++) {
					world.send(i, MPIMessage::INPUT_MSG_MPI_PARAMS_TAG, mpiparams);
				}
			}

			if(test) {
				BOOST_LOG_TRIVIAL(info) << "TEST SUITE EXECUTION IS ENABLED.";
				TestController testCtl;
				if(testCtl.runTestSuite()) {
					BOOST_LOG_TRIVIAL(info) << "TEST SUITE EXECUTION SUCCESSFUL.";
				} else {
					BOOST_LOG_TRIVIAL(info) << "TEST SUITE EXECUTION FAILED!";
				}
			} else {
				// -------------------  G R A P H     F I L E S     P R O C E S S I N G -------------------------
				for(unsigned int i = 0; i < fileList.size(); i++) {
					fs::path filePath = fileList.at(i);

					if(not profile) { // default mode: use specified parameters
						BOOST_LOG_TRIVIAL(info) << "Alpha value is " << std::setprecision(2) << fixed << alpha << "\n";
						BOOST_LOG_TRIVIAL(info) << "Number of iterations is " << numberOfIterations << "\n";
						processInputFile(filePath, outputFolder, executionId, debug, alpha, l, firstImprovementOnOneNeig,
								numberOfIterations, timeLimit, mpiparams.machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlavesPerMaster,
								myRank, functionType, seed, CCEnabled, RCCEnabled, k, strategy, searchType, iterMaxILS, perturbationLevelMax, splitGraph, cuda);
					} else {
						BOOST_LOG_TRIVIAL(info) << "Profile mode on." << endl;
						for(double alpha2 = 0.0F; alpha2 < 1.1F; alpha2 += 0.1F) {
							BOOST_LOG_TRIVIAL(info) << "Processing problem with alpha = " << std::setprecision(2) << alpha2 << endl;
							processInputFile(filePath, outputFolder, executionId, debug, alpha2, l, firstImprovementOnOneNeig,
									numberOfIterations, timeLimit, mpiparams.machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlavesPerMaster,
									myRank, functionType, seed, CCEnabled, RCCEnabled, k, strategy, searchType, iterMaxILS, perturbationLevelMax, splitGraph, cuda);
						}
					}
				}
			}

			// ------------------ M P I    T E R M I N A T I O N ---------------------
			terminateMPIProcessesIfAny(np, mpiparams.machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlavesPerMaster);
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
			terminateMPIProcessesIfAny(np, machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlavesPerMaster);
			BOOST_LOG_TRIVIAL(info) << "Done.";
			return 1;
		}
	} else { // slave processes
		// First message received from leader contains the number of masters
		// and is used for a process to discover if it is a master or a VND slave
		MPIInitParams mpiparams;
		world.recv(MPIMessage::LEADER_ID, MPIMessage::INPUT_MSG_MPI_PARAMS_TAG, mpiparams);
		BOOST_LOG_TRIVIAL(info) << "number of masters received: " << mpiparams.numberOfMasters;
		BOOST_LOG_TRIVIAL(info) << "number of VND slaves per master received: " << mpiparams.numberOfSearchSlavesPerMaster;
		BOOST_LOG_TRIVIAL(info) << "Machine-process allocation strategy received: " << mpiparams.machineProcessAllocationStrategy;
		BOOST_LOG_TRIVIAL(info) << "Local search type received: " << mpiparams.searchType;
		//cout << "number of masters received: " << mpiparams.numberOfMasters << "\n";

		// Controls the number of messages received
		unsigned int messageCount = 0;
		// common slave variables
		ClusteringProblemFactory problemFactory;
		SimpleTextGraphFileReader reader = SimpleTextGraphFileReader();
		SignedGraphPtr g;
		unsigned int previousId = 0;

		try {
			if(MPIUtil::isMaster(mpiparams.machineProcessAllocationStrategy, myRank, mpiparams.numberOfMasters, mpiparams.numberOfSearchSlavesPerMaster)) {  // master process
				BOOST_LOG_TRIVIAL(info) << "Process " << myRank << " ready [master process].";
				//cout << "Process " << myRank << " ready [master process]." << endl;
				while(true) {
					mpi::status stat = world.probe(MPIMessage::LEADER_ID, mpi::any_tag);
					if(stat.tag() == MPIMessage::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": terminate message received.";
						return 0;
					} else if(stat.tag() == MPIMessage::INPUT_MSG_PARALLEL_GRASP_TAG) {
						// Receives a message with GRASP parameters and triggers local GRASP execution
						InputMessageParallelGrasp imsgpg;
						mpi::status msg = world.recv(MPIMessage::LEADER_ID, mpi::any_tag, imsgpg);
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

						// triggers the local GRASP routine
						// Chooses between the sequential or parallel search algorithm
						NeighborhoodSearch* neigborhoodSearch;
						NeighborhoodSearchFactory nsFactory(mpiparams.machineProcessAllocationStrategy, imsgpg.numberOfMasters, imsgpg.numberOfSearchSlaves);
						if(mpiparams.searchType == MPIUtil::PARALLEL_SEARCH) {
							neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::PARALLEL);
						} else {
							neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
						}
						GainFunctionFactory functionFactory(g.get());
						ConstructClustering defaultConstruct(functionFactory.build(imsgpg.gainFunctionType), seed, imsgpg.alpha);
						ConstructClustering noConstruct(functionFactory.build(imsgpg.gainFunctionType), seed, imsgpg.alpha, &imsgpg.CCclustering);
						ConstructClustering* construct = &defaultConstruct;
						if((imsgpg.problemType == ClusteringProblem::RCC_PROBLEM) and (imsgpg.k < 0)) {
								construct = &noConstruct;
						}
						VariableNeighborhoodDescent vnd(*neigborhoodSearch, seed, imsgpg.l, imsgpg.firstImprovementOnOneNeig,
								imsgpg.timeLimit);
						// Execution additional info
						ExecutionInfo info(imsgpg.executionId, imsgpg.fileId, imsgpg.outputFolder, myRank);
						Grasp resolution;
						CUDAGrasp CUDAgrasp;
						Clustering bestClustering;
						// global vertex id array used in split graph
						std::vector<long> globalVertexId;
						// time spent on processing
						double timeSpent = 0.0;

						// If split graph is enabled (= vertexList not empty), creates a subgraph induced by the vertex set
						if(imsgpg.vertexList.size() > 0) {
							BOOST_LOG_TRIVIAL(info) << "Split graph partial graph";
							// SubGraph subg = (g->graph).create_subgraph();
							SignedGraph sg(g->graph, imsgpg.vertexList);
							/*
							for(std::vector<long>::iterator it = imsgpg.vertexList.begin(); it != imsgpg.vertexList.end(); it++) {
								add_vertex(*it, subg);
							}
							sg.graph = subg; */
							BOOST_LOG_TRIVIAL(info) << "Processing subgraph with n =  " << num_vertices(sg.graph) << ", " << "e =  " << num_edges(sg.graph);

							GainFunctionFactory functionFactory(&sg);
							ConstructClustering defaultConstruct(functionFactory.build(imsgpg.gainFunctionType), seed, imsgpg.alpha);
							ConstructClustering noConstruct(functionFactory.build(imsgpg.gainFunctionType), seed, imsgpg.alpha, &imsgpg.CCclustering);
							ConstructClustering* construct = &defaultConstruct;
							if((imsgpg.problemType == ClusteringProblem::RCC_PROBLEM) and (imsgpg.k < 0)) {
									construct = &noConstruct;
							}

							if(imsgpg.cudaEnabled) {
								bestClustering = CUDAgrasp.executeGRASP(construct, &vnd, &sg, imsgpg.iter,
													problemFactory.build(imsgpg.problemType, imsgpg.k), info);
								timeSpent = CUDAgrasp.getTotalTimeSpent();
							} else {
								bestClustering = resolution.executeGRASP(construct, &vnd, &sg, imsgpg.iter,
													problemFactory.build(imsgpg.problemType, imsgpg.k), info);
								timeSpent = resolution.getTotalTimeSpent();
							}
						} else {
							if(imsgpg.cudaEnabled) {
								bestClustering = CUDAgrasp.executeGRASP(construct, &vnd, g.get(), imsgpg.iter,
													problemFactory.build(imsgpg.problemType, imsgpg.k), info);
								timeSpent = CUDAgrasp.getTotalTimeSpent();
							} else {
								bestClustering = resolution.executeGRASP(construct, &vnd, g.get(), imsgpg.iter,
													problemFactory.build(imsgpg.problemType, imsgpg.k), info);
								timeSpent = resolution.getTotalTimeSpent();
							}
						}

						// Sends the result back to the leader process
						OutputMessage omsg(bestClustering, CUDAgrasp.getNumberOfTestedCombinations(), timeSpent, globalVertexId);
						world.send(MPIMessage::LEADER_ID, MPIMessage::OUTPUT_MSG_PARALLEL_GRASP_TAG, omsg);
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": GRASP Output Message sent to leader.";

					} else if(stat.tag() == MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG) {
						// Receives a message with ILS parameters and triggers local ILS execution
						InputMessageParallelILS imsgpils;
						mpi::status msg = world.recv(MPIMessage::LEADER_ID, mpi::any_tag, imsgpils);
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel ILS]: Received message from leader.";
						messageCount++;

						if(imsgpils.l == 0) {
							BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty ILS message received. Terminating.\n";
							return 1;
						}
						// builds a new graph object only if it has changed (saves time)
						if(previousId != imsgpils.id) {
							g.reset();
							// Reads the graph from the specified text file
							g = reader.readGraphFromFile(imsgpils.graphInputFilePath);
							previousId = imsgpils.id;
						}

						// triggers the local ILS routine
						// Chooses between the sequential or parallel search algorithm
						NeighborhoodSearch* neigborhoodSearch;
						NeighborhoodSearchFactory nsFactory(mpiparams.machineProcessAllocationStrategy,
								imsgpils.numberOfMasters, imsgpils.numberOfSearchSlaves);
						if(mpiparams.searchType == MPIUtil::PARALLEL_SEARCH) { 
							neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::PARALLEL);
						} else {
							neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
						}
						GainFunctionFactory functionFactory(g.get());
						ConstructClustering defaultConstruct(functionFactory.build(imsgpils.gainFunctionType), seed, imsgpils.alpha);
						ConstructClustering noConstruct(functionFactory.build(imsgpils.gainFunctionType), seed, imsgpils.alpha, &imsgpils.CCclustering);
						ConstructClustering* construct = &defaultConstruct;
						if((imsgpils.problemType == ClusteringProblem::RCC_PROBLEM) and (imsgpils.k < 0)) {
								construct = &noConstruct;
						}
						VariableNeighborhoodDescent vnd(*neigborhoodSearch, seed, imsgpils.l, imsgpils.firstImprovementOnOneNeig,
								imsgpils.timeLimit);
						// Additional execution info
						ExecutionInfo info(imsgpils.executionId, imsgpils.fileId, imsgpils.outputFolder, myRank);
						resolution::ils::ILS resolution;
						resolution::ils::CUDAILS CUDAILS;
						Clustering bestClustering;
						// global vertex id array used in split graph
						std::vector<long> globalVertexId;
						// time spent on processing
						double timeSpent = 0.0;
						// info about the processed graph
						long n = num_vertices(g->graph), e = num_edges(g->graph);

						// If split graph is enabled, creates a subgraph induced by the vertex set
						if(imsgpils.isSplitGraph) {
							BOOST_LOG_TRIVIAL(info) << "Split graph partial graph";
							if(imsgpils.vertexList.size() == 0) {
								BOOST_LOG_TRIVIAL(info) << "Empty subgraph, returning zero imbalance.";
								std::vector<unsigned int> clusterProcessOrigin;
								ClusterArray cArray(g->getN(), Clustering::NO_CLUSTER);
								bestClustering = Clustering(cArray, *g, problemFactory.build(imsgpils.problemType, imsgpils.k),
										0.0, 0.0, clusterProcessOrigin);
								n = 0;
								e = 0;
							} else {
								SignedGraph sg(g->graph, imsgpils.vertexList);
								/*
								sg.graph = (g->graph).create_subgraph(); // imsgpils.vertexList.begin(), imsgpils.vertexList.end());
								for(std::vector<long>::iterator it = imsgpils.vertexList.begin(); it != imsgpils.vertexList.end(); it++) {
									add_vertex(*it, sg.graph);
								} */
								n = num_vertices(sg.graph);
								e = num_edges(sg.graph);
								BOOST_LOG_TRIVIAL(info) << "Processing subgraph with n =  " << n << ", " << "e =  " << e;

								GainFunctionFactory functionFactory(&sg);
								ConstructClustering defaultConstruct(functionFactory.build(imsgpils.gainFunctionType), seed, imsgpils.alpha);
								ConstructClustering noConstruct(functionFactory.build(imsgpils.gainFunctionType), seed, imsgpils.alpha, &imsgpils.CCclustering);
								ConstructClustering* construct = &defaultConstruct;
								if((imsgpils.problemType == ClusteringProblem::RCC_PROBLEM) and (imsgpils.k < 0)) {
										construct = &noConstruct;
								}

								if(imsgpils.cudaEnabled) {
									bestClustering = CUDAILS.executeILS(construct, &vnd, &sg, imsgpils.iter,
											imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
											problemFactory.build(imsgpils.problemType, imsgpils.k), info);
									timeSpent = CUDAILS.getTotalTimeSpent();
								} else {
									bestClustering = resolution.executeILS(construct, &vnd, &sg, imsgpils.iter,
											imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
											problemFactory.build(imsgpils.problemType, imsgpils.k), info);
									timeSpent = resolution.getTotalTimeSpent();
								}

								// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
								std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(sg.graph);
								for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
									globalVertexId.push_back(sg.graph.local_to_global(*it));
								}
							}
						} else {
							if(imsgpils.cudaEnabled) {
								bestClustering = CUDAILS.executeILS(construct, &vnd, g.get(), imsgpils.iter,
														imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
														problemFactory.build(imsgpils.problemType, imsgpils.k), info);
								timeSpent = CUDAILS.getTotalTimeSpent();
							} else {
								bestClustering = resolution.executeILS(construct, &vnd, g.get(), imsgpils.iter,
														imsgpils.iterMaxILS, imsgpils.perturbationLevelMax,
														problemFactory.build(imsgpils.problemType, imsgpils.k), info);
								timeSpent = resolution.getTotalTimeSpent();
							}
						}

						// Sends the result back to the leader process
						OutputMessage omsg(bestClustering, CUDAILS.getNumberOfTestedCombinations(), timeSpent, globalVertexId);
						omsg.num_vertices = n;
						omsg.num_edges = e;
						world.send(MPIMessage::LEADER_ID, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": ILS Output Message sent to leader (size = " <<
								(sizeof(omsg)/1024.0) << " KB).";
					}
				}
			} else {  // Parallel VND slave
				BOOST_LOG_TRIVIAL(info) << "Process " << myRank << " ready [VND slave process].";
				while(true) {
					InputMessageParallelVND imsgvns;
                    mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, imsgvns);
					if(stat.tag() == MPIMessage::TERMINATE_MSG_TAG) {
						BOOST_LOG_TRIVIAL(info) << "Process " << myRank << ": terminate message received.";
						return 0;
					} else if(stat.tag() == MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG) {
						// ignores interrupt message from previous VND processing
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": ignoring VND interrupt message received.";
						continue;
					} else if(stat.tag() == MPIMessage::INPUT_MSG_PARALLEL_VND_TAG) {
						// Receives a parallel VND message with parameters and triggers local VND execution
						// receives a message of type ParallelGrasp::INPUT_MSG_PARALLEL_VND_TAG or a terminate msg
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << " [Parallel VND]: Received message from master." << endl;
						messageCount++;

						if(imsgvns.l == 0) {
							BOOST_LOG_TRIVIAL(fatal) << "ERROR: Empty VND message received. Terminating.\n";
							return 1;
						}
						// builds a new graph object only if it has changed (saves time)
						if(previousId != imsgvns.id) {
							g.reset();
							g = reader.readGraphFromFile(imsgvns.graphInputFilePath);
							previousId = imsgvns.id;
						}
						// global vertex id array used in split graph
						std::vector<long> globalVertexId;

						// triggers the local partial VND search to be done between initial and final cluster indices
						// TODO: include firstImprovementOnOneNeig as a parameter inside the MPI message
						bool firstImprovementOnOneNeig = true;
						ParallelNeighborhoodSearch pnSearch(mpiparams.machineProcessAllocationStrategy, mpiparams.numberOfMasters, mpiparams.numberOfSearchSlavesPerMaster);
						Clustering bestClustering = pnSearch.searchNeighborhood(imsgvns.l, g.get(), &imsgvns.clustering,
								problemFactory.build(imsgvns.problemType, imsgvns.k), imsgvns.timeSpentSoFar, imsgvns.timeLimit, seed, myRank,
								imsgvns.initialClusterIndex, imsgvns.finalClusterIndex, firstImprovementOnOneNeig, imsgvns.k);

						// Sends the result back to the master process
						int leader = stat.source();
						OutputMessage omsg(bestClustering, pnSearch.getNumberOfTestedCombinations(), 0.0, globalVertexId);
						world.send(leader, MPIMessage::OUTPUT_MSG_PARALLEL_VND_TAG, omsg);
						BOOST_LOG_TRIVIAL(debug) << "Process " << myRank << ": VND Output Message sent to master number "
								<< leader << ". I(P) = " << bestClustering.getImbalance().getValue();
					}
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

void CommandLineInterfaceController::terminateMPIProcessesIfAny(int np, int machineProcessAllocationStrategy,
		int numberOfMasters, int numberOfSearchSlavesPerMaster) {
	BOOST_LOG_TRIVIAL(info) << "Terminating MPI slave processes...";
	mpi::communicator world;
	if(np > 1) {
		InputMessageParallelGrasp imsgpgrasp;
		InputMessageParallelVND imsgpvns;
		for(int i = 1; i < np; i++) {
			if(MPIUtil::isMaster(machineProcessAllocationStrategy, i, numberOfMasters, numberOfSearchSlavesPerMaster)) {
				// GRASP or ILS slave
				world.send(i, MPIMessage::TERMINATE_MSG_TAG, imsgpgrasp);
			} else {  // VND slave
				world.send(i, MPIMessage::TERMINATE_MSG_TAG, imsgpvns);
			}
			BOOST_LOG_TRIVIAL(info) << "Terminate message sent to process " << i << "\n";
		}
	}
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
