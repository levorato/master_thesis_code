/*
 * TestController.cpp
 *
 *  Created on: Oct 13, 2015
 *      Author: mario
 */

#include "include/TestController.h"

#include "graph/include/SimpleTextGraphFileReader.h"
#include "graph/include/Graph.h"
#include "../../resolution/grasp/include/ParallelGrasp.h"
#include "../../resolution/ils/include/ILS.h"
#include "../../resolution/ils/include/ParallelILS.h"
#include "graph/include/Clustering.h"
#include "graph/include/Imbalance.h"
#include "problem/include/ClusteringProblem.h"
#include "problem/include/CCProblem.h"
#include "util/include/TimeDateUtil.h"
#include "util/include/EnumUtil.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "problem/include/ClusteringProblemFactory.h"
#include "../../resolution/construction/include/GainFunctionFactory.h"
#include "../../resolution/construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"
#include "../../resolution/construction/include/ConstructClustering.h"
#include "util/include/RandomUtil.h"
#include "../../resolution/vnd/include/CUDAVariableNeighborhoodDescent.h"
#include "../../resolution/vnd/include/CUDANeighborhoodSearch.h"
#include "../../resolution/grasp/include/CUDAGrasp.h"
#include "../../resolution/construction/include/CUDAConstructClustering.h"
#include "../../resolution/construction/include/CUDAImbalanceGainFunction.h"
#include "../../resolution/ils/include/CUDAILS.h"

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
#include <boost/mpi/communicator.hpp>

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

std::istream& operator>>(std::istream& in, TestController::StategyName& strategy)
{
    std::string token;
    in >> token;
    if (token == "GRASP")
    	strategy = TestController::GRASP;
    else if (token == "ILS")
    	strategy = TestController::ILS;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

std::istream& operator>>(std::istream& in, TestController::SearchName& search)
{
    std::string token;
    in >> token;
    if (token == "SEQUENTIAL")
        search = TestController::SEQUENTIAL_SEARCH;
    else if (token == "PARALLEL")
        search = TestController::PARALLEL_SEARCH;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

TestController::TestController() {
	// TODO Auto-generated constructor stub

}

TestController::~TestController() {
	// TODO Auto-generated destructor stub
}

void TestController::testSubgraphCreationUndirectedGraph(string executionId, unsigned long seed) {

	// Application Parameters:
	// -l 1 --iter=1 --alpha=1.0
	// --input-file "/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g"
	// --output-folder "../output/cuda-ils/sdot-split-8p-l1-i1-vpartition" --gain-function-type=0 --time-limit=1800
	// --rcc=false --strategy ILS --split=true
	mpi::communicator world;
	int np = world.size();
	int myRank = 0;
	int numberOfIterations = 1, l = 1;
	double alpha = 1.0;
	long k = 0;
	bool debug = false, profile = false, firstImprovementOnOneNeig = true, CCEnabled = true, RCCEnabled = false;
	string inputFileDir;
	string outputFolder("../output/test-suite/testSubgraphCreationUndirectedGraph");
	int timeLimit = 1800;
	int functionType = GainFunction::IMBALANCE;
	int numberOfMasters = np - 1;
	int totalNumberOfVNDSlaves = 0;
	string jobid;
	TestController::StategyName strategy = TestController::ILS;
	TestController::SearchName searchType = TestController::SEQUENTIAL_SEARCH;
	int iterMaxILS = 5, perturbationLevelMax = 30;  // for Slashdot and completely Random instances
	bool splitGraph = true;
	bool cuda = true;
	// fs::path filePath("/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g");
	fs::path filePath("/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/directed/slashdot.g");

	if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
		// Reads the graph from the specified text file
		UndirectedGraph g = readGraphFromFile(filePath.string());
		string fileId = filePath.filename().string();
		boost::timer::cpu_timer timer;
		boost::timer::cpu_times start_time, end_time;
		double timeSpent = 0.0;
		string filename;

		// Divides the set of vertices V (size = n) into 'np' random subsets of the same size
		BOOST_LOG_TRIVIAL(info) << "Partitioning vertices...";
		std::vector<long> vertexList;
		long n = num_vertices(g);
		for(long i = 0; i < n; i++) {
			vertexList.push_back(i);
		}
		// using built-in random generator
		std::random_shuffle(vertexList.begin(), vertexList.end());
		long desiredCardinality = long(floor(n / (double)np));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " vertices.";
		std::vector< std::vector< long > > verticesInCluster(np, std::vector< long >());
		for(long i = 0, k = 0; i < n; i++) {
			verticesInCluster[k].push_back(vertexList[i]);
			if((((i + 1) % desiredCardinality) == 0) and (k < np - 1)) {  k++;  }
		}
		BOOST_LOG_TRIVIAL(info) << "Cardinality of last partition is " << verticesInCluster[np - 1].size() << " vertices.";
		BOOST_LOG_TRIVIAL(info) << "Full graph data: n =  " << num_vertices(g) << ", " << "e =  " << num_edges(g);
		BOOST_LOG_TRIVIAL(info) << "[Master process] Creating undirected subgraphs...";

		for(long k = 0; k < np; k++) {
			// time measure of subgraph creation
			timer.start();
			start_time = timer.elapsed();

			// SignedGraph sg(g->graph, verticesInCluster[k]);
			UndirectedGraph& sg = g.create_subgraph(verticesInCluster[k].begin(), verticesInCluster[k].end());
			BOOST_LOG_TRIVIAL(info) << "Created subgraph " << k << " with n =  " << num_vertices(sg) << ", " << "e =  " << num_edges(sg);

			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			BOOST_LOG_TRIVIAL(info) << "Time spent: " << std::fixed << std::setprecision(2) << timeSpent << " s";
		}


	} else {
		BOOST_LOG_TRIVIAL(fatal) << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
	}


}

void TestController::testSubgraphCreationPerformance(string executionId, unsigned long seed) {

	// Application Parameters:
	// -l 1 --iter=1 --alpha=1.0
	// --input-file "/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g"
	// --output-folder "../output/cuda-ils/sdot-split-8p-l1-i1-vpartition" --gain-function-type=0 --time-limit=1800
	// --rcc=false --strategy ILS --split=true
	mpi::communicator world;
	int np = world.size();
	int myRank = 0;
	int numberOfIterations = 1, l = 1;
	double alpha = 1.0;
	long k = 0;
	bool debug = false, profile = false, firstImprovementOnOneNeig = true, CCEnabled = true, RCCEnabled = false;
	string inputFileDir;
	string outputFolder("../output/test-suite/testSubgraphCreationPerformance");
	int timeLimit = 1800;
	int functionType = GainFunction::IMBALANCE;
	int numberOfMasters = np - 1;
	int totalNumberOfVNDSlaves = 0;
	string jobid;
	TestController::StategyName strategy = TestController::ILS;
	TestController::SearchName searchType = TestController::SEQUENTIAL_SEARCH;
	int iterMaxILS = 5, perturbationLevelMax = 30;  // for Slashdot and completely Random instances
	bool splitGraph = true;
	bool cuda = true;
	// fs::path filePath("/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g");
	fs::path filePath("/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/directed/slashdot.g");

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
		int machineProcessAllocationStrategy = MPIUtil::MASTER_AND_VND_SLAVES_TOGETHER;
		int numberOfSearchSlaves = 0;
		NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves);
		if(searchType == TestController::PARALLEL_SEARCH) {
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

		// Divides the set of vertices V (size = n) into 'np' random subsets of the same size
		BOOST_LOG_TRIVIAL(info) << "Partitioning vertices...";
		std::vector<long> vertexList;
		long n = g->getN();
		for(long i = 0; i < n; i++) {
			vertexList.push_back(i);
		}
		// using built-in random generator
		std::random_shuffle(vertexList.begin(), vertexList.end());
		long desiredCardinality = long(floor(n / (double)np));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " vertices.";
		std::vector< std::vector< long > > verticesInCluster(np, std::vector< long >());
		for(long i = 0, k = 0; i < n; i++) {
			verticesInCluster[k].push_back(vertexList[i]);
			if((((i + 1) % desiredCardinality) == 0) and (k < np - 1)) {  k++;  }
		}
		BOOST_LOG_TRIVIAL(info) << "Cardinality of last partition is " << verticesInCluster[np - 1].size() << " vertices.";
		BOOST_LOG_TRIVIAL(info) << "Full graph data: n =  " << num_vertices(g->graph) << ", " << "e =  " << num_edges(g->graph);
		BOOST_LOG_TRIVIAL(info) << "[Master process] Creating subgraphs...";

		for(long k = 0; k < np; k++) {
			// time measure of subgraph creation
			timer.start();
			start_time = timer.elapsed();

			SignedGraph sg(&(g->graph), verticesInCluster[k]);
			BOOST_LOG_TRIVIAL(info) << "Created subgraph " << k << " with n =  " << num_vertices(sg.graph) << ", " << "e =  " << num_edges(sg.graph);

			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			BOOST_LOG_TRIVIAL(info) << "Time spent: " << std::fixed << std::setprecision(2) << timeSpent << " s";
		}


	} else {
		BOOST_LOG_TRIVIAL(fatal) << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
	}
}

void TestController::testSplitGraphParallelILS(string executionId, unsigned long seed) {

	// Application Parameters:
	// -l 1 --iter=1 --alpha=1.0
	// --input-file "/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g"
	// --output-folder "../output/cuda-ils/sdot-split-8p-l1-i1-vpartition" --gain-function-type=0 --time-limit=1800
	// --rcc=false --strategy ILS --split=true
	mpi::communicator world;
	int np = world.size();
	int myRank = 0;
	int numberOfIterations = 1, l = 1;
	double alpha = 1.0;
	long k = 0;
	bool debug = false, profile = false, firstImprovementOnOneNeig = true, CCEnabled = true, RCCEnabled = false;
	string inputFileDir;
	string outputFolder("../output/test-suite/testSplitGraphParallelILS");
	int timeLimit = 1800;
	int functionType = GainFunction::IMBALANCE;
	int numberOfMasters = np - 1;
	int totalNumberOfVNDSlaves = 0;
	string jobid;
	TestController::StategyName strategy = TestController::ILS;
	TestController::SearchName searchType = TestController::SEQUENTIAL_SEARCH;
	int iterMaxILS = 5, perturbationLevelMax = 30;  // for Slashdot and completely Random instances
	bool splitGraph = true;
	bool cuda = true;
	fs::path filePath("/home/mlevorato/mestrado-cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0/bigger/slashdot-undirected-size10000-part0.g");

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
		int machineProcessAllocationStrategy = MPIUtil::MASTER_AND_VND_SLAVES_TOGETHER;
		int numberOfSearchSlaves = 0;
		NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, numberOfMasters, numberOfSearchSlaves);
		if(searchType == TestController::PARALLEL_SEARCH) {
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
		// medicao de tempo do CC
		timer.start();
		start_time = timer.elapsed();

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
	} else {
		BOOST_LOG_TRIVIAL(fatal) << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
	}
}

// http://www.concentric.net/~Ttwang/tech/inthash.htm
unsigned long TestController::mix(unsigned long a, unsigned long b, unsigned long c)
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

int TestController::runTestSuite() {

	// used for debugging purpose
	std::set_terminate( handler );
	mpi::communicator world;

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

	// job id is obtained through command line parameter from PBS Scheduler
	//cout << "Job id is " << jobid << "\n";
	// initializes the logging subsystem
	string jobid = TimeDateUtil::generateRandomId();
	
	unsigned int numberOfSearchSlavesPerMaster = 0;
	int machineProcessAllocationStrategy = 0;
	// id used for output folders
	string executionId = jobid;

	try {
		// TEST INVOCATIONS GO HERE
		this->testSubgraphCreationUndirectedGraph(executionId, seed);
		this->testSubgraphCreationPerformance(executionId, seed);

		BOOST_LOG_TRIVIAL(info) << "Test suite execution done.";
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

void TestController::handler()
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

UndirectedGraph TestController::readGraphFromString(const string& graphContents) {
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	int n = 0, e = 0;
	// defines the format of the input file
	int formatType = 0;
	std::vector< string > lines;

	// captura a primeira linha do arquivo contendo as informacoes
	// de numero de vertices e arestas do grafo
	char_separator<char> sep("\r\n");
	char_separator<char> sep2(" \t");
	tokenizer< char_separator<char> > tokens(graphContents, sep);
	lines.assign(tokens.begin(),tokens.end());

	try {
		string line = lines.at(0);
		lines.erase(lines.begin());
		trim(line);
		BOOST_LOG_TRIVIAL(trace) << "Line: " << line << endl;

		if(line.find("people") != string::npos) {  // xpress files
			BOOST_LOG_TRIVIAL(trace) << "Format type is 0" << endl;
			string firstLine = lines.at(0);
			string number = line.substr(line.find("people:") + 7);
			trim(number);
			n = boost::lexical_cast<int>(number);
			BOOST_LOG_TRIVIAL(trace) << "n value is " << n << endl;
			formatType = 0;
		} else if(line.find("Vertices") != string::npos) {
			BOOST_LOG_TRIVIAL(trace) << "Format type is 1" << endl;
			string number = line.substr(line.find("Vertices") + 8);
			trim(number);
			n = boost::lexical_cast<int>(number);
			BOOST_LOG_TRIVIAL(trace) << "n value is " << n << endl;
			while(lines.at(0).find("Arcs") == string::npos && lines.at(0).find("Edges") == string::npos) {
				lines.erase(lines.begin());
			}
			lines.erase(lines.begin());
			formatType = 1;
		} else {
			tokenizer< char_separator<char> > tokens2(line, sep2);
			std::vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			n = boost::lexical_cast<int>(vec.at(0));
			if(vec.size() == 1) { // .dat files
				BOOST_LOG_TRIVIAL(trace) << "Format type is 3" << endl;
				formatType = 3;
			} else {
				BOOST_LOG_TRIVIAL(trace) << "Format type is 2" << endl;
				e = boost::lexical_cast<int>(vec.at(1));
				formatType = 2;
			}
		}
	} catch( boost::bad_lexical_cast const& e ) {
	    std::cerr << "Error: input string was not valid" << std::endl;
	    cerr << e.what() << "\n";
	    BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
	}

	UndirectedGraph g(n);
	BOOST_LOG_TRIVIAL(debug) << "Successfully created undirected signed graph with " << n << " vertices." << std::endl;

	// captura as arestas do grafo com seus valores
	if(formatType == 2 || formatType == 1) {
		long imbalance = 0;
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep2);
			std::vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());

			if (vec.size() < 3) continue;
			//if(vec.at(2).rfind('\n') != string::npos)
			// BOOST_LOG_TRIVIAL(trace) << vec.at(0) << vec.at(1) << vec.at(2) << "/" << std::endl;

			try {
				int a = boost::lexical_cast<int>(vec.at(0));
				int b = boost::lexical_cast<int>(vec.at(1));
				double value = 0.0;
				if(vec.at(2) != "*") {
					sscanf(vec.at(2).c_str(), "%lf", &value);
					if(formatType == 2) {
						// vertex number must be in the interval 0 <= i < n
						if(a >= n or b >= n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
						}
						add_edge(a, b, value, g);
						// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
					} else {
						if(a > n or b > n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
						}
						add_edge(a - 1, b - 1, value, g);
						// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					}
				} else {  // special notation for directed edges (add 1 to imbalance)
					imbalance++;
				}
			} catch( boost::bad_lexical_cast const& ) {
				BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
			}

		}
	} else if(formatType == 0) {  // xpress files
		char_separator<char> sep3(" (),\t\r\n[]");
		string temp = graphContents.substr(graphContents.find("Mrel:") + 5);
		tokenizer< char_separator<char> > tokens2(temp, sep3);
		std::vector<string> vec2;
		vec2.assign(tokens2.begin(),tokens2.end());
		while(vec2.at(0).length() == 0) {
			vec2.erase(vec2.begin());
		}
		while(vec2.back().length() == 0) {
			vec2.pop_back();
		}

		int size = vec2.size();
		if (size % 3 != 0) {
			BOOST_LOG_TRIVIAL(fatal) << "Error: invalid XPRESS file format!" << std::endl;
		}
		for(int i = 0; i + 2 < size; i = i + 3) {
			try {
				// std::cout << "Processing line " << vec2.at(i) << " " << vec2.at(i+1) << " " << vec2.at(i+2) << std::endl;
				int a = boost::lexical_cast<int>(vec2.at(i));
				int b = boost::lexical_cast<int>(vec2.at(i + 1));
				double value = 0.0;
				sscanf(vec2.at(i + 2).c_str(), "%lf", &value);
				// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;

				if(a > n or b > n) {
					BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
					cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
				}
				add_edge(a - 1, b - 1, value, g);
				// g->addEdge(b - 1, a - 1, value);
			} catch( boost::bad_lexical_cast const& ) {
				BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
			}
		}
	} else {  // formatType == 3, .dat files
		char_separator<char> sep3(" \t");
		unsigned int a = 0;
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep3);
			std::vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			// cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

			for(unsigned int b = 0; b < vec.size(); b++) {
				try {
					double value = 0.0;
	                sscanf(vec.at(b).c_str(), "%lf", &value);
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					// the following is to avoid duplicate couting of arcs in the objective function
	                if(a >= n or b >= n) {
						BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
						cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
					}
	                add_edge(a, b, value, g);
				} catch( boost::bad_lexical_cast const& ) {
					BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
				}
			}
			a++;
		}
	}
	// g->printGraph();
	BOOST_LOG_TRIVIAL(info) << "Successfully read graph file.";

	return g;
}

UndirectedGraph TestController::readGraphFromFile(const string& filepath) {
	BOOST_LOG_TRIVIAL(info) << "Reading input file: '" << filepath << "' ...";
	return readGraphFromString(get_file_contents(filepath.c_str()));
}

std::string TestController::get_file_contents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return(contents.str());
  }
  throw(errno);
}

} // namespace controller
