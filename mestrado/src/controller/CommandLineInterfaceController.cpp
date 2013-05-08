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
#include "../graph/include/Clustering.h"
#include "../problem/include/ClusteringProblem.h"
#include "../problem/include/CCProblem.h"

#include <boost/program_options.hpp>

using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
using namespace std;
using namespace clusteringgraph;
using namespace resolution::grasp;
using namespace problem;

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

int CommandLineInterfaceController::processArgumentsAndExecute(int argc, char *argv[]) {
    cout << "Correlation clustering problem solver" << endl << endl;

	try {
        float alpha = 0.5F;
        int numberOfIterations = 500, k = -1, l = 1;
        bool debug = false, RCC = false;
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

        if (vm.count("input-file"))
        {
            cout << "Input files are: "
                 << vm["input-file"].as< vector<string> >() << "\n";
        } else {
        	cout << "Please specify at least one input file.";
        	return 1;
        }

        if (k != -1) {
        	cout << "k value is " << k << ". RCC is enabled." << endl;
        	RCC = true;
        }
        cout << "Alpha value is " << std::setprecision(2) << fixed << alpha << "\n";
        cout << "Neighborhood size (l) is " << l << "\n";
        cout << "Number of iterations is " << numberOfIterations << "\n";
        cout << "Resolution strategy is " << strategy << endl;

        // Reads the graph from the specified text file
        SimpleTextGraphFileReader reader;
        SignedGraph* g = reader.readGraphFromFile(vm["input-file"].as< vector<string> >().at(0));
        if(debug) {
        	g->printGraph();
        }
        // Triggers the execution of the GRASP algorithm
        Grasp resolution;
        // TODO resolver problema de referencia do objeto CCProblem
        CCProblem* problem = new problem::CCProblem();
        Clustering* c = resolution.executeGRASP(g, numberOfIterations, alpha, l, problem);
        c->printClustering();
    }
    catch(std::exception& e)
    {
        cout << "Abnormal program termination. Stracktrace: " << endl;
    	cout << e.what() << "\n";
        return 1;
    }
    return 0;
}

} // namespace controller
