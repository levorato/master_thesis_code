/*
 * CommandLineInterfaceController.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: mario
 */

#include "include/CommandLineInterfaceController.h"
#include "include/SimpleTextGraphFileReader.h"

#include <boost/program_options.hpp>

using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
using namespace std;

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

int CommandLineInterfaceController::processArguments(int argc, char *argv[]) {
    cout << "Correlation clustering problem solver" << endl << endl;

	try {
        float alpha;
        int numberOfIterations;
        CommandLineInterfaceController::StategyName strategy = CommandLineInterfaceController::GRASP;

        po::options_description desc("Available options:");
        desc.add_options()
            ("help", "show program options")
            ("alpha", po::value<float>(&alpha)->default_value(0.5F),
                  "alpha - randomness factor")
            ("iterations,it", po::value<int>(&numberOfIterations)->default_value(500),
                  "number of iterations")
            ("input-file", po::value< vector<string> >(), "input file")
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
        }

        cout << "Alpha value is " << std::setprecision(2) << fixed << alpha << "\n";
        cout << "Number of iterations is " << numberOfIterations << "\n";
        cout << "Resolution strategy is " << strategy << endl;

        SimpleTextGraphFileReader::readGraphFromFile(vm["input-file"].at(0));
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
