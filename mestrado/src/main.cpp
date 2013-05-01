//============================================================================
// Name        : main.cpp
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Program entry point. Includes main function and program options.
//============================================================================

#include <boost/program_options.hpp>

using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
using namespace std;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

namespace Resolution { enum StategyName {GRASP, GRASP_PR};};

std::istream& operator>>(std::istream& in, Resolution::StategyName& strategy)
{
    std::string token;
    in >> token;
    if (token == "GRASP")
        strategy = Resolution::GRASP;
    else if (token == "GRASP_PR")
        strategy = Resolution::GRASP_PR;
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

int main(int ac, char* av[])
{
    cout << "Correlation clustering problem solver" << endl << endl;

	try {
        float alpha;
        int numberOfIterations;
        enum Resolution::StategyName strategy;

        po::options_description desc("Available options:");
        desc.add_options()
            ("help", "show program options")
            ("alpha", po::value<float>(&alpha)->default_value(0.5F),
                  "alpha - ramdomness factor")
            ("iterations,it", po::value<int>(&numberOfIterations)->default_value(500),
                  "number of iterations")
            ("input-file", po::value< vector<string> >(), "input file")
            ("strategy",
                         po::value<Resolution::StategyName>(&strategy)->default_value(Resolution::GRASP),
                         "Resolution strategy to be used. Accepted values: GRASP, GRASP_PR.")
        ;

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
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
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
    return 0;
}
