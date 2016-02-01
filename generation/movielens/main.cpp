#include <boost/program_options.hpp>
using namespace boost::program_options;

#include <iostream>
using namespace std;

#include "MovieLensSGConverter.h"
using namespace generation;

int main(int argc, char* argv[])
{
    try {
        string folder, fileFilter = "ratings.dat";
        
        options_description desc("Convert MovieLens dataset file (ratings.dat) to unweighted signed graph.");
        desc.add_options()
        // First parameter describes option name/short name
        // The second is parameter to option
        // The third is description
        ("help,h", "print usage message")
        ("folder,f", value<string>(&folder), "the folder containing the ratings.dat files")
        ("filefilter,l", value<string>(&fileFilter)->default_value("ratings.dat"), 
         "the filename for MovieLens ratings dataset files (default: ratings.dat)")
        ;
    
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {  
            cout << desc << "\n";
            return 0;
        }

		cout << "Input folder is " << folder << endl;
        cout << "File filter is " << fileFilter << endl;
		
		MovieLensSGConverter converter;
		converter.processMovieLensFolder(folder, fileFilter);
    }
    catch(exception& e) {
        cerr << e.what() << "\n";
    }
}
