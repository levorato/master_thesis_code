/*
 * SimpleTextGraphFileReader.cpp
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#include "include/SimpleTextGraphFileReader.h"
#include "../graph/include/Graph.h"
#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <streambuf>
#include <cerrno>
#include <cstdio>

#include <boost/algorithm/string.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/log/trivial.hpp>
#include <boost/functional/hash.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;

namespace controller {

SimpleTextGraphFileReader::SimpleTextGraphFileReader() {
	// TODO Auto-generated constructor stub

}

SimpleTextGraphFileReader::~SimpleTextGraphFileReader() {
	// TODO Auto-generated destructor stub
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromString(const string& graphContents) {
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	int n = 0, e = 0;
	// defines the format of the input file
	int formatType = 0;
	vector< string > lines;

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
			vector<string> vec;
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

	SignedGraphPtr g = boost::make_shared<SignedGraph>(n);
	BOOST_LOG_TRIVIAL(trace) << "Successfully created signed graph with " << n << " vertices." << std::endl;

	// captura as arestas do grafo com seus valores
	if(formatType == 2 || formatType == 1) {
		long imbalance = 0;
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep2);
			vector<string> vec;
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
						g->addEdge(a, b, value);
						// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
					} else {
						g->addEdge(a - 1, b - 1, value);
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
		vector<string> vec2;
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
				g->addEdge(a - 1, b - 1, value);
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
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			// cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

			for(unsigned int b = 0; b < vec.size(); b++) {
				try {
					double value = 0.0;
	                sscanf(vec.at(b).c_str(), "%lf", &value);
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					// the following is to avoid duplicate couting of arcs in the objective function
					g->addEdge(a, b, value);
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

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromFile(const string& filepath) {
	BOOST_LOG_TRIVIAL(info) << "Reading input file: '" << filepath << "' ..." << endl;
	SignedGraphPtr g = readGraphFromString(get_file_contents(filepath.c_str()));
	g->setGraphFileLocation(filepath);
	boost::hash<std::string> string_hash;
	g->setId(string_hash(filepath));

	return g;
}

std::string SimpleTextGraphFileReader::get_file_contents(const char *filename)
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

} /* namespace controller */
