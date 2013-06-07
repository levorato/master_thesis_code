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

#include <boost/algorithm/string.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_matrix.hpp>

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
	char_separator<char> sep2(" ");
	tokenizer< char_separator<char> > tokens(graphContents, sep);
	lines.assign(tokens.begin(),tokens.end());

	try {
		string line = lines.at(0);
		lines.erase(lines.begin());
		trim(line);
		cout << "Line: " << line << endl;

		if(line.find("people") != string::npos) {
			cout << "Format type is 0" << endl;
			string number = line.substr(line.find(":") + 1);
			trim(number);
			n = boost::lexical_cast<int>(number);
			cout << "n value is " << n << endl;
			// ignore the next 2 lines of the file
			lines.erase(lines.begin());
			lines.erase(lines.begin());
			// reads the first line
			// removes the Mrel: [  from the first graph file line
			string firstLine = lines.at(0);
			lines.erase(lines.begin());
			lines.push_back(firstLine.substr(firstLine.find("[") + 1));
			formatType = 0;
		} else if(line.find("Vertices") != string::npos) {
			cout << "Format type is 1" << endl;
			string number = line.substr(line.find("Vertices") + 8);
			trim(number);
			n = boost::lexical_cast<int>(number);
			cout << "n value is " << n << endl;
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
			if(vec.size() == 1) {
				cout << "Format type is 3" << endl;
				formatType = 3;
			} else {
				cout << "Format type is 2" << endl;
				e = boost::lexical_cast<int>(vec.at(1));
				formatType = 2;
			}
		}
	} catch( boost::bad_lexical_cast const& ) {
	    std::cerr << "Error: input string was not valid" << std::endl;
	}

	SignedGraphPtr g = boost::make_shared<SignedGraph>(n);
	std::cout << "Successfully created signed graph with " << n << " vertices." << std::endl;

	// captura as arestas do grafo com seus valores
	if(formatType == 2 || formatType == 1) {
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep2);
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());

			if (vec.size() < 3) continue;
			if(vec.at(2).rfind('\n') != string::npos)
				std::cout << vec.at(0) << vec.at(1) << vec.at(2) << "/" << std::endl;

			try {
				int a = boost::lexical_cast<int>(vec.at(0));
				int b = boost::lexical_cast<int>(vec.at(1));
				int value = boost::lexical_cast<int>(vec.at(2));
				if(formatType == 2) {
					g->addEdge(a, b, value);
					// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
				} else {
					g->addEdge(a - 1, b - 1, value);
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
				}
			} catch( boost::bad_lexical_cast const& ) {
				std::cerr << "Error: input string was not valid" << std::endl;
			}
		}
	} else if(formatType == 0) {
		char_separator<char> sep3(" (),");
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			if(line.find("]") != string::npos)  continue;
			tokenizer< char_separator<char> > tokens2(line, sep3);
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			// cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

			if (vec.size() < 3) continue;
			if(vec.at(2).rfind('\n') != string::npos)
				std::cout << vec.at(0) << vec.at(1) << vec.at(2) << "/" << std::endl;

			try {
				int a = boost::lexical_cast<int>(vec.at(0));
				int b = boost::lexical_cast<int>(vec.at(1));
				int value = boost::lexical_cast<int>(vec.at(2));
				// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
				g->addEdge(a - 1, b - 1, value);
				// g->addEdge(b - 1, a - 1, value);
			} catch( boost::bad_lexical_cast const& ) {
				std::cerr << "Error: input string was not valid" << std::endl;
			}
		}
	} else {  // formatType == 3, .dat files
		char_separator<char> sep3(" ");
		int a = 0;
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
					int value = boost::lexical_cast<int>(vec.at(b));
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					g->addEdge(a, b, value);
				} catch( boost::bad_lexical_cast const& ) {
					std::cerr << "Error: input string was not valid" << std::endl;
				}
			}
			a++;
		}
	}
	g->printGraph();
	g->setGraphAsText(graphContents);

	return g;
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromFile(const string& filepath) {
	cout << "Reading input file: '" << filepath << "' ..." << endl;
	return readGraphFromString(get_file_contents(filepath.c_str()));
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
