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

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_matrix.hpp>

using namespace std;
using namespace boost;

namespace controller {

SimpleTextGraphFileReader::SimpleTextGraphFileReader() {
	// TODO Auto-generated constructor stub

}

SimpleTextGraphFileReader::~SimpleTextGraphFileReader() {
	// TODO Auto-generated destructor stub
}

SignedGraph* SimpleTextGraphFileReader::readGraphFromFile(string filepath) {

	int n = 0, e = 0;
	ifstream in(filepath.c_str());
	if (!in.is_open()) return NULL;

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;

	vector< string > vec;
	string line;

	// captura a primeira linha do arquivo contendo as informacoes
	// de numero de vertices e arestas do grafo
	getline(in,line);
	char_separator<char> sep(" \r\n");
	tokenizer< char_separator<char> > tokens(line, sep);
	vec.assign(tokens.begin(),tokens.end());
	std::cout << line << std::endl;
	try {
	    n = boost::lexical_cast<int>(vec.at(0));
	    e = boost::lexical_cast<int>(vec.at(1));
	} catch( boost::bad_lexical_cast const& ) {
	    std::cerr << "Error: input string was not valid" << std::endl;
	}

	SignedGraph* g = new SignedGraph(n);
	std::cout << "Successfully created signed graph with " << n << " vertices." << std::endl;

	// captura as arestas do grafo com seus valores
	while (getline(in,line))
	{
		char_separator<char> sep(" \r\n");
		tokenizer< char_separator<char> > tokens(line, sep);
		vec.assign(tokens.begin(),tokens.end());
		if (vec.size() < 3) continue;
		if(vec.at(2).rfind('\n') != string::npos)
		std::cout << vec.at(0) << vec.at(1) << vec.at(2) << "/" << std::endl;

		try {
			int a = boost::lexical_cast<int>(vec.at(0));
			int b = boost::lexical_cast<int>(vec.at(1));
			int value = boost::lexical_cast<int>(vec.at(2));
			// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
			g->addEdge(a, b, value);
		} catch( boost::bad_lexical_cast const& ) {
			std::cerr << "Error: input string was not valid" << std::endl;
		}
	}

	return g;
}

} /* namespace controller */
