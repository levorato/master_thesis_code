/*
 * SimpleTextGraphFileReader.cpp
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#include "include/SimpleTextGraphFileReader.h"
#include "../graph/Graph.h"
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

*SignedGraph SimpleTextGraphFileReader::readGraphFromFile(string filepath) {

	int n = 0, e = 0;
	ifstream in(filepath.c_str());
	if (!in.is_open()) return NULL;

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;

	vector< string > vec;
	string line;

	// captura a primeira linha do arquivo contendo as informacoes
	// de numero de vertices e arestas do grafo
	getline(in,line);
	Tokenizer tok(line);
	vec.assign(tok.begin(),tok.end());
	try {
	    n = boost::lexical_cast<int>(vec.at(0));
	    e = boost::lexical_cast<int>(vec.at(1));
	} catch( boost::bad_lexical_cast const& ) {
	    std::cout << "Error: input string was not valid" << std::endl;
	}

	g = new SignedGraph(n);

	// captura as arestas do grafo com seus valores
	while (getline(in,line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(),tok.end());

		if (vec.size() < 3) continue;

		try {
			int a = boost::lexical_cast<int>(vec.at(0));
			int b = boost::lexical_cast<int>(vec.at(1));
			int value = boost::lexical_cast<int>(vec.at(2));
		} catch( boost::bad_lexical_cast const& ) {
			std::cout << "Error: input string was not valid" << std::endl;
		}

		g->addEdge(a, b, value);
	}

	g->printGraph();

	return g;
}

} /* namespace controller */
