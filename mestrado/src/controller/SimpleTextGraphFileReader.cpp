/*
 * SimpleTextGraphFileReader.cpp
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#include "include/SimpleTextGraphFileReader.h"
#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator

#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

namespace controller {

SimpleTextGraphFileReader::SimpleTextGraphFileReader() {
	// TODO Auto-generated constructor stub

}

SimpleTextGraphFileReader::~SimpleTextGraphFileReader() {
	// TODO Auto-generated destructor stub
}

static Graph SimpleTextGraphFileReader::readGraphFromFile(string filepath) {

	Graph g;
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
	n = stoi(vec.at(0));
	e = stoi(vec.at(1));

	// captura as arestas do grafo com seus valores
	while (getline(in,line))
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(),tok.end());

		if (vec.size() < 3) continue;

		int a = stoi(vec.at(0));
		int b = stoi(vec.at(1));
		int value = stoi(vec.at(2));

		add_edge(a, b, value);
	}

	std::cout << "vertex set: ";
	boost::print_vertices(g, name);
	std::cout << std::endl;

	std::cout << "edge set: ";
	boost::print_edges(g, name);
	std::cout << std::endl;

	std::cout << "out-edges: " << std::endl;
	boost::print_graph(g, name);
	std::cout << std::endl;

	return g;
}

} /* namespace controller */
