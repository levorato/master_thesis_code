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

#include <boost/smart_ptr.hpp>
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

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromString(const string& graphContents) {
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	int n = 0, e = 0;
	vector< string > lines;

	// captura a primeira linha do arquivo contendo as informacoes
	// de numero de vertices e arestas do grafo
	char_separator<char> sep("\r\n");
	char_separator<char> sep2(" ");
	tokenizer< char_separator<char> > tokens(graphContents, sep);
	lines.assign(tokens.begin(),tokens.end());

	try {
		string line = lines.back();
		lines.pop_back();
		tokenizer< char_separator<char> > tokens2(line, sep2);
		vector<string> vec;
		vec.assign(tokens2.begin(),tokens2.end());

	    n = boost::lexical_cast<int>(vec.at(0));
	    e = boost::lexical_cast<int>(vec.at(1));
	} catch( boost::bad_lexical_cast const& ) {
	    std::cerr << "Error: input string was not valid" << std::endl;
	}

	SignedGraphPtr g = boost::make_shared<SignedGraph>(n);
	std::cout << "Successfully created signed graph with " << n << " vertices." << std::endl;

	// captura as arestas do grafo com seus valores
	while (not lines.empty())
	{
		string line = lines.back();
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
			// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
			g->addEdge(a, b, value);
		} catch( boost::bad_lexical_cast const& ) {
			std::cerr << "Error: input string was not valid" << std::endl;
		}
	}
	g->setGraphAsText(graphContents);

	return g;
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromFile(const string& filepath) {
	cout << "Reading input file: '" << filepath << "' ..." << endl;
	ifstream in(filepath.c_str());
	if (!in.is_open()) throw "Cannot open file " + filepath;
	std::string str((std::istreambuf_iterator<char>(in)),
	                 std::istreambuf_iterator<char>());
	return readGraphFromString(str);
}

} /* namespace controller */
