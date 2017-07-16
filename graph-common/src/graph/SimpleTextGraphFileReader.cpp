/*
 * SimpleTextGraphFileReader.cpp
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#include "include/SimpleTextGraphFileReader.h"
#include "include/ParallelBGLSignedGraph.h"
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
#include <boost/graph/graphviz.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;

namespace controller {

SimpleTextGraphFileReader::~SimpleTextGraphFileReader() {
	// TODO Auto-generated destructor stub
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromFilepath(const string& filepath, const bool& parallelgraph) {
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	long n = 0, e = 0;
	// defines the format of the input file
	int formatType = 0;
	std::ifstream infile(filepath.c_str(), std::ios::in | std::ios::binary);
	SignedGraphPtr g;
	if (infile) {
		BOOST_LOG_TRIVIAL(info) << "Reading input file line by line, avoiding full file reading.";
		// captura a primeira linha do arquivo contendo as informacoes
		// de numero de vertices e arestas do grafo
		string line;
		std::getline(infile, line);

		char_separator<char> sep2(" \t");
		// BOOST_LOG_TRIVIAL(trace) << "Line read: " << line;
		try {
			if(line.find("%%MatrixMarket") != string::npos) {  // matrix market files
				BOOST_LOG_TRIVIAL(info) << "Format type is 4 (matrix market)" << endl;
				std::getline(infile, line);
				trim(line);
				// captura as dimensoes da matriz e o numero de arestas
				tokenizer< char_separator<char> > tokens2(line, sep2);
				vector<string> vec;
				vec.assign(tokens2.begin(),tokens2.end());
				n = boost::lexical_cast<long>(vec.at(0));
				formatType = 4;
			} else if(line.find("people") != string::npos) {  // xpress files
				BOOST_LOG_TRIVIAL(trace) << "Format type is 0 (xpress)" << endl;
				string number = line.substr(line.find("people:") + 7);
				trim(number);
				n = boost::lexical_cast<long>(number);
				BOOST_LOG_TRIVIAL(trace) << "n value is " << n << endl;
				formatType = 0;
			} else if(line.find("Vertices") != string::npos) {
				BOOST_LOG_TRIVIAL(trace) << "Format type is 1 (pajek)" << endl;
				string number = line.substr(line.find("Vertices") + 8);
				trim(number);
				n = boost::lexical_cast<long>(number);
				BOOST_LOG_TRIVIAL(trace) << "n value is " << n << endl;
				formatType = 1;
			} else {
				tokenizer< char_separator<char> > tokens2(line, sep2);
				vector<string> vec;
				vec.assign(tokens2.begin(),tokens2.end());
				n = boost::lexical_cast<long>(vec.at(0));
				if(vec.size() == 1) { // .dat files
					BOOST_LOG_TRIVIAL(trace) << "Format type is 3 (dat)" << endl;
					formatType = 3;
				} else {
					BOOST_LOG_TRIVIAL(trace) << "Format type is 2 (.g)" << endl;
					BOOST_LOG_TRIVIAL(trace) << "vec.at(0) = " << vec.at(0);
					//BOOST_LOG_TRIVIAL(trace) << "vec.at(1) = " << vec.at(1);
					// e = boost::lexical_cast<long>(vec.at(1));
					formatType = 2;
					// BOOST_LOG_TRIVIAL(trace) << "Num of edges is e = " << e;
				}
			}
		} catch( boost::bad_lexical_cast const& e ) {
			std::cerr << "Error: input string was not valid" << std::endl;
			cerr << e.what() << "\n";
			BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
		}

		BOOST_LOG_TRIVIAL(trace) << "Creating graph with size n = " << n;

		// All processes synchronize at this point, then the graph is complete
		// synchronize(grf.process_group());
		BOOST_LOG_TRIVIAL(trace) << "Successfully created distributed signed graph with " << n << " vertices.";

		// Only process 0 loads the graph, which is distributed automatically
		g = boost::make_shared<ParallelBGLSignedGraph>(n, graph);

		if(formatType == 1) {
			std::getline(infile, line);
			while(line.find("Arcs") == string::npos and line.find("Edges") == string::npos) {
				std::getline(infile, line);
			}
		}

		// captura as arestas do grafo com seus valores
		if(formatType == 2 || formatType == 1) {
			long imbalance = 0;
			while (std::getline(infile, line)) {
				trim(line);
				tokenizer< char_separator<char> > tokens2(line, sep2);
				vector<string> vec;
				vec.assign(tokens2.begin(),tokens2.end());

				if (vec.size() < 3) continue;
				//if(vec.at(2).rfind('\n') != string::npos)
				// BOOST_LOG_TRIVIAL(trace) << vec.at(0) << "; " << vec.at(1) << "; " << vec.at(2) << "/" << std::endl;

				try {
					long a = boost::lexical_cast<long>(vec.at(0));
					long b = boost::lexical_cast<long>(vec.at(1));
					double value = 0.0;
					if(vec.at(2) != "*") {
						sscanf(vec.at(2).c_str(), "%lf", &value);
						if(formatType == 2) {
							// vertex number must be in the interval 0 <= i < n
							if(a >= n or b >= n) {
								BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
								cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
							}
							g->addEdge(a, b, value);
							// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
						} else {
							if(a > n or b > n) {
								BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
								cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
							}
							g->addEdge(a - 1, b - 1, value);
							// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
						}
					} else {  // special notation for directed edges (add 1 to imbalance)
						imbalance++;
					}
				} catch( boost::bad_lexical_cast const& ) {
					BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
					BOOST_LOG_TRIVIAL(error) << "vec.at(0) = " << vec.at(0);
					BOOST_LOG_TRIVIAL(error) << "vec.at(1) = " << vec.at(1);
					BOOST_LOG_TRIVIAL(error) << "vec.at(2) = " << vec.at(2);
				}

			}
		} else if(formatType == 0) {  // xpress files
			char_separator<char> sep3(" (),\t[]");
			while(std::getline(infile, line)) {
				if(line.find("Mrel:") != string::npos) {
					line = line.substr(line.find("Mrel:") + 5);
				}
				trim(line);
				if(line.length() == 0) {
					continue;
				}
				tokenizer< char_separator<char> > tokens2(line, sep3);
				vector<string> vec2;
				vec2.assign(tokens2.begin(),tokens2.end());

				long size = vec2.size();
				if (size % 3 != 0) {
					BOOST_LOG_TRIVIAL(fatal) << "Error: invalid XPRESS file format!" << std::endl;
				}
				for(long i = 0; i + 2 < size; i = i + 3) {
					try {
						// std::cout << "Processing line " << vec2.at(i) << " " << vec2.at(i+1) << " " << vec2.at(i+2) << std::endl;
						long a = boost::lexical_cast<long>(vec2.at(i));
						long b = boost::lexical_cast<long>(vec2.at(i + 1));
						double value = 0.0;
						sscanf(vec2.at(i + 2).c_str(), "%lf", &value);
						// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;

						if(a > n or b > n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
						}
						g->addEdge(a - 1, b - 1, value);
						// g->addEdge(b - 1, a - 1, value);
					} catch( boost::bad_lexical_cast const& ) {
						BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
						BOOST_LOG_TRIVIAL(error) << "vec.at(0) = " << vec2.at(i);
						BOOST_LOG_TRIVIAL(error) << "vec.at(1) = " << vec2.at(1+1);
						BOOST_LOG_TRIVIAL(error) << "vec.at(2) = " << vec2.at(i+2);
					}
				}
			}
		} else if(formatType == 3) {  // formatType == 3, .dat files
			char_separator<char> sep3(" \t");
			unsigned long a = 0;
			while (std::getline(infile, line)) {
				trim(line);
				tokenizer< char_separator<char> > tokens2(line, sep3);
				vector<string> vec;
				vec.assign(tokens2.begin(),tokens2.end());
				// cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

				for(unsigned long b = 0; b < vec.size(); b++) {
					try {
						double value = 0.0;
						sscanf(vec.at(b).c_str(), "%lf", &value);
						// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
						// the following is to avoid duplicate couting of arcs in the objective function
						if(a >= n or b >= n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
						}
						g->addEdge(a, b, value);
					} catch( boost::bad_lexical_cast const& ) {
						BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
						BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
						BOOST_LOG_TRIVIAL(error) << "vec.at(b) = " << vec.at(b);
					}
				}
				a++;
			}
		} else {  // formatType == 4, matrix market files
			while (std::getline(infile, line)) {
				trim(line);
				tokenizer< char_separator<char> > tokens2(line, sep2);
				vector<string> vec;
				vec.assign(tokens2.begin(),tokens2.end());

				if (vec.size() < 2) continue;
				try {
					long a = boost::lexical_cast<long>(vec.at(0));
					long b = boost::lexical_cast<long>(vec.at(1));
					double value = 1.0;
					if(a > n or b > n) {
						BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
						cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
					}
					g->addEdge(a - 1, b - 1, value);
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
				} catch( boost::bad_lexical_cast const& ) {
					BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
					BOOST_LOG_TRIVIAL(error) << "vec.at(0) = " << vec.at(0);
					BOOST_LOG_TRIVIAL(error) << "vec.at(1) = " << vec.at(1);
					BOOST_LOG_TRIVIAL(error) << "vec.at(2) = " << vec.at(2);
				}

			}
		}
		infile.close();
		BOOST_LOG_TRIVIAL(info) << "Successfully read graph file.";
	} else {
		BOOST_LOG_TRIVIAL(error) << "Failed to read graph file.";
	}
	g->setGlobalN(n);
	BOOST_LOG_TRIVIAL(trace) << "Successfully created local signed graph with " << num_vertices(*(g->graph)) << " vertices.";
	return g;
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromString(const string& graphContents) {
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	long n = 0, e = 0;
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

		if(line.find("%%MatrixMarket") != string::npos) {  // matrix market files
			BOOST_LOG_TRIVIAL(trace) << "Format type is 4 (matrix market)" << endl;
			string line = lines.at(0);
			trim(line);
			lines.erase(lines.begin());
			// captura as dimensoes da matriz e o numero de arestas
			tokenizer< char_separator<char> > tokens2(line, sep2);
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			n = boost::lexical_cast<long>(vec.at(0));
			formatType = 4;
		} else if(line.find("people") != string::npos) {  // xpress files
			BOOST_LOG_TRIVIAL(trace) << "Format type is 0" << endl;
			string firstLine = lines.at(0);
			string number = line.substr(line.find("people:") + 7);
			trim(number);
			n = boost::lexical_cast<long>(number);
			BOOST_LOG_TRIVIAL(trace) << "n value is " << n << endl;
			formatType = 0;
		} else if(line.find("Vertices") != string::npos) {
			BOOST_LOG_TRIVIAL(trace) << "Format type is 1" << endl;
			string number = line.substr(line.find("Vertices") + 8);
			trim(number);
			n = boost::lexical_cast<long>(number);
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
			n = boost::lexical_cast<long>(vec.at(0));
			if(vec.size() == 1) { // .dat files
				BOOST_LOG_TRIVIAL(trace) << "Format type is 3" << endl;
				formatType = 3;
			} else {
				BOOST_LOG_TRIVIAL(trace) << "Format type is 2" << endl;
				// e = boost::lexical_cast<long>(vec.at(1));
				formatType = 2;
			}
		}
	} catch( boost::bad_lexical_cast const& e ) {
	    std::cerr << "Error: input string was not valid" << std::endl;
	    cerr << e.what() << "\n";
	    BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
	}

	SignedGraphPtr g = boost::make_shared<SignedGraph>(n, this->graph);
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
				long a = boost::lexical_cast<long>(vec.at(0));
				long b = boost::lexical_cast<long>(vec.at(1));
				double value = 0.0;
				if(vec.at(2) != "*") {
					sscanf(vec.at(2).c_str(), "%lf", &value);
					if(formatType == 2) {
						// vertex number must be in the interval 0 <= i < n
						if(a >= n or b >= n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
						}
						g->addEdge(a, b, value);
						// std::cout << "Adding edge (" << a << ", " << b << ") = " << value << std::endl;
					} else {
						if(a > n or b > n) {
							BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
							cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
						}
						g->addEdge(a - 1, b - 1, value);
						// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					}
				} else {  // special notation for directed edges (add 1 to imbalance)
					imbalance++;
				}
			} catch( boost::bad_lexical_cast const& ) {
				BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
				BOOST_LOG_TRIVIAL(error) << "vec.at(0) = " << vec.at(0);
				BOOST_LOG_TRIVIAL(error) << "vec.at(1) = " << vec.at(1);
				BOOST_LOG_TRIVIAL(error) << "vec.at(2) = " << vec.at(2);
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

		long size = vec2.size();
		if (size % 3 != 0) {
			BOOST_LOG_TRIVIAL(fatal) << "Error: invalid XPRESS file format!" << std::endl;
		}
		for(long i = 0; i + 2 < size; i = i + 3) {
			try {
				// std::cout << "Processing line " << vec2.at(i) << " " << vec2.at(i+1) << " " << vec2.at(i+2) << std::endl;
				long a = boost::lexical_cast<long>(vec2.at(i));
				long b = boost::lexical_cast<long>(vec2.at(i + 1));
				double value = 0.0;
				sscanf(vec2.at(i + 2).c_str(), "%lf", &value);
				// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;

				if(a > n or b > n) {
					BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
					cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
				}
				g->addEdge(a - 1, b - 1, value);
				// g->addEdge(b - 1, a - 1, value);
			} catch( boost::bad_lexical_cast const& ) {
				BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
				BOOST_LOG_TRIVIAL(error) << "vec.at(0) = " << vec2.at(i);
				BOOST_LOG_TRIVIAL(error) << "vec.at(1) = " << vec2.at(1+1);
				BOOST_LOG_TRIVIAL(error) << "vec.at(2) = " << vec2.at(i+2);
			}
		}
	} else if(formatType == 3) {  // formatType == 3, .dat files
		char_separator<char> sep3(" \t");
		unsigned long a = 0;
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep3);
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());
			// cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

			for(unsigned long b = 0; b < vec.size(); b++) {
				try {
					double value = 0.0;
	                sscanf(vec.at(b).c_str(), "%lf", &value);
					// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
					// the following is to avoid duplicate couting of arcs in the objective function
	                if(a >= n or b >= n) {
						BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
						cerr << "Error: invalid edge. Vertex number must be less than n (" << n << ").";
					}
					g->addEdge(a, b, value);
				} catch( boost::bad_lexical_cast const& ) {
					BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
				}
			}
			a++;
		}
	} else {  // formatType == 4, matrix market files
		long imbalance = 0;
		while (not lines.empty()) {
			string line = lines.back();
			trim(line);
			lines.pop_back();
			tokenizer< char_separator<char> > tokens2(line, sep2);
			vector<string> vec;
			vec.assign(tokens2.begin(),tokens2.end());

			if (vec.size() < 2) continue;
			try {
				long a = boost::lexical_cast<long>(vec.at(0));
				long b = boost::lexical_cast<long>(vec.at(1));
				double value = 1.0;
				if(a > n or b > n) {
					BOOST_LOG_TRIVIAL(error) << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
					cerr << "Error: invalid edge. Vertex number must be less or equal to n (" << n << ").";
				}
				g->addEdge(a - 1, b - 1, value);
				// std::cout << "Adding edge (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
			} catch( boost::bad_lexical_cast const& ) {
				BOOST_LOG_TRIVIAL(fatal) << "Error: input string was not valid" << std::endl;
			}
		}
	}
	// g->printGraph();
	BOOST_LOG_TRIVIAL(info) << "Successfully read graph file.";

	return g;
}

SignedGraphPtr SimpleTextGraphFileReader::readGraphFromFile(const string& filepath, const bool& parallelgraph) {
	BOOST_LOG_TRIVIAL(info) << "Reading input file: '" << filepath << "' ..." << endl;
	SignedGraphPtr g;
	if(not parallelgraph) {
		g = readGraphFromString(get_file_contents(filepath.c_str()));
	} else {
		g = readGraphFromFilepath(filepath, parallelgraph);
		// All processes synchronize at this point, then the graph is complete
		BOOST_LOG_TRIVIAL(info) << "Synchronizing process for global graph creation...";
		synchronize(g->graph->process_group());  // this is the synchronization from master process (id 0)
	}
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

bool SimpleTextGraphFileReader::exportGraphToGraphVizFile(SignedGraph &g, const std::string outputFolder,
		const std::string &filename) {
	namespace fs = boost::filesystem;
	if (!fs::exists(fs::path(outputFolder))) {
		fs::create_directories(fs::path(outputFolder));
	}
	stringstream sfilename;
	sfilename << outputFolder << "/" << filename << ".dot";
	fs::path newFile(sfilename.str());
	ofstream os;
	os.open(newFile.c_str(), ios::out | ios::trunc);
	if (!os) {
		BOOST_LOG_TRIVIAL(fatal) << "Can't open output file!" << endl;
		throw "Cannot open output file.";
	}

	long n = g.getN();
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g.graph));
	ParallelGraph::edge_descriptor e;

	os << "graph graphname {\n";
	// For each vertex i
	for(long i = 0; i < n; i++) {
		ParallelGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(vertex(i, *(g.graph)), *(g.graph)); f != l; ++f) {
			e = *f;
			double weight = ew[e].weight;
			long j = target(*f, *(g.graph)).local;
			if(i < j) {
				os << i << " -- " << j << "[ weight = " << weight << " ];\n";
			}
		}
	}
	os << "}\n";

	os.close();
	return true;
}

} /* namespace controller */
