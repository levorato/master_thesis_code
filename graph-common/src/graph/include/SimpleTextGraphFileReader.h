/*
 * SimpleTextGraphFileReader.h
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#ifndef SIMPLETEXTGRAPHFILEREADER_H_
#define SIMPLETEXTGRAPHFILEREADER_H_

#include "Graph.h"
#include "ParallelBGLSignedGraph.h"
#include <string>

using namespace clusteringgraph;

namespace controller {

class SimpleTextGraphFileReader {
public:
	SimpleTextGraphFileReader(clusteringgraph::ParallelGraph *pgraph) :
		graph(pgraph) {
	}
	virtual ~SimpleTextGraphFileReader();

	SignedGraphPtr readGraphFromFile(const std::string &filepath, const bool& parallelgraph);

	SignedGraphPtr readGraphFromString(const std::string &graphContents);

	bool exportGraphToGraphVizFile(SignedGraph &g, const std::string outputFolder,
			const std::string &filename);

protected:
	std::string get_file_contents(const char *filename);
	SignedGraphPtr readGraphFromFilepath(const string& filepath, const bool& parallelgraph);

	clusteringgraph::ParallelGraph *graph;
};

} /* namespace controller */
#endif /* SIMPLETEXTGRAPHFILEREADER_H_ */
