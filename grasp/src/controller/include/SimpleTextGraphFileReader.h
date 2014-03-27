/*
 * SimpleTextGraphFileReader.h
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#ifndef SIMPLETEXTGRAPHFILEREADER_H_
#define SIMPLETEXTGRAPHFILEREADER_H_

#include "../../graph/include/Graph.h"
#include <string>

using namespace clusteringgraph;

namespace controller {

class SimpleTextGraphFileReader {
public:
	SimpleTextGraphFileReader();
	virtual ~SimpleTextGraphFileReader();

	SignedGraphPtr readGraphFromFile(const std::string &filepath);

	SignedGraphPtr readGraphFromString(const std::string &graphContents);

private:
	std::string get_file_contents(const char *filename);
};

} /* namespace controller */
#endif /* SIMPLETEXTGRAPHFILEREADER_H_ */
