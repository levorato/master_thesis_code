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

	SignedGraphPtr readGraphFromFile(std::string filepath);
};

} /* namespace controller */
#endif /* SIMPLETEXTGRAPHFILEREADER_H_ */
