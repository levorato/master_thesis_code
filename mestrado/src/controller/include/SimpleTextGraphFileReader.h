/*
 * SimpleTextGraphFileReader.h
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#ifndef SIMPLETEXTGRAPHFILEREADER_H_
#define SIMPLETEXTGRAPHFILEREADER_H_

#include "../../graph/GraphTypes.h"
#include <string>

namespace controller {

class SimpleTextGraphFileReader {
public:
	SimpleTextGraphFileReader();
	virtual ~SimpleTextGraphFileReader();

	static *SignedGraph readGraphFromFile(std::string filepath);
};

} /* namespace controller */
#endif /* SIMPLETEXTGRAPHFILEREADER_H_ */
