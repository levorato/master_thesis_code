/*
 * SimpleTextGraphFileReader.h
 *
 *  Created on: May 1, 2013
 *      Author: mario
 */

#ifndef SIMPLETEXTGRAPHFILEREADER_VALIDATION_H_
#define SIMPLETEXTGRAPHFILEREADER_VALIDATION_H_

#include "./Graph.h"
#include <string>



namespace clusteringgraph {
namespace validation {

class SimpleTextGraphFileReader {
public:
	SimpleTextGraphFileReader();
	virtual ~SimpleTextGraphFileReader();

	SignedGraphPtr readGraphFromFile(const std::string &filepath);

	SignedGraphPtr readGraphFromString(const std::string &graphContents);

	bool exportGraphToGraphVizFile(SignedGraph &g, const std::string outputFolder,
			const std::string &filename);

private:
	std::string get_file_contents(const char *filename);
};

}

} /* namespace controller */
#endif /* SIMPLETEXTGRAPHFILEREADER_H_ */
