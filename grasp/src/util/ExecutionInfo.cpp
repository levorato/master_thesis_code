/*
 * ExecutionInfo.cpp
 *
 *  Created on: 24/04/2014
 *      Author: czt0
 */

#include "include/ExecutionInfo.h"

namespace util {

using namespace std;

ExecutionInfo::ExecutionInfo(string id, string fileid, string outfolder, int rank) :
	executionId(id), fileId(fileid), outputFolder(outfolder), processRank(rank) {
	// TODO Auto-generated constructor stub

}

ExecutionInfo::~ExecutionInfo() {
	// TODO Auto-generated destructor stub
}

} /* namespace util */
