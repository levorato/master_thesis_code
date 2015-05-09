/*
 * ProcessUtil.h
 *
 *  Created on: May 8, 2015
 *      Author: mlevorato
 */

#ifndef SRC_UTIL_INCLUDE_PROCESSUTIL_H_
#define SRC_UTIL_INCLUDE_PROCESSUTIL_H_

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <boost/log/trivial.hpp>

namespace util {

class ProcessUtil {
public:
	ProcessUtil();
	virtual ~ProcessUtil();

	static int exec(char* cmd) {
		BOOST_LOG_TRIVIAL(info) << "Command-line call: '" << cmd << "'";
		int errorCount = 0, successCount = 0;

	    FILE* pipe = popen(cmd, "r");
	    if (!pipe) { errorCount++; }
	    char buffer[128];
	    std::string result = "";
	    while(!feof(pipe)) {
	    	if(fgets(buffer, 128, pipe) != NULL)
	    		result += buffer;
	    }
	    pclose(pipe);
	    if(result.find("Graph Information:") != std::string::npos) {
	    	successCount++;
		}
	    BOOST_LOG_TRIVIAL(info) << "stdout: " << result;
	    return ((errorCount > 0) or (successCount <= 0));
	}
};

} /* namespace util */

#endif /* SRC_UTIL_INCLUDE_PROCESSUTIL_H_ */
