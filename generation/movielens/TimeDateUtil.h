/*
 * TimeDateUtil.h
 *
 *  Created on: 23/05/2013
 *      Author: Mario Levorato
 */

#ifndef TIMEDATEUTIL_H_
#define TIMEDATEUTIL_H_

#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>
#include <map>
#include <stdexcept>
#include <iostream>

#include <boost/log/trivial.hpp>

namespace util {

using namespace std;
namespace logging = boost::log;

class TimeDateUtil {
public:
	TimeDateUtil();
	virtual ~TimeDateUtil();

	static std::string FormatTime(boost::posix_time::ptime now);
	static std::string getDateAndTime();
};

class EnumUtil {
public:
	EnumUtil();
	virtual ~EnumUtil();
};

} /* namespace util */
#endif /* TIMEDATEUTIL_H_ */
