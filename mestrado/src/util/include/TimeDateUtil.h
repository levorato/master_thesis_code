/*
 * TimeDateUtil.h
 *
 *  Created on: 23/05/2013
 *      Author: czt0
 */

#ifndef TIMEDATEUTIL_H_
#define TIMEDATEUTIL_H_

#include <boost/date_time/posix_time/posix_time.hpp>

namespace util {

class TimeDateUtil {
public:
	TimeDateUtil();
	virtual ~TimeDateUtil();

	static std::wstring FormatTime(boost::posix_time::ptime now);
};

} /* namespace util */
#endif /* TIMEDATEUTIL_H_ */
