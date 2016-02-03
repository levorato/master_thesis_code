/*
 * TimeDateUtil.cpp
 *
 *  Created on: 23/05/2013
 *      Author: Mario Levorato
 */

#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>
#include <boost/lexical_cast.hpp>

#include "TimeDateUtil.h"

namespace util {

TimeDateUtil::TimeDateUtil() {
	// TODO Auto-generated constructor stub

}

TimeDateUtil::~TimeDateUtil() {
	// TODO Auto-generated destructor stub
}

std::string TimeDateUtil::FormatTime(boost::posix_time::ptime now)
{
  using namespace boost::posix_time;
  static std::locale loc(std::wcout.getloc(),
                         new wtime_facet(L"%Y%m%d_%H%M%S"));

  std::basic_stringstream<wchar_t> wss;
  wss.imbue(loc);
  wss << now;

  std::wstring ws = wss.str();

  return std::string ( ws.begin(), ws.end() );
}

std::string TimeDateUtil::getDateAndTime() {
	using namespace boost::posix_time;
	ptime now = second_clock::local_time();

	return TimeDateUtil::FormatTime(now);
}


} /* namespace util */
