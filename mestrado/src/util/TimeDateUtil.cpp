/*
 * TimeDateUtil.cpp
 *
 *  Created on: 23/05/2013
 *      Author: czt0
 */

#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>
#include "include/TimeDateUtil.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace util {

TimeDateUtil::TimeDateUtil() {
	// TODO Auto-generated constructor stub

}

TimeDateUtil::~TimeDateUtil() {
	// TODO Auto-generated destructor stub
}

std::wstring TimeDateUtil::FormatTime(boost::posix_time::ptime now)
{
  using namespace boost::posix_time;
  static std::locale loc(std::wcout.getloc(),
                         new wtime_facet(L"%Y%m%d_%H%M%S"));

  std::basic_stringstream<wchar_t> wss;
  wss.imbue(loc);
  wss << now;

  boost::uuids::uuid tag = boost::uuids::random_generator()();
  wss << "_" << tag;

  return wss.str();
}

std::string TimeDateUtil::getTimeAndDateAsString() {
	using namespace boost::posix_time;
	ptime now = second_clock::local_time();

	std::wstring ws(TimeDateUtil::FormatTime(now));
	return std::string ( ws.begin(), ws.end() );;
}

} /* namespace util */
