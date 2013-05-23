/*
 * TimeDateUtil.cpp
 *
 *  Created on: 23/05/2013
 *      Author: czt0
 */

#include "include/TimeDateUtil.h"

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
  return wss.str();
}

} /* namespace util */
