/*
 * EnumUtil.cpp
 *
 *  Created on: 11/07/2013
 *      Author: mlevorato
 */

#include "include/EnumUtil.h"

namespace util {
EnumUtil::EnumUtil() {
	// TODO Auto-generated constructor stub

}

EnumUtil::~EnumUtil() {
	// TODO Auto-generated destructor stub
}

template<> LogSeverityEnumParser::EnumParser() {
	enumMap["fatal"] = logging::trivial::fatal;
	enumMap["error"] = logging::trivial::error;
	enumMap["warning"] = logging::trivial::warning;
	enumMap["info"] = logging::trivial::info;
	enumMap["debug"] = logging::trivial::debug;
	enumMap["trace"] = logging::trivial::trace;
}

}
