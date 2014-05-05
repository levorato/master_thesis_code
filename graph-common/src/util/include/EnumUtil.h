/*
 * EnumUtil.h
 *
 *  Created on: 11/07/2013
 *      Author: mlevorato
 */

#ifndef ENUMUTIL_H_
#define ENUMUTIL_H_

#include <string>
#include <map>
#include <stdexcept>

#include <boost/log/trivial.hpp>

using namespace std;
namespace logging = boost::log;

namespace util {

class EnumUtil {
public:
	EnumUtil();
	virtual ~EnumUtil();
};

template <typename T>
class EnumParser
{
    map <string, T> enumMap;
public:
    EnumParser(){
    	enumMap["fatal"] = logging::trivial::fatal;
    	enumMap["error"] = logging::trivial::error;
    	enumMap["warning"] = logging::trivial::warning;
    	enumMap["info"] = logging::trivial::info;
    	enumMap["debug"] = logging::trivial::debug;
    	enumMap["trace"] = logging::trivial::trace;
    };

    T ParseSomeEnum(const string &value)
    {
    	typename map <string, T>::const_iterator iValue = enumMap.find(value);
        if (iValue  == enumMap.end()) {
        	cout << "enum not found\n";
            throw runtime_error("item de enum nao encontrado");
        }
        return iValue->second;
    }
};

typedef	EnumParser<logging::trivial::severity_level> LogSeverityEnumParser;

} // namespace util

#endif /* ENUMUTIL_H_ */
