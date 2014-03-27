/*
 * MPIMessage.cpp
 *
 *  Created on: May 30, 2013
 *      Author: mario
 */

#include "include/MPIMessage.h"

namespace util {

MPIMessage::MPIMessage() {
	// TODO Auto-generated constructor stub

}

MPIMessage::~MPIMessage() {
	// TODO Auto-generated destructor stub
}

} /* namespace util */

// important serialization polymorphic class export
//BOOST_SERIALIZATION_ASSUME_ABSTRACT( util::InputMessage )
//BOOST_CLASS_EXPORT_GUID( util::InputMessageParallelGrasp, "InputMessageParallelGrasp" )
//BOOST_CLASS_EXPORT_GUID( util::InputMessageParallelVNS, "InputMessageParallelVNS" )
