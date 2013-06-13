/*
 * MPIMessage.h
 *
 *  Created on: May 30, 2013
 *      Author: mario
 */

#ifndef MPIMESSAGE_H_
#define MPIMESSAGE_H_

#include <string>
#include <sstream>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "../../graph/include/Clustering.h"

using namespace std;
using namespace clusteringgraph;

namespace util {

class InputMessage {
public:
	static const int TAG = 50;
	string graphInputFileContents;
	double alpha;
	int l;
	int iter;
	int problemType;
	string fileId;
	long timeLimit;

	InputMessage() : graphInputFileContents(),
			alpha(0.0F), l(1), iter(500), problemType(0), fileId("noId"), timeLimit(1800) {

	}

	InputMessage(string graphContents, int it, double a, int neigh,
			int pType, string id, long t) : graphInputFileContents(graphContents),
					alpha(a), l(neigh), iter(it), problemType(pType), fileId(id),
					timeLimit(t) {

	}

	string toString() {
		stringstream ss;
		ss << "Alpha: " << alpha << "; l = " << l << "; iter = " << iter << "; fileId = " <<
				fileId << "; " << graphInputFileContents << "\n\n";
		return ss.str();
	}

private:
	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & graphInputFileContents;
		ar & alpha;
		ar & l;
		ar & iter;
		ar & problemType;
		ar & fileId;
		ar & timeLimit;
	}
};

class OutputMessage {
public:
	static const int TAG = 60;
	Clustering clustering;

	OutputMessage() : clustering() {

	}

	OutputMessage(Clustering &c) : clustering(c) {

	}

private:
	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clustering;
	}
};

class MPIMessage {
public:
	MPIMessage();
	virtual ~MPIMessage();
};

} /* namespace util */
#endif /* MPIMESSAGE_H_ */
