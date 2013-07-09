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

};

class InputMessageParallelGrasp : InputMessage {
public:
	static const int TAG = 50;
	string graphInputFileContents;
	double alpha;
	int l;
	int iter;
	int gainFunctionType;
	int problemType;
	string fileId;
	string outputFolder;
	long timeLimit;

	InputMessageParallelGrasp() : graphInputFileContents(),
			alpha(0.0F), l(1), iter(500), gainFunctionType(0), problemType(0),
			fileId("noId"), outputFolder(""), timeLimit(1800) {

	}

	InputMessageParallelGrasp(string graphContents, int it, double a, int neigh,
			int pType, int gfType, string id, string folder, long t) :
				graphInputFileContents(graphContents),
					alpha(a), l(neigh), iter(it), gainFunctionType(gfType),
					problemType(pType), fileId(id),
					outputFolder(folder), timeLimit(t) {

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
		ar & gainFunctionType;
		ar & problemType;
		ar & fileId;
		ar & outputFolder;
		ar & timeLimit;
	}
};

class InputMessageParallelVNS : InputMessage {
public:
	static const int TAG = 70;
	int l;
	string graphInputFileContents;
	Clustering clustering;
	int problemType;
	double timeSpentSoFar;
	double timeLimit;
	unsigned long initialClusterIndex;
	unsigned long finalClusterIndex;
	int numberOfSlaves;
	int numberOfSearchSlaves;

	InputMessageParallelVNS() : l(1), graphInputFileContents(), clustering(),
			problemType(0), timeSpentSoFar(0.0), timeLimit(3600.0), initialClusterIndex(0),
			finalClusterIndex(0), numberOfSlaves(0), numberOfSearchSlaves(0) {

	}

	InputMessageParallelVNS(int neig, string graphContents, Clustering c,
			int pType, double timeSoFar, double tl, unsigned long startIdx,
			unsigned long endIdx, int slaves, int searchSlaves) : l(neig),
			graphInputFileContents(graphContents), clustering(c),
			problemType(pType), timeSpentSoFar(timeSoFar), timeLimit(tl),
			initialClusterIndex(startIdx) , finalClusterIndex(endIdx),
			numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves) {

	}

private:
	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & l;
		ar & graphInputFileContents;
		ar & clustering;
		ar & problemType;
		ar & timeSpentSoFar;
		ar & timeLimit;
		ar & initialClusterIndex;
		ar & finalClusterIndex;
		ar & numberOfSlaves;
		ar & numberOfSearchSlaves;
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
