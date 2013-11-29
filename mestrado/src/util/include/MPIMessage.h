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
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include "../../graph/include/Clustering.h"

using namespace std;
using namespace clusteringgraph;
namespace bser = boost::serialization;

namespace util {

/**
 * Abstract input message class for use with MPI.
 */
class InputMessage {
public:
	// the identifier of the graph
	unsigned int id;
	string graphInputFileContents;
	int l;
	// the number of grasp slave processes
	unsigned int numberOfSlaves;
	// the number of vns slave processes
	unsigned int numberOfSearchSlaves;

	InputMessage() : id(0), graphInputFileContents(), l(0), numberOfSlaves(0),
			numberOfSearchSlaves(0) {

	}

	InputMessage(unsigned int i, string graphContents, int nl, unsigned int slaves, unsigned int searchSlaves) :
		id(i), graphInputFileContents(graphContents), l(nl), numberOfSlaves(slaves),
		numberOfSearchSlaves(searchSlaves) {

	}

	virtual ~InputMessage(){};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & id;
		ar & graphInputFileContents;
		ar & l;
		ar & numberOfSlaves;
		ar & numberOfSearchSlaves;
	}
};


class InputMessageParallelGrasp : public InputMessage {
public:
	static const int TAG = 50;
	double alpha;
	int iter;
	int gainFunctionType;
	int problemType;
	string executionId;
	string fileId;
	string outputFolder;
	long timeLimit;
	bool firstImprovementOnOneNeig;

	InputMessageParallelGrasp() : InputMessage(),
			alpha(0.0F), iter(500), gainFunctionType(0), problemType(0),
			fileId("noId"), outputFolder(""), timeLimit(1800), firstImprovementOnOneNeig(false) {

	}

	InputMessageParallelGrasp(unsigned int i, string graphContents, int it, double a, int neigh,
			int pType, int gfType, string eid, string fid, string folder, long t, unsigned int slaves,
			unsigned int searchSlaves, bool fiOneNeig) :
				InputMessage(i, graphContents, neigh, slaves, searchSlaves),
					alpha(a), iter(it), gainFunctionType(gfType),
					problemType(pType), executionId(eid), fileId(fid),
					outputFolder(folder), timeLimit(t), firstImprovementOnOneNeig(fiOneNeig) {

	}

	string toString() {
		stringstream ss;
		ss << "Alpha: " << alpha << "; l = " << l << "; iter = " << iter << "; fileId = " <<
				fileId << "; " << graphInputFileContents << "\n";
		return ss.str();
	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar.template register_type< InputMessage >();
		// ar.template register_type< InputMessageParallelGrasp >();
		//ar.template register_type< InputMessage >();
		// ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InputMessage);
		// invoke serialization of the base class
		ar & boost::serialization::base_object<InputMessage>(*this);
		boost::serialization::void_cast_register<InputMessageParallelGrasp, InputMessage>(
			static_cast<InputMessageParallelGrasp *>(NULL),
			static_cast<InputMessage *>(NULL)
		);
		// save/load class member variables
		ar & alpha;
		ar & iter;
		ar & gainFunctionType;
		ar & problemType;
		ar & executionId;
		ar & fileId;
		ar & outputFolder;
		ar & timeLimit;
		ar & firstImprovementOnOneNeig;
	}
};

class InputMessageParallelVNS : public InputMessage {
public:
	static const int TAG = 70;
	Clustering clustering;
	int problemType;
	double timeSpentSoFar;
	double timeLimit;
	unsigned long initialClusterIndex;
	unsigned long finalClusterIndex;
	unsigned long k; /* number of max clusters (RCC Problem only) */

	InputMessageParallelVNS() : InputMessage(), clustering(),
			problemType(0), timeSpentSoFar(0.0), timeLimit(3600.0), initialClusterIndex(0),
			finalClusterIndex(0), k(0) {

	}

	InputMessageParallelVNS(unsigned int i, int neig, string graphContents, Clustering c,
			int pType, double timeSoFar, double tl, unsigned long startIdx,
			unsigned long endIdx, unsigned int slaves, unsigned int searchSlaves, unsigned long _k) :
			InputMessage(i, graphContents, neig, slaves, searchSlaves),
			clustering(c),
			problemType(pType), timeSpentSoFar(timeSoFar), timeLimit(tl),
			initialClusterIndex(startIdx) , finalClusterIndex(endIdx), k(_k) {

	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar.template register_type< InputMessage >();
		// ar.template register_type< InputMessageParallelVNS >();
		// ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InputMessage);
		// invoke serialization of the base class
		ar & boost::serialization::base_object<InputMessage>(*this);
		boost::serialization::void_cast_register<InputMessageParallelVNS, InputMessage>(
			static_cast<InputMessageParallelVNS *>(NULL),
			static_cast<InputMessage *>(NULL)
		);
		// save/load class member variables
		ar & clustering;
		ar & problemType;
		ar & timeSpentSoFar;
		ar & timeLimit;
		ar & initialClusterIndex;
		ar & finalClusterIndex;
		ar & k;
	}
};



class OutputMessage {
public:
	static const int TAG = 60;
	Clustering clustering;
	long numberOfTestedCombinations;

	OutputMessage() : clustering(), numberOfTestedCombinations(0) {

	}

	OutputMessage(Clustering &c, long nc) : clustering(c), numberOfTestedCombinations(nc) {

	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clustering;
		ar & numberOfTestedCombinations;
	}
};

class MPIMessage {
public:
	// Message with the number of grasp slaves
	static const int INPUT_MSG_NUM_SLAVES_TAG = 40;
	// Parallel Grasp message tag
	static const int INPUT_MSG_PARALLEL_GRASP_TAG = 50;
	static const int OUTPUT_MSG_PARALLEL_GRASP_TAG = 60;
	// Parallel VNS message tag
	static const int INPUT_MSG_PARALLEL_VNS_TAG = 70;
	static const int OUTPUT_MSG_PARALLEL_VNS_TAG = 80;
	static const int INTERRUPT_MSG_PARALLEL_VNS_TAG = 85;
	// Other tags
	static const int TERMINATE_MSG_TAG = 90;
	static const int LEADER_ID = 0;

	MPIMessage();
	virtual ~MPIMessage();
};

} /* namespace util */

#endif /* MPIMESSAGE_H_ */
