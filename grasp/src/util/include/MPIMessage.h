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
	string graphInputFilePath;
	int l;
	// the number of master processes
	unsigned int numberOfMasters;
	// the number of VND slave processes
	unsigned int numberOfSearchSlaves;

	InputMessage() : id(0), graphInputFilePath(), l(0), numberOfMasters(0),
			numberOfSearchSlaves(0) {

	}

	InputMessage(unsigned int i, string graphFilePath, int nl, unsigned int masters, unsigned int searchSlaves) :
		id(i), graphInputFilePath(graphFilePath), l(nl), numberOfMasters(masters),
		numberOfSearchSlaves(searchSlaves) {

	}

	virtual ~InputMessage(){};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & id;
		ar & graphInputFilePath;
		ar & l;
		ar & numberOfMasters;
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

	InputMessageParallelGrasp(unsigned int i, string graphFilePath, int it, double a, int neigh,
			int pType, int gfType, string eid, string fid, string folder, long t, unsigned int masters,
			unsigned int searchSlaves, bool fiOneNeig) :
				InputMessage(i, graphFilePath, neigh, masters, searchSlaves),
					alpha(a), iter(it), gainFunctionType(gfType),
					problemType(pType), executionId(eid), fileId(fid),
					outputFolder(folder), timeLimit(t), firstImprovementOnOneNeig(fiOneNeig) {

	}

	string toString() {
		stringstream ss;
		ss << "Alpha: " << alpha << "; l = " << l << "; iter = " << iter << "; fileId = " <<
				fileId << "; " << graphInputFilePath << "\n";
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

class InputMessageParallelILS : public InputMessage {
public:
	static const int TAG = 60;
	double alpha;
	int iter;
	int gainFunctionType;
	int problemType;
	string executionId;
	string fileId;
	string outputFolder;
	long timeLimit;
	bool firstImprovementOnOneNeig;
	int iterMaxILS;
	int perturbationLevelMax;

	InputMessageParallelILS() : InputMessage(),
			alpha(0.0F), iter(400), gainFunctionType(0), problemType(0),
			fileId("noId"), outputFolder(""), timeLimit(1800), firstImprovementOnOneNeig(false),
			iterMaxILS(3), perturbationLevelMax(7) {

	}

	InputMessageParallelILS(unsigned int i, string graphFilePath, int it, double a, int neigh,
			int pType, int gfType, string eid, string fid, string folder, long t, unsigned int masters,
			unsigned int searchSlaves, bool fiOneNeig, int maxilsiter, int maxpertlevel) :
				InputMessage(i, graphFilePath, neigh, masters, searchSlaves),
					alpha(a), iter(it), gainFunctionType(gfType),
					problemType(pType), executionId(eid), fileId(fid),
					outputFolder(folder), timeLimit(t), firstImprovementOnOneNeig(fiOneNeig),
					iterMaxILS(maxilsiter), perturbationLevelMax(maxpertlevel) {

	}

	string toString() {
		stringstream ss;
		ss << "Alpha: " << alpha << "; l = " << l << "; iter = " << iter << "; fileId = " <<
				fileId << "; " << graphInputFilePath << "\n";
		return ss.str();
	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar.template register_type< InputMessage >();
		// ar.template register_type< InputMessageParallelILS >();
		//ar.template register_type< InputMessage >();
		// ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InputMessage);
		// invoke serialization of the base class
		ar & boost::serialization::base_object<InputMessage>(*this);
		boost::serialization::void_cast_register<InputMessageParallelILS, InputMessage>(
			static_cast<InputMessageParallelILS *>(NULL),
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
		ar & iterMaxILS;
		ar & perturbationLevelMax;
	}
};


class InputMessageParallelVND : public InputMessage {
public:
	static const int TAG = 70;
	Clustering clustering;
	int problemType;
	double timeSpentSoFar;
	double timeLimit;
	unsigned long initialClusterIndex;
	unsigned long finalClusterIndex;
	unsigned long k; /* number of max clusters (RCC Problem only) */

	InputMessageParallelVND() : InputMessage(), clustering(),
			problemType(0), timeSpentSoFar(0.0), timeLimit(3600.0), initialClusterIndex(0),
			finalClusterIndex(0), k(0) {

	}

	InputMessageParallelVND(unsigned int i, int neig, string graphFilePath, Clustering c,
			int pType, double timeSoFar, double tl, unsigned long startIdx,
			unsigned long endIdx, unsigned int masters, unsigned int searchSlaves, unsigned long _k) :
			InputMessage(i, graphFilePath, neig, masters, searchSlaves),
			clustering(c),
			problemType(pType), timeSpentSoFar(timeSoFar), timeLimit(tl),
			initialClusterIndex(startIdx) , finalClusterIndex(endIdx), k(_k) {

	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar.template register_type< InputMessage >();
		// ar.template register_type< InputMessageParallelVND >();
		// ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InputMessage);
		// invoke serialization of the base class
		ar & boost::serialization::base_object<InputMessage>(*this);
		boost::serialization::void_cast_register<InputMessageParallelVND, InputMessage>(
			static_cast<InputMessageParallelVND *>(NULL),
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
	// Message with the number of heuristic masters
	static const int INPUT_MSG_NUM_MASTERS_TAG = 40;
	// Parallel Grasp message tag
	static const int INPUT_MSG_PARALLEL_GRASP_TAG = 50;
	static const int OUTPUT_MSG_PARALLEL_GRASP_TAG = 55;
	// Parallel ILS message tag
	static const int INPUT_MSG_PARALLEL_ILS_TAG = 60;
	static const int OUTPUT_MSG_PARALLEL_ILS_TAG = 65;
	// Parallel VND message tag
	static const int INPUT_MSG_PARALLEL_VND_TAG = 70;
	static const int OUTPUT_MSG_PARALLEL_VND_TAG = 80;
	static const int INTERRUPT_MSG_PARALLEL_VND_TAG = 85;
	// Other tags
	static const int TERMINATE_MSG_TAG = 90;
	static const int LEADER_ID = 0;

	MPIMessage();
	virtual ~MPIMessage();
};

} /* namespace util */

#endif /* MPIMESSAGE_H_ */
