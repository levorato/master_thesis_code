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
	// if the local search is CUDA-enabled
	bool cudaEnabled;

	InputMessage() : id(0), graphInputFilePath(), l(0), numberOfMasters(0),
			numberOfSearchSlaves(0), cudaEnabled(false) {

	}

	InputMessage(unsigned int i, string graphFilePath, int nl, unsigned int masters, unsigned int searchSlaves, bool cuda) :
		id(i), graphInputFilePath(graphFilePath), l(nl), numberOfMasters(masters),
		numberOfSearchSlaves(searchSlaves), cudaEnabled(cuda) {

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
		ar & cudaEnabled;
	}
};


class InputMessageParallelGrasp : public InputMessage {
public:
	static const int TAG = 50;
	double alpha;
	int iter;
	int gainFunctionType;
	int problemType;
	long k;  // number of clusters of solution (for RCC Problem)
	string executionId;
	string fileId;
	string outputFolder;
	long timeLimit;
	bool firstImprovementOnOneNeig;
	Clustering CCclustering;  // best solution found by CC problem (for use on RCC problem solve)
	std::vector<long> vertexList;  // list of vertices of the subgraph (used on split graph feature)

	InputMessageParallelGrasp() : InputMessage(),
			alpha(0.0F), iter(500), gainFunctionType(0), problemType(0), k(0),
			fileId("noId"), outputFolder(""), timeLimit(1800), firstImprovementOnOneNeig(false),
			CCclustering(), vertexList() {

	}

	InputMessageParallelGrasp(unsigned int i, string graphFilePath, int it, double a, int neigh,
			int pType, int gfType, string eid, string fid, string folder, long t, unsigned int masters,
			unsigned int searchSlaves, bool fiOneNeig, long numberOfClustersInSolution = 0, bool cuda = false, Clustering* cl = NULL) :
				InputMessage(i, graphFilePath, neigh, masters, searchSlaves, cuda),
					alpha(a), iter(it), gainFunctionType(gfType),
					problemType(pType), k(numberOfClustersInSolution), executionId(eid), fileId(fid),
					outputFolder(folder), timeLimit(t), firstImprovementOnOneNeig(fiOneNeig), CCclustering(), vertexList() {
			if(cl != NULL) {
				CCclustering = *cl;
			}
	}

	void setVertexList(std::vector<long> vlist) {
		vertexList = vlist;
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
		ar & k;
		ar & executionId;
		ar & fileId;
		ar & outputFolder;
		ar & timeLimit;
		ar & firstImprovementOnOneNeig;
		ar & CCclustering;
		ar & vertexList;
	}
};

class InputMessageParallelILS : public InputMessage {
public:
	static const int TAG = 60;
	double alpha;
	int iter;
	int gainFunctionType;
	int problemType;
	long k;  // number of clusters of solution (for RCC Problem)
	string executionId;
	string fileId;
	string outputFolder;
	long timeLimit;
	bool firstImprovementOnOneNeig;
	int iterMaxILS;
	int perturbationLevelMax;
	Clustering CCclustering;  // best solution found by CC problem (for use on RCC problem solve)
	bool isSplitGraph;
	std::vector<long> vertexList;  // list of vertices of the subgraph (used on split graph feature)
	bool isParallelGraph;  // indicates that parallel graph from boost parallel bgl is being used
	bool runILS;  // indicates that local ILS must run in the process
	bool redistributeVertices;

	InputMessageParallelILS() : InputMessage(),
			alpha(0.0F), iter(400), gainFunctionType(0), problemType(0), k(0),
			fileId("noId"), outputFolder(""), timeLimit(1800), firstImprovementOnOneNeig(false),
			iterMaxILS(3), perturbationLevelMax(7), CCclustering(), isSplitGraph(true), vertexList(),
			isParallelGraph(true), runILS(true), redistributeVertices(true) {

	}

	InputMessageParallelILS(unsigned int i, string graphFilePath, int it, double a, int neigh,
			int pType, int gfType, string eid, string fid, string folder, long t, unsigned int masters,
			unsigned int searchSlaves, bool fiOneNeig, int maxilsiter, int maxpertlevel,
			long numberOfClustersInSolution = 0, bool cuda = false, Clustering* cl = NULL,
			bool parallelgraph = true, bool runILSproc = true, bool redistVertices = true) :
				InputMessage(i, graphFilePath, neigh, masters, searchSlaves, cuda),
					alpha(a), iter(it), gainFunctionType(gfType),
					problemType(pType), k(numberOfClustersInSolution), executionId(eid), fileId(fid),
					outputFolder(folder), timeLimit(t), firstImprovementOnOneNeig(fiOneNeig),
					iterMaxILS(maxilsiter), perturbationLevelMax(maxpertlevel), CCclustering(),
					isSplitGraph(true), vertexList(), isParallelGraph(parallelgraph), runILS(runILSproc),
					redistributeVertices(redistVertices) {
			if(cl != NULL) {
				CCclustering = *cl;
			}
	}

	void setClustering(Clustering* cl) {
		if(cl != NULL) {
			CCclustering = *cl;
		}
	}

	void setVertexList(std::vector<long> vlist) {
		vertexList = vlist;
		isSplitGraph = true;
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
		ar & k;
		ar & executionId;
		ar & fileId;
		ar & outputFolder;
		ar & timeLimit;
		ar & firstImprovementOnOneNeig;
		ar & iterMaxILS;
		ar & perturbationLevelMax;
		ar & CCclustering;
		ar & isSplitGraph;
		ar & vertexList;
		ar & isParallelGraph;
		ar & runILS;
		ar & redistributeVertices;
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
	long k; /* number of max clusters (RCC Problem only) */

	InputMessageParallelVND() : InputMessage(), clustering(),
			problemType(0), timeSpentSoFar(0.0), timeLimit(3600.0), initialClusterIndex(0),
			finalClusterIndex(0), k(0) {

	}

	InputMessageParallelVND(unsigned int i, int neig, string graphFilePath, Clustering c,
			int pType, double timeSoFar, double tl, unsigned long startIdx,
			unsigned long endIdx, unsigned int masters, unsigned int searchSlaves, long _k, bool cuda) :
			InputMessage(i, graphFilePath, neig, masters, searchSlaves, cuda),
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
	// local processing time of each process
	double timeSpent;
	// cluster array used in split graph solutions, contains the vertices ids in the full / global graph
	std::vector<long> globalVertexId;
	// number of vertices and edges in the processor subgraph
	long num_vertices, num_edges;

	OutputMessage() : clustering(), numberOfTestedCombinations(0), timeSpent(0.0), globalVertexId(),
			num_vertices(0), num_edges(0) {

	}

	OutputMessage(Clustering &c, long nc, double time, std::vector<long> gVertexId, long n, long m) : clustering(c),
			numberOfTestedCombinations(nc), timeSpent(time), globalVertexId(gVertexId), num_vertices(n), num_edges(m) {

	}

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clustering;
		ar & numberOfTestedCombinations;
		ar & timeSpent;
		ar & globalVertexId;
		ar & num_vertices;
		ar & num_edges;
	}
};

class MPIMessage {
public:
	// Message with the number of heuristic masters
	static const int INPUT_MSG_MPI_PARAMS_TAG = 40;
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
