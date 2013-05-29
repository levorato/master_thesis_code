/*
 * ParallelGrasp.h
 *
 *  Created on: May 27, 2013
 *      Author: mario
 */

#ifndef PARALLELGRASP_H_
#define PARALLELGRASP_H_

#include "Grasp.h"
#include "../../../problem/include/CCProblem.h"
#include <mpi.h>

using namespace problem;

namespace resolution {
namespace grasp {

struct InputMessage {

	static const int TAG = 50;
	string graphInputFileContents;
	float alpha;
	int l;
	int iter;
	int problemType;
	string fileId;

	InputMessage() : graphInputFileContents(),
			alpha(0.0F), l(1), iter(500), problemType(0), fileId("noId") {

	}

	InputMessage(string graphContents, int it, float a, int neigh,
			int pType, string id) : graphInputFileContents(graphContents),
					alpha(a), l(neigh), iter(it), problemType(pType), fileId(id) {

	}

	string toString() {
		stringstream ss;
		ss << "Alpha: " << alpha << "; l = " << l << "; iter = " << iter << "; fileId = " <<
				fileId << "; " << graphInputFileContents << "\n\n";
		return ss.str();
	}
};

class OutputMessage {
public:
	static const int TAG = 60;
	string clusteringAsText;
	float objectiveFunctionValue;

	OutputMessage() : clusteringAsText("No clustering data available."), objectiveFunctionValue(0.0F) {

	}

	OutputMessage(string cluster, float of) : clusteringAsText(cluster),
			objectiveFunctionValue(of) {

	}
};

class ParallelGrasp : resolution::grasp::Grasp {
public:
	static const int TAG = 50;
	static const int LEADER_ID = 0;

	ParallelGrasp();
	virtual ~ParallelGrasp();

	/**
	 * Triggers the parallel execution of the GRASP algorithm using MPI.
	 */
	ClusteringPtr executeGRASP(SignedGraph *g, int iter, float alpha, int l,
			ClusteringProblem& problem, string& fileId, int& np, int& myRank);
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELGRASP_H_ */
