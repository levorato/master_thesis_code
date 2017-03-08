/*
 * ConstructClustering.h
 *
 *  Created on: Apr 23, 2014
 *      Author: mlevorato
 */

#ifndef CONSTRUCTCLUSTERING_H_
#define CONSTRUCTCLUSTERING_H_

#include "../../construction/include/GainFunction.h"
#include "graph/include/Graph.h"
#include "graph/include/Clustering.h"
#include "problem/include/ClusteringProblem.h"

#include <fstream>

namespace resolution {
namespace construction {

class ConstructClustering {
public:
	ConstructClustering(string& initPartitionFile, const bool& initPartitionsFromFileForAllItersEnabled, GainFunction* f, const unsigned long& seed, const double& alpha, Clustering* cl = NULL);
	virtual ~ConstructClustering();

	/**
	 * Constructs a clustering in a greedy randomized fashion,
	 * starting from the empty set.
	 * This is the first phase of the metaheuristic.
	 * @param g graph to be used as the base
	 * @param problem the ClusteringProblem object for the objective function calculation
	 * @return Clustering Cc
	 */
	Clustering constructClustering(SignedGraph *g, ClusteringProblem& problem,
			const int& myRank);


///////////////////////////////////////////////////////////////////////////
// New

  /*
   * Constructs a clustering based on a given initial partition file.
   * This filename is given by the command parameter 'init-partition-file' and contains clustering
   * information estimated by a community detection algorithm. Mainly, each line of this file
   * represents a node id (i.e. 1st line represents nodeId1, etc.).
   *
   * In order to construct a clustering object, we need to do the following steps:
   * 1) Create an empty clustering object
   * 2) Get the variable 'clusterArray' of the this empty clustering object
   * 3) Read the initial partition file and get cluster information into 'clusterArray'
   * 4) Construct the variable 'clusterSize' based on the information in 'clusterArray'
   * 5) Compute imbalance value based in the clustering information and assing to the clustering object
   */
	Clustering constructInitialClusteringFromFile(SignedGraph *g, ClusteringProblem& problem,
			const int& myRank);

/////////////////////////////////////////////////////////////////////////////


	double getAlpha() {
		return _alpha;
	}

	int getGainFunctionType() {
		return gainFunction->getType();
	}

	Clustering* getCCclustering() {
		return CCclustering;
	}

	string getInitPartitionFileName(){
		return _initPartitionFileName;
	}

	bool getInitPartitionsFromFileForAllItersEnabled(){
		return _initPartitionsFromFileForAllItersEnabled;
	}

private:
	GainFunction* gainFunction;
	unsigned long randomSeed;
	/**
	 * alpha parameter belonging to the interval [0, 1]
	 * if alpha < 0, the constructClustering method will always choose the first vertex
	 * in the gainFunction list, that is, the one that minimizes the objective (VOTE algorithm).
	 */
	double _alpha;
	Clustering* CCclustering;

	string _initPartitionFileName;
	bool _initPartitionsFromFileForAllItersEnabled;
};

} /* namespace construction */
} /* namespace resolution */

#endif /* CONSTRUCTCLUSTERING_H_ */
