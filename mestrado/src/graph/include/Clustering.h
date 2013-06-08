/*
 * Clustering.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <vector>
#include <map>
#include <boost/config.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "../../util/serialization/dynamic_bitset.hpp"
#include "Graph.h"

using namespace boost;
using namespace std;

namespace clusteringgraph {

/**
 *  uses dynamic_bitset for bool array, a high performance and space saving structure
 *  based on real bits
 *  the following array is initially empty and needs to be dynamically intialized.
 */
typedef dynamic_bitset<> BoolArray;
// Defines the cluster list
// the list is made of boolean arrays, indicating that node i is in the cluster
typedef vector<BoolArray> ClusterList;
typedef struct {
	float value;
	int clusterNumber;
} GainCalculation;

/**
 * This class models a set of clusters of a graph. Its main data structure is
 * the ClusterList, a vector of boolean arrays, where each vector represents
 * a cluster and each boolean array marks if a node is in the cluster or not.
 */
class Clustering {
public:
	static const int NEW_CLUSTER = -1;

	Clustering();

	/**
	 * Creates a Clustering object with n nodes based on the clusterList.
	 */
	Clustering(const Clustering& clustering, int n);

	virtual ~Clustering();

	/**
	 * Adds i to a new cluster to the clustering configuration.
	 * Recalculates the objective function associated with
	 * the clustering, based on the modification.
	 */
	void addCluster(SignedGraph& g, const int& i);

	/**
	 * Returns the n-th cluster of the list.
	 */
	BoolArray& getCluster(int clusterNumber);

	/**
	 * Adds a node i in cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void addNodeToCluster(SignedGraph& g, int i, int k);

	/**
	 * Removes a node i from cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void removeNodeFromCluster(SignedGraph& g, int i, int k);

	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering();

	/**
	 * Prints the clustering config on the ostream os.
	 */
	void printClustering(ostream& os);

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	int getNumberOfClusters();

	/**
	 * Removes the k-th cluster. Attention!!! This method DOES NOT
	 * recalculate the objective function associated with the clustering,
	 * based on the modification.
	 */
	void removeCluster(SignedGraph& g, int k);

	/**
	 * Calculates the size of the k-th cluster.
	 */
	int clusterSize(int k);

	/**
	 * Returns the gain of a given vertex (based on the modularity matrix)
	 * and the number of the cluster where the insertion of the vertex
	 * brings the best gain possible (return type is GainCalculation).
	 */
	GainCalculation& gain(SignedGraph& graph, const int &a);

	/**
	 * Returns the value of the objective function corresponding to this cluster.
	 */
	float getObjectiveFunctionValue();

	void setObjectiveFunctionValue(float f);

	/**
	 * Calculates the delta of the objective function caused by the
	 * insertion of node i in cluster k.
	 */
	float calculateDeltaObjectiveFunction(SignedGraph& g, BoolArray& cluster, const int& i);

	void calculateGainList(SignedGraph &g, list<int>& nodeList);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering& c);

	string toString();

private:
	/** the cluster list, with dimensions k x n */
	ClusterList clusterList;
	/** the value of the objective function corresponding to this cluster */
	float objectiveFunctionValue;
	/** the map of nodes gain value */
	map<int, GainCalculation> gainMap;

	void print(std::ostream& os, ClusterList& l);

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clusterList;
		ar & objectiveFunctionValue;
	}
};

// See Class ClusteringProblem.
class GainFunctionComparison
{
	SignedGraph* graph;
	Clustering* clustering;
public:
  GainFunctionComparison(SignedGraph *g, Clustering* c)
    { graph = g;	clustering = c; }
    bool operator () ( const int& a, const int& b ) const
    {
      return clustering->gain(*graph, a).value < clustering->gain(*graph, b).value;
    }
};

typedef shared_ptr<Clustering> ClusteringPtr;

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
