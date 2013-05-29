/*
 * Clustering.h
 *
 *  Created on: May 6, 2013
 *      Author: mario
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <vector>
#include <boost/config.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph.h"

using namespace boost;
using namespace std;

namespace clusteringgraph {

// uses dynamic_bitset for bool array, a high performance and space saving structure
// based on real bits
// the following array is initially empty and needs to be dynamically intialized
typedef boost::dynamic_bitset<> BoolArray;
typedef shared_ptr<BoolArray> BoolArrayPtr;
// Defines the cluster list and its pointer
// the list is made of boolean arrays, indicating that node i is in the cluster
// TODO verificar se eh necessario armazenar o ponteiro para o array ao inves do array em si
typedef vector<BoolArrayPtr> ClusterList;
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

	/**
	 * Creates a Clustering object with n nodes and no clusters.
	 */
	Clustering(int n);
	/**
	 * Creates a Clustering object with n nodes based on the clusterList.
	 */
	Clustering(const Clustering& clustering, int n);
	/**
	 * Creates a Clustering object based on the clusterList.
	 */
	Clustering(ClusterList& clusterList);

	virtual ~Clustering();

	/**
	 * Adds a new cluster to the clustering configuration.
	 */
	void addCluster(int vertexList[], unsigned int arraySize);

	/**
	 * Returns the n-th cluster of the list.
	 */
	BoolArray* getCluster(int clusterNumber);

	/**
	 * Adds a node i in cluster k.
	 */
	void addNodeToCluster(int i, int k);

	/**
	 * Removes a node i from cluster k.
	 */
	void removeNodeFromCluster(int i, int k);

	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering();

	/**
	 * Prints the clustering config on the ostream os.
	 */
	void printClustering(ostream& os);

	/**
	 * Returns the number of nodes in the graph.
	 */
	int getNumberOfNodes();

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	int getNumberOfClusters();

	/**
	 * Removes the k-th cluster.
	 */
	void removeCluster(int k);

	/**
	 * Calculates the size of the k-th cluster.
	 */
	int clusterSize(int k);

	/**
	 * Returns the gain of a given vertex (based on the modularity matrix)
	 * and the number of the cluster where the insertion of the vertex
	 * brings the best gain possible (return type is GainCalculation).
	 */
	GainCalculation gain(SignedGraph& graph, const int &a);

	/**
	 * Returns the value of the objective function corresponding to this cluster.
	 */
	float getObjectiveFunctionValue();

	void setObjectiveFunctionValue(float f);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering& c);

private:
	/** number of nodes in the graph (n) */
	int numberOfNodes;
	/** the cluster list, with dimensions k x n */
	ClusterList clusterList;
	/** the value of the objective function corresponding to this cluster */
	float objectiveFunctionValue;

	void print(std::ostream& os, ClusterList& l);
};

// TODO implement the gain function according to the gain function
// gc(i) specified in the article.
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
