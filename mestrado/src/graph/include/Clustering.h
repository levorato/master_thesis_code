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
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph.h"

using namespace boost;
using namespace std;

namespace clusteringgraph {

// uses dynamic_bitset for bool array, a high performance and space saving structure
// based on real bits
// the following array is initially empty and needs to be dynamically intialized
typedef boost::dynamic_bitset<> BoolArray;
// Defines the cluster list and its pointer
// the list is made of boolean arrays, indicating that node i is in the cluster
// TODO verificar se eh necessario armazenar o ponteiro para o array ao inves do array em si
typedef vector<BoolArray> ClusterList;
typedef shared_ptr<ClusterList> ClusterListPtr;
// the modularity matrix: a matrix of float
typedef multi_array<float, 2> ModularityMatrix;

/**
 * This class models a set of clusters of a graph. Its main data structure is
 * the ClusterList, a vector of boolean arrays, where each vector represents
 * a cluster and each boolean array marks if a node is in the cluster or not.
 */
class Clustering {
public:
	/**
	 * Creates a Clustering object with n nodes and no clusters.
	 */
	Clustering(int n);
	/**
	 * Creates a Clustering object with n nodes based on the clusterList.
	 */
	Clustering(ClusterList* clusterList, int numberOfNodes);
	/**
	 * Creates a Clustering object based on the clusterList.
	 */
	Clustering(ClusterList* clusterList);

	virtual ~Clustering();

	/**
	 * Adds a new cluster to the clustering configuration.
	 */
	void addCluster(int vertexList[], unsigned int arraySize);

	/**
	 * Returns the n-th cluster of the list.
	 */
	BoolArray getCluster(int clusterNumber);

	/**
	 * Calculates the modularity matrix for this clustering.
	 */
	void calculateModularityMatrix(SignedGraph* g);

	float gain(const int &a);

	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering();

	/**
	 * Returns the number of nodes in the graph.
	 */
	int getNumberOfNodes();

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	int getNumberOfClusters();

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering *c);

private:
	/** number of nodes in the graph (n) */
	int numberOfNodes;
	/** the cluster list, with dimensions k x n */
	ClusterListPtr clusterListPtr;
	/** the modularity matrix */
	ModularityMatrix modularityMatrix;

	void print(std::ostream& os, ClusterList *l);
};

// TODO implement the gain function according to the gain function
// gc(i) specified in the article.
// See Class ClusteringProblem.
class GainFunctionComparison
{
  Clustering* clustering;
public:
  GainFunctionComparison(Clustering* c)
    { clustering = c; }
    bool operator () ( const int& a, const int& b ) const
    {
      return clustering->gain(a) < clustering->gain(b);
    }
};

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
