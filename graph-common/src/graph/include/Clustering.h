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

#include <boost/shared_ptr.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "Graph.h"
#include "Imbalance.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;
using namespace boost::numeric::ublas;

namespace clusteringgraph {

/**
 * This class models a set of clusters of a graph. Its main data structure is
 * the ClusterList, a vector of boolean arrays, where each vector represents
 * a cluster and each boolean array marks if a node is in the cluster or not.
 */
class Clustering {
public:
	static const long NEW_CLUSTER = -1;
	static const long NO_CLUSTER = -2;

	/* Matrices used by RCC relaxed imbalance calculation */
	matrix<double> positiveSum, negativeSum;

	/**
	 * Creates an empty Clustering object.
	 */
	Clustering();

	/**
	 * Creates an empty Clustering object.
	 */
	Clustering(SignedGraph& g);

	/**
	 * Creates a Clustering object based on another Clustering object.
	 */
	Clustering(const Clustering& clustering);

	/**
	 * Creates a Clustering object based on the clusterArray.
	 */
	Clustering(const ClusterArray cArray, SignedGraph& g, ClusteringProblem &p);

	virtual ~Clustering();

	/**
	 * Adds i to a new cluster to the clustering configuration.
	 * Recalculates the objective function associated with
	 * the clustering, based on the modification.
	 */
	void addCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i);

	/**
	 * Returns the biggest cluster (index) of the cluster (the one with more elements).
	 */
	int getBiggestClusterIndex();

	/**
	 * Adds a node i in cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void addNodeToCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k);

	/**
	 * Removes a node i from cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void removeNodeFromCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k);

	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering(unsigned long n);

	/**
	 * Prints the clustering config on the ostream os.
	 */
	void printClustering(ostream& os, unsigned long n);

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	unsigned long getNumberOfClusters() const;

	/**
	 * Removes the k-th cluster. Attention!!! This method DOES NOT
	 * recalculate the objective function associated with the clustering,
	 * based on the modification.
	 */
	void removeCluster(SignedGraph& g, unsigned long k);

	/**
	 * Calculates the size of the k-th cluster.
	 */
	unsigned long getClusterSize(unsigned long k);

	void setClusterSize(unsigned long k, unsigned long size);

	Imbalance getObjectiveFunctionValue();

	void setObjectiveFunctionValue(Imbalance f);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering& c);

	string toString(unsigned long n);

	/**
	 * Returns the value of the objective function corresponding to this cluster.
	 */
	Imbalance& getImbalance() {
		return imbalance;
	}

	void setImbalance(Imbalance imbalance) {
		this->imbalance = imbalance;
	}

	const ClusterArray& getClusterArray() {
		return clusterArray;
	}

private:
	/** the cluster array, with dimension n */
	ClusterArray clusterArray;
	/** the size of each cluster */
	std::vector<unsigned long> clusterSize;
	/** the value of the objective function corresponding to this cluster */
	Imbalance imbalance;
	/** The problem type this clustering refers to (CC, RCC) */
	int problemType;

	void print(std::ostream& os, ClusterArray& a, unsigned long n);

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clusterArray;
		ar & clusterSize;
		ar & imbalance;
		ar & problemType;
		ar & positiveSum;
		ar & negativeSum;
	}
};

typedef shared_ptr<Clustering> ClusteringPtr;

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
