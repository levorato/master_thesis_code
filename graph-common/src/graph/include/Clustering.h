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
#include "ParallelBGLSignedGraph.h"
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
	 * Creates a Clustering object based on the clusterArray and
	 * calculates full objective function.
	 */
	Clustering(ClusterArray &cArray, SignedGraph& g, ClusteringProblem &p);

	/**
	 * Creates a Clustering object based on the clusterArray and
	 * objective function supplied.
	 */
	Clustering(ClusterArray &cArray, SignedGraph& g, ClusteringProblem &p,
			double positiveImbalance, double negativeImbalance);

	/**
	 * Creates a Clustering object based on the clusterArray, cluster process origin and
	 * objective function supplied.
	 */
	Clustering(ClusterArray &cArray, SignedGraph& g, ClusteringProblem &p,
			double positiveImbalance, double negativeImbalance,
			std::vector<unsigned int>& clusterProcessOrigin,
			std::vector<Imbalance>& inProcessImbalance);

	virtual ~Clustering();

	/**
	 * Adds i to a new cluster to the clustering configuration.
	 * Recalculates the objective function associated with
	 * the clustering, based on the modification.
	 */
	void addCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i, bool updateImbalance = true);

	/**
	 * Returns the biggest cluster (index) of the cluster (the one with more elements).
	 */
	int getBiggestClusterIndex() const;

	/**
	 * Adds a node i in cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void addNodeToCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k, bool updateImbalance = true);

	/**
	 * Removes a node i from cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void removeNodeFromCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k, bool updateImbalance = true);

	/**
	 * Prints the clustering config on the screen.
	 */
	void printClustering(unsigned long n) const;

	/**
	 * Prints the clustering config on the ostream os.
	 */
	void printClustering(ostream& os, unsigned long n) const;

	/**
	 * Returns the number of clusters in this clustering configuration.
	 */
	const unsigned long getNumberOfClusters() const;

	/**
	 * Removes the k-th cluster. Attention!!! This method DOES NOT
	 * recalculate the objective function associated with the clustering,
	 * based on the modification.
	 */
	void removeCluster(SignedGraph& g, unsigned long k);

	/**
	 * Calculates the size of the k-th cluster.
	 */
	unsigned long getClusterSize(unsigned long k) const;

	void setClusterSize(unsigned long k, unsigned long size);

	Imbalance getObjectiveFunctionValue() const;

	void setObjectiveFunctionValue(Imbalance f);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering& c) const;

	string toString(unsigned long n) const;

	/**
	 * Returns the value of the objective function corresponding to this cluster.
	 */
	const Imbalance& getImbalance() const {
		return imbalance;
	}

	void setImbalance(Imbalance imbalance) {
		this->imbalance = imbalance;
	}

	const ClusterArray& getClusterArray() const {
		return clusterArray;
	}

	/**
	 * Return to which process a cluster belongs to (used in splitgraph clustering only).
	 */
	const std::vector<unsigned int>& getClusterProcessOrigin() const {
		return processOrigin;
	}

	/**
	 * Sets the process number to which a cluster belongs to.
	 */
	void setProcessOrigin(unsigned long clusterNumber, unsigned int processNumber) {
		processOrigin[clusterNumber] = processNumber;
	}

	/**
	 * Returns a list containing the internal process imbalance for each process in splitgraph.
	 */
	const std::vector<Imbalance>& getInternalProcessImbalance() const {
		return internalProcessImbalance;
	}

	/**
	 * Sets the internal imbalance of a specific process.
	 */
	void setProcessImbalance(unsigned int processNumber, Imbalance imb) {
		internalProcessImbalance[processNumber] = imb;
	}

	/**
	 * Deletes all the clusters which belong to a specific process (spligraph mode).
	 */
	void removeAllClustersFromProcess(SignedGraph *g, unsigned int processNumber);

private:
	/** the cluster array, with dimension n */
	ClusterArray clusterArray;
	/** the size of each cluster */
	std::vector<unsigned long> clusterSize;
	/** the value of the objective function corresponding to this cluster */
	Imbalance imbalance;
	/** The problem type this clustering refers to (CC, RCC) */
	int problemType;
	/** To which process a cluster belongs to (used in splitgraph clustering only) */
	std::vector<unsigned int> processOrigin;
	/** Stores the internal imbalance (positive and negative) for a given process (splitgraph only) */
	std::vector<Imbalance> internalProcessImbalance;

	void print(std::ostream& os, const ClusterArray& a, unsigned long n) const;

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
		ar & processOrigin;
		ar & internalProcessImbalance;
	}
};

typedef shared_ptr<Clustering> ClusteringPtr;

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
