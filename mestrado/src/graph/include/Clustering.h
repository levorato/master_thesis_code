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

#include "Graph.h"
#include "Imbalance.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

/**
 * This class models a set of clusters of a graph. Its main data structure is
 * the ClusterList, a vector of boolean arrays, where each vector represents
 * a cluster and each boolean array marks if a node is in the cluster or not.
 */
class Clustering {
public:
	static const int NEW_CLUSTER = -1;

	/**
	 * Creates a Clustering object.
	 */
	Clustering();

	/**
	 * Creates a Clustering object based on the clusterList.
	 */
	Clustering(const Clustering& clustering);

	virtual ~Clustering();

	/**
	 * Adds i to a new cluster to the clustering configuration.
	 * Recalculates the objective function associated with
	 * the clustering, based on the modification.
	 */
	BoolArray addCluster(SignedGraph& g, const ClusteringProblem& p,
			const unsigned long& i);

	/**
	 * Returns the n-th cluster of the list.
	 */
	BoolArray& getCluster(unsigned long clusterNumber);

	/**
	 * Adds a node i in cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void addNodeToCluster(SignedGraph& g, const ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k);

	/**
	 * Removes a node i from cluster k. Recalculates the objective
	 * function associated with the clustering, based on the
	 * modification.
	 */
	void removeNodeFromCluster(SignedGraph& g, const ClusteringProblem& p,
			const unsigned long& i, const unsigned long& k);

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
	unsigned long getNumberOfClusters();

	/**
	 * Removes the k-th cluster. Attention!!! This method DOES NOT
	 * recalculate the objective function associated with the clustering,
	 * based on the modification.
	 */
	void removeCluster(SignedGraph& g, unsigned long k);

	/**
	 * Calculates the size of the k-th cluster.
	 */
	unsigned long clusterSize(unsigned long k);

	Imbalance getObjectiveFunctionValue();

	void setObjectiveFunctionValue(Imbalance f);

	/**
	 * Verifies if this clustering object equals another clustering object.
	 * @return bool
	 */
	bool equals(Clustering& c);

	string toString();

	/**
	 * Returns the value of the objective function corresponding to this cluster.
	 */
	const Imbalance& getImbalance() const {
		return imbalance;
	}

	void setImbalance(const Imbalance& imbalance) {
		this->imbalance = imbalance;
	}

	const ClusterList& getClusterList() {
		return clusterList;
	}

private:
	/** the cluster list, with dimensions k x n */
	ClusterList clusterList;
	/** the value of the objective function corresponding to this cluster */
	Imbalance imbalance;
	/** The problem type this clustering refers to (CC, RCC) */
	int problemType;

	void print(std::ostream& os, ClusterList& l);

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & clusterList;
		ar & imbalance;
		ar & problemType;
	}
};

typedef shared_ptr<Clustering> ClusteringPtr;

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
