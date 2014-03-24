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
	static const int NEW_CLUSTER = -1;
	/* Matrices used by RCC relaxed imbalance calculation */
	matrix<double> positiveSum, negativeSum;

	/**
	 * Creates a Clustering object with n vertices.
	 */
	Clustering(int n);

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
	void addCluster(SignedGraph& g, ClusteringProblem& p,
			const unsigned long& i);

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
	void printClustering();

	/**
	 * Prints the clustering config on the ostream os.
	 */
	void printClustering(ostream& os);

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
	 * Returns the size of the k-th cluster.
	 */
	unsigned long getClusterSize(unsigned long k);

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
	Imbalance& getImbalance() {
		return imbalance;
	}

	void setImbalance(Imbalance imbalance) {
		this->imbalance = imbalance;
	}

	/**
	 * Returns the clusterId number to which a vertex vertexId belongs.
	 * TODO: return real cluster position!!!
	 */
	unsigned long getVertexClusterIdNumber(int vertexId) const {
		return vertexInClusterId.at(vertexId);
	}

	const std::vector<unsigned long>& getVertexInClusterVector() const {
		return vertexInClusterId;
	}

	const std::vector<unsigned long>& getClusterIdList() const {
		return clusterIdList;
	}

private:
	/** indicates to which cluster a vertex belongs */
	std::vector<unsigned long> vertexInClusterId;
	/** contains the list of IDs of each existing cluster */
	std::vector<unsigned long> clusterIdList;
	/** cluster id counter - always incremented when a new cluster is created */
	unsigned long clusterIdCount;
	/** contains the number of vertices in each cluster */
	std::vector<unsigned long> clusterSize;
	/** the value of the objective function corresponding to this cluster */
	Imbalance imbalance;
	/** The problem type this clustering refers to (CC, RCC) */
	int problemType;

	void print(std::ostream& os, ClusterList& l);

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & vertexInClusterId;
		ar & clusterIdList;
		ar & clusterIdCount,
		ar & clusterSize,
		ar & imbalance;
		ar & problemType;
		ar & positiveSum;
		ar & negativeSum;
	}
};

typedef shared_ptr<Clustering> ClusteringPtr;

} /* namespace clusteringgraph */
#endif /* CLUSTERING_H_ */
