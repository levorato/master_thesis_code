/*
 * ProcessClustering.h
 *
 *  Created on: Mar 10, 2016
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_PROCESSCLUSTERING_H_
#define SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_PROCESSCLUSTERING_H_

#include "graph/include/Clustering.h"

namespace resolution {
namespace ils {

class ProcessClustering {
public:
	ProcessClustering(SignedGraph *g, ClusteringProblem &p, ClusterArray& splitgraphClusterArray) :
		splitgraphClustering(splitgraphClusterArray, *g, p), interProcessImbalanceMatrix() {

	}
	ProcessClustering(SignedGraph *g, ClusteringProblem &p, ClusterArray& splitgraphClusterArray, ImbalanceMatrix& processClusterImbMatrix) :
		splitgraphClustering(splitgraphClusterArray, *g, p), interProcessImbalanceMatrix(processClusterImbMatrix) {

	}
	virtual ~ProcessClustering();

	const Clustering& getSplitgraphClustering() const {
		return splitgraphClustering;
	}

	const ClusterArray& getClusterArray() const {
		return splitgraphClustering.getClusterArray();
	}

	const unsigned long getNumberOfClusters() const {
		return splitgraphClustering.getNumberOfClusters();
	}

	const unsigned long getClusterSize(unsigned long k) const {
		return splitgraphClustering.getClusterSize(k);
	}

	const ImbalanceMatrix& getInterProcessImbalanceMatrix() {
		return interProcessImbalanceMatrix;
	}

	std::vector< std::vector< long > > getListOfVerticesInEachProcess(SignedGraph *g) const {
		std::vector< std::vector< long > > verticesInCluster(splitgraphClustering.getNumberOfClusters(), std::vector< long >());
		long n = g->getN();
		const ClusterArray& splitgraphClusterArray = splitgraphClustering.getClusterArray();
		for(long i = 0; i < n; i++) {
			long k = splitgraphClusterArray[i];
			verticesInCluster[k].push_back(i);
		}
		return verticesInCluster;
	}



private:

	Clustering splitgraphClustering;
	ImbalanceMatrix interProcessImbalanceMatrix;

};

} /* namespace ils */
} /* namespace resolution */

#endif /* SRC_RESOLUTION_ILS_SPLITGRAPH_INCLUDE_PROCESSCLUSTERING_H_ */
