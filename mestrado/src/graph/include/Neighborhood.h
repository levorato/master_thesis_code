/*
 * Neighborhood.h
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include <vector>
#include <boost/config.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph.h"
#include "Clustering.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

typedef shared_ptr<Clustering> ClusteringPtr;

// Defines the neighborhood list (of cluster lists) and its pointer
class NeighborhoodList : vector<ClusterListPtr> {
public:
	NeighborhoodList(Clustering *c, int n);

	/**
	 * Generates a l-neighborhood for this clustering.
	 * @return NeighborhoodList*
	 */
	NeighborhoodList* generateNeighborhood(int l);

	Clustering* findLocalOptimum(SignedGraph *g, ClusteringProblem* problem);

private:
	int numberOfNodes;
	ClusteringPtr clusteringPtr;
	DigraphPtr digraphPtr;
};
typedef shared_ptr<NeighborhoodList> NeighborhoodListPtr;

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
