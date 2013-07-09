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
#include <boost/make_shared.hpp>
#include <cassert>

#include "Graph.h"
#include "Clustering.h"
#include "../../problem/include/ClusteringProblem.h"

using namespace boost;
using namespace std;
using namespace problem;

namespace clusteringgraph {

// Defines the neighborhood list
class NeighborhoodSearch {
public:
        NeighborhoodSearch() {

        }

        /**
         * Traverses the l-neighborhood for the given clustering, searching
         * for the best value of objective function.
         * @return NeighborhoodList*
         */
        virtual ClusteringPtr searchNeighborhood(int l, SignedGraph* g,
                        Clustering* clustering, const ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        unsigned long numberOfSlaves, int myRank, unsigned long numberOfSearchSlaves) = 0;

protected:

        /**
         * Generates a new cluster by swithcing the node i from cluster k1 to k3, and
         * node j from cluster k2 to k4. Includes the cases of cluster deletion and
         * new cluster creation, when the clusters contain only one element.
         * The parameter n represents the number of nodes in the graph / clustering.
         * If parameter k3 == -1, inserts node i in a new cluster (alone).
         * If parameter k4 == -1, inserts node j in a new cluster (anlone).
         */
        ClusteringPtr process2optCombination(SignedGraph& g, Clustering* clustering, int k1, int k2,
                        int k3, int k4, int n, int i, int j) {

                // cout << "2-opt-comb: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << ", " << i << ", " << j << endl;
                // clustering->printClustering();
                ClusteringPtr cTemp = make_shared < Clustering > (*clustering);
                int nc = cTemp->getNumberOfClusters();
                // the offset caused by cluster deletions
                // removes node i from cluster1 and inserts in cluster3
                // TODO check if the removal of node i destroys cluster1
                // cout << "k3" << endl;
                cTemp->removeNodeFromCluster(g, i, k1);
                // recalculates the number of clusters, as one of them may have been removed
                int newnc1 = cTemp->getNumberOfClusters();
                if(newnc1 < nc) {
                        // cluster k1 has been removed
                        if(k2 >= k1) { k2--; assert(k2 >= 0); }
                        if(k3 >= k1) { k3--; assert(k3 >= 0); }
                        if(k4 >= k1) { k4--; /* assert(k4 >= 0); */ }
                }
                if (k3 > k1) {
                        // inserts i in existing cluster k3
                        cTemp->addNodeToCluster(g, i, k3);
                } else {
                        // inserts i in a new cluster (alone)
                        cTemp->addCluster(g, i);
                }
                // cout << "k4" << endl;
                // removes node j from cluster2 and inserts in cluster4
                cTemp->removeNodeFromCluster(g, j, k2);
                int newnc2 = cTemp->getNumberOfClusters();
                if(newnc2 < newnc1) {
                        // cout << "cluster k2 has been removed" << endl;
                        if(k4 >= k2) { k4--; assert(k4 >= 0); }
                }
                // cout << "Node removed" << endl;
                if (k4 > k2) {
                        // inserts j in existing cluster k4
                        cTemp->addNodeToCluster(g, j, k4);
                } else {
                        // inserts j in a new cluster (alone)
                        cTemp->addCluster(g, j);
                }
                // cout << "Return" << endl;
                return cTemp;
        }
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
