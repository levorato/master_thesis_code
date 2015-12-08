/*
 * Neighborhood.h
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include <vector>
#include <boost/unordered_set.hpp>
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
        NeighborhoodSearch() : numberOfTestedCombinations(0), h_isNeighborCluster(),
    		h_vertexClusterPosSum(), h_vertexClusterNegSum(), bestClustering() {

        }

        virtual ~NeighborhoodSearch() {

        }

        /**
         * Traverses the l-neighborhood for the given clustering, searching
         * for the best value of objective function.
         * @return NeighborhoodList*
         */
        virtual Clustering searchNeighborhood(int l, SignedGraph* g,
                        Clustering* clustering, ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        int myRank, bool firstImprovementOnOneNeig) = 0;

		long getNumberOfTestedCombinations() {
			return numberOfTestedCombinations;
		}

		/**
		 * Resets internal auxiliary matrices (used between Metaheuristic iterations)
		 * and updates the best clustering based on the constructive phase result.
		 */
		void reset(Clustering coClustering) {
			this->h_isNeighborCluster.clear();
			this->h_vertexClusterNegSum.clear();
			this->h_vertexClusterPosSum.clear();
			bestClustering = coClustering;
		}

		virtual bool isParallel() {
			return false;
		}

protected:

    	/**
    	 * Searches the l-neighborhood of the given clustering, by removing a vertex
    	 * from a cluster (vertex index is between initialSearchIndex and finalSearchIndex)
    	 * and adding it to any other cluster available (and also a new cluster). In the case
    	 * of 2-opt, the initial and final cluster indices only apply to the first vertex in
    	 * the search (not the second one).
    	 * @param k the maximum number of clusters of RCC Problem (optional).
    	 */
    	virtual Clustering searchNeighborhood(int l, SignedGraph* g,
    					Clustering* clustering, ClusteringProblem& problem,
    					double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
    					int myRank,	unsigned long initialSearchIndex, unsigned long finalSearchIndex,
    					bool firstImprovementOnOneNeig, unsigned long k) = 0;

        virtual Clustering search1opt(SignedGraph* g,
                        Clustering* clustering, ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        int myRank, unsigned long initialSearchIndex,
                		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k);

        virtual Clustering search1optCCProblem(SignedGraph* g,
                        Clustering* clustering, ClusteringProblem& problem,
                        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
                        int myRank, unsigned long initialSearchIndex,
                		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k,
                		ClusterArray& myCluster,
                		std::vector<long>& isNeighborCluster,
                		std::vector<double>& vertexClusterPosSum,
                		std::vector<double>& vertexClusterNegSum);

		virtual Clustering search2opt(SignedGraph* g,
						Clustering* clustering, ClusteringProblem* problem,
						double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
						int myRank, unsigned long initialSearchIndex,
						unsigned long finalSearchIndex, bool firstImprovement, unsigned long k);

		virtual Clustering search2optCCProblem(SignedGraph* g,
		                Clustering* clustering, ClusteringProblem& problem,
		                double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
		                int myRank, unsigned long initialSearchIndex,
		        		unsigned long finalSearchIndex, bool firstImprovement, unsigned long k,
		        		ClusterArray& myCluster,
						std::vector<long>& isNeighborCluster,
						std::vector<double>& vertexClusterPosSum,
						std::vector<double>& vertexClusterNegSum);

		virtual void updateVertexClusterSumArrays(SignedGraph* g, std::vector<double>& vertexClusterPosSumArray,
				std::vector<double>& vertexClusterNegSumArray, std::vector<long>& isNeighborClusterArray,	Clustering* clustering);

		virtual void updateVertexClusterSumArraysDelta(SignedGraph* g, std::vector<double>& vertexClusterPosSumArray,
				std::vector<double>& vertexClusterNegSumArray, std::vector<long>& isNeighborClusterArray,
				Clustering& clustering, uint nc, uint new_nc, int i, int k1, int k2);

        /**
         * Generates a new cluster by swithcing the node i from cluster k1 to k3, and
         * node j from cluster k2 to k4. Includes the cases of cluster deletion and
         * new cluster creation, when the clusters contain only one element.
         * The parameter n represents the number of nodes in the graph / clustering.
         * If parameter k3 == -1, inserts node i in a new cluster (alone).
         * If parameter k4 == -1, inserts node j in a new cluster (anlone).
         * Returns false if processing failed (e.g. cluster created in RCC Problem).
         */
        bool process2optCombination(SignedGraph& g, Clustering& clustering,
        		ClusteringProblem* problem,
        		int k1, int k2, int k3, int k4,
        		int n, int i, int j) {

			 // cout << "2-opt-comb: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << ", " << i << ", " << j << endl;
			 // clustering->printClustering();
			 // ClusteringPtr cTemp = make_shared < Clustering > (*clustering);
			 int nc = clustering.getNumberOfClusters();
			 // increments number of tested combinations
			 numberOfTestedCombinations++;
			 // the offset caused by cluster deletions
			 // removes node i from cluster1 and inserts in cluster3
			 // TODO check if the removal of node i destroys cluster1
			 // cout << "k3" << endl;
			 clustering.removeNodeFromCluster(g, *problem, i, k1);
			 // recalculates the number of clusters, as one of them may have been removed
			 int newnc1 = clustering.getNumberOfClusters();
			 if(newnc1 < nc) {
					 // cluster k1 has been removed
					 if(k2 >= k1) { k2--; assert(k2 >= 0); }
					 if(k3 >= k1) { k3--; assert(k3 >= 0); }
					 if(k4 >= k1) { k4--; /* assert(k4 >= 0); */ }
			 }
			 if (k3 > k1) {
					 // inserts i in existing cluster k3
					 clustering.addNodeToCluster(g, *problem, i, k3);
			 } else {
					 // inserts i in a new cluster (alone)
				 	 if(problem->getType() == ClusteringProblem::RCC_PROBLEM) {
						 return false;
					 }
					 clustering.addCluster(g, *problem, i);
			 }
			  // cout << "k4" << endl;
			 // removes node j from cluster2 and inserts in cluster4
			 clustering.removeNodeFromCluster(g, *problem, j, k2);
			 int newnc2 = clustering.getNumberOfClusters();
			 if(newnc2 < newnc1) {
					 // cout << "cluster k2 has been removed" << endl;
					 if(k4 >= k2) { k4--; assert(k4 >= 0); }
			 }
			 if (k4 > k2) {
					 // inserts j in existing cluster k4
					 clustering.addNodeToCluster(g, *problem, j, k4);
			 } else {
					 // inserts j in a new cluster (alone)
				 	 if(problem->getType() == ClusteringProblem::RCC_PROBLEM) {
				 		 return false;
				 	 }
					 clustering.addCluster(g, *problem, j);
			 }
			 return true;
		}

	/* Number of tested combinations during neighborhood search */
	unsigned long numberOfTestedCombinations;
	/* Data structures reused during each metaheuristic iteration */
	// Arrays that store the sum of edge weights between vertex i and all clusters
	std::vector<long> h_isNeighborCluster;
	std::vector<double> h_vertexClusterPosSum;
	std::vector<double> h_vertexClusterNegSum;
	// Best clustering solution so far
	Clustering bestClustering;
};

} /* namespace clusteringgraph */
#endif /* NEIGHBORHOOD_H_ */
