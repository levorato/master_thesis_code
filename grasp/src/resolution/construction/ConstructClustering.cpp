/*
 * ConstructClustering.cpp
 *
 *  Created on: Apr 23, 2014
 *      Author: mlevorato
 */

#include "include/ConstructClustering.h"
#include "../construction/include/VertexSet.h"
#include <boost/log/trivial.hpp>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/round.hpp>
#include "problem/include/RCCProblem.h"
#include <cassert>

namespace resolution {
namespace construction {

using namespace problem;

ConstructClustering::ConstructClustering(string& initPartitionFile, const bool& initPartitionsFromFileForAllItersEnabled, GainFunction* f,
		const unsigned long& seed, const double& alpha, Clustering* cl) :
		_initPartitionFileName(initPartitionFile), _initPartitionsFromFileForAllItersEnabled(initPartitionsFromFileForAllItersEnabled), gainFunction(f), randomSeed(seed), _alpha(alpha), CCclustering(cl) {

}

ConstructClustering::~ConstructClustering() {

}


/////////////////////////////////////////////////////////////////////////////////////////////////
// New

Clustering ConstructClustering::constructInitialClusteringFromFile(SignedGraph *g,
		ClusteringProblem& problem, const int& myRank) {

	long n = g->getN(); // graph size
	Clustering Cc(*g); // Cc = empty
	ClusterArray clArray = Cc.getClusterArray();

///////////////////////////

	//Read initial partition file to fill in the initial partitions
	char line[256]; // Each line in the file represents a cluster id. size of 256 may be too much
	ifstream initFileStream;
	initFileStream.open(_initPartitionFileName.c_str());

	if(initFileStream.is_open()){
		long indx=0;
	   	while ( initFileStream.getline(line, sizeof line) ){
	     	  //cout << line << '\n';
		  long clId = atol(line);
		  clArray[indx]=clId-1; // minus 1: cluster id starts at 1 in membership file
		  indx++;
	    	}	
	        initFileStream.close();
		
    // here "indx" is supposed to be equal to the number of lines in the file
		long nbLineInInitPartitionFile = indx;
		assert(nbLineInInitPartitionFile == n); // if not, quit

	} else {
		cout << "can not open the initial partition file:" << _initPartitionFileName << endl;
	}

  /////////////////////////
	// construct a clustering object from ClusterArray

	/*
  clArray[0] = 1;
	clArray[1] = 0;
	clArray[2] = 0;
  */

	Cc.setClusterArray(clArray, n);

	long maxClusterId = 0;
	for(long i=0; i<n;i++){
    if(clArray[i]>maxClusterId){
      maxClusterId = clArray[i];
	  }
	}
	long numberOfClusters = maxClusterId+1; // Plus 1: cluster ids start at 0
	
  std::vector< std::vector<long> > clusters(numberOfClusters, std::vector<long>());
	for(long i = 0; i < n; i++) {
    // fill in the vector of vector based on clustering information given by initial partition file
	  clusters[clArray[i]].push_back(i); 
	}

	// compute clusters' size. These clusters' size are stored in the variable 'clusterSize'
  // of the class 'Clustering'
	Cc.constructClusterSize(numberOfClusters, clusters);

  // compute the imbalance value based on the new partition information and store it
	Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	//cout << "Imbalance=" << problem.objectiveFunction(*g, Cc).getValue() << endl;

	return Cc;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


Clustering ConstructClustering::constructClustering(SignedGraph *g,
		ClusteringProblem& problem, const int& myRank) {
	Clustering Cc(*g); // Cc = empty
	VertexSet lc(randomSeed, g->getN()); // L(Cc) = V(G)
	BOOST_LOG_TRIVIAL(debug)<< "Construct clustering...\n";
	// 0. Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// if CC Clustering object is not null, skips the construct clustering and returns CC Clustering
	if(CCclustering != NULL) {
		CCclustering->setImbalance(problem.objectiveFunction(*g, *CCclustering));
		return *CCclustering;
	}

	ClusterArray mycluster = Cc.getClusterArray();
	int nc = Cc.getNumberOfClusters();
	for(int e = 0; e < mycluster.size(); e++) {
		if (mycluster[e] == Clustering::NO_CLUSTER) {
			mycluster[e] = nc;
		}
	}

	while (lc.size() > 0) { // lc != empty
		GainCalculation gainCalculation;
		if (_alpha == 1.0) {
			// alpha = 1.0 (completely random): no need to calculate all gains (saves time)
			int i = lc.chooseRandomVertex(lc.size()).vertex;
			gainCalculation = gainFunction->calculateIndividualGain(problem, Cc, i);
		} else {
			// 1. Compute L(Cc): order the elements of the VertexSet class (lc)
			// according to the value of the gain function
			gainFunction->calculateGainList(problem, Cc, lc.getVertexList());

			// 2. Choose i randomly among the first (alpha x |lc|) elements of lc
			// (alpha x |lc|) is a rounded number
			// IMPORTANT: if alpha < 0, the constructClustering method will always
			//             choose the first vertex in the gainFunction list, that is,
			//             the one that minimizes the objective (VOTE algorithm).
			if (_alpha <= 0.0) {
				// chooses first element (min) without ordering (saves time)
				gainCalculation = lc.chooseFirstElement(gainFunction);
			} else {
				lc.sort(gainFunction);
				gainCalculation = lc.chooseRandomVertex(
						boost::math::iround(_alpha * lc.size()));
			}
			// std::cout << "Random vertex between 0 and " << boost::math::iround(alpha * lc.size()) << " is " << i << std::endl;
		}
		// 3. Cc = C union {i}
		// Adds the vertex i to the partial clustering C, in a way so defined by
		// its gain function. The vertex i can be augmented to C either as a
		// separate cluster {i} or as a member of an existing cluster c in C.
		// cout << "Selected vertex is " << i << endl;
		// The gain function ensures no clusters of size > k are created if RCC Problem
		if (gainCalculation.clusterNumber == Clustering::NEW_CLUSTER) {
			// inserts i as a separate cluster
			Cc.addCluster(*g, problem, gainCalculation.vertex);
		} else {
			// inserts i into existing cluster
			Cc.addNodeToCluster(*g, problem, gainCalculation.vertex,
					gainCalculation.clusterNumber);
		}

		// 4. lc = lc - {i}
		// the choosing vertex i automatically removes it from the list
		// Removal already done by the chooseVertex methods above
		// Cc->printClustering();
	}
	// Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	// => Finally: Stops the timer and stores the elapsed time
        timer.stop();
        boost::timer::cpu_times end_time = timer.elapsed();
        double timeSpent = (end_time.wall - start_time.wall)
                         / double(1000000000);
	timer.resume();
	if(timeSpent >= 3600) {
		BOOST_LOG_TRIVIAL(info)<< "GRASP construct clustering done (before post-processing). Time spent: " << timeSpent;
	}
	// Post-processing phase, only for RCC Problem
	if (problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		unsigned long k = rp.getK();
		int n = g->getN();
		// Repartitions the clustering so that initial clustering has exactly k clusters
		while (Cc.getNumberOfClusters() < k) {
			ClusterArray myCluster = Cc.getClusterArray();
			// BOOST_LOG_TRIVIAL(debug)<< "Cc less than k clusters detected.";
			// Splits a cluster in two new clusters
			// 1. Selects the cluster that has more nodes
			int c = Cc.getBiggestClusterIndex();
			// 2. Randomly selects half of the nodes of cl to move to a new cluster
			//    2.1. Populates a vertex set with the vertices in cluster cl
			VertexSet set(randomSeed);
			for (int i = 0; i < n; i++) {
				if (myCluster[i] == c) {
					set.addVertex(i);
				}
			}
			//    2.2. Randomly chooses count/2 vertices, removing them from cluster cl
			//         and adding to new cluster cln
			int half = boost::math::iround(set.size() / 2.0);
			GainCalculation v = set.chooseRandomVertex(set.size());
			Cc.removeNodeFromCluster(*g, problem, v.vertex, c, false);
			Cc.addCluster(*g, problem, v.vertex, false);
			int new_c = Cc.getNumberOfClusters() - 1;
			for (int i = 1; i < half; i++) {
				v = set.chooseRandomVertex(set.size());
				// BOOST_LOG_TRIVIAL(debug)<< "Selecionou " << v.vertex << ", size = " << set.size();
				Cc.removeNodeFromCluster(*g, problem, v.vertex, c, false);
				Cc.addNodeToCluster(*g, problem, v.vertex, new_c, false);
			}
		}
		BOOST_LOG_TRIVIAL(debug)<< "RCC post-processing completed.";
		// Cc.printClustering();
	}
	Cc.setImbalance(problem.objectiveFunction(*g, Cc));
	BOOST_LOG_TRIVIAL(debug)<< "Initial clustering completed. Obj = " << Cc.getImbalance().getValue();
	timer.stop();
        end_time = timer.elapsed();
        timeSpent = (end_time.wall - start_time.wall)
                         / double(1000000000);
        if(timeSpent >= 3600) {
                BOOST_LOG_TRIVIAL(info)<< "Construct clustering done and post-processing done. Total time: " << timeSpent;
        }

	// Cc.printClustering();
	return Cc;
}

} /* namespace construction */
} /* namespace resolution */
