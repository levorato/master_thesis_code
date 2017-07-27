/*
 * SequentialNeighborhoodSearch.cpp
 *
 *  Created on: May 15, 2013
 *      Author: mario
 */

#include "include/SequentialNeighborhoodSearch.h"
#include "../problem/include/CCProblem.h"
#include "../problem/include/RCCProblem.h"
#include <boost/log/trivial.hpp>

namespace clusteringgraph {

using namespace boost;

SequentialNeighborhoodSearch::SequentialNeighborhoodSearch() {

}

Clustering SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, bool firstImprovementOnOneNeig) {
	unsigned long nc = clustering->getNumberOfClusters();
	numberOfTestedCombinations = 0;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
	}
	return SequentialNeighborhoodSearch::searchNeighborhood(l, g, clustering, problem,
			timeSpentSoFar, timeLimit, randomSeed, myRank, 0, g->getN() - 1, firstImprovementOnOneNeig, k);
}

Clustering SequentialNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialSearchIndex,
		unsigned long finalSearchIndex, bool firstImprovementOnOneNeig, unsigned long k) {

	if(not (initialSearchIndex < g->getN())) {
		BOOST_LOG_TRIVIAL(error) << "Assertion FAILED! initialSearchIndex < g->getN()";
	}
	if(not (finalSearchIndex < g->getN())) {
		BOOST_LOG_TRIVIAL(error) << "Assertion FAILED! finalSearchIndex < g->getN()";
	}

	if (l == 1) {  // 1-opt
		// Sequential search does first improvement in 1-opt, depending on parameter value
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, firstImprovementOnOneNeig, k);
	} else {  // 2-opt is always first improvement
		return this->search2opt(g, clustering, &problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, true, k);
	}
}

} /* namespace clusteringgraph */

namespace clusteringgraph {

SequentialNeighborhoodSearch::~SequentialNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

} /* namespace clusteringgraph */
