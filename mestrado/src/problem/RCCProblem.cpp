/*
 * RCCProblem.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/RCCProblem.h"

#include <cmath>
#include <algorithm>
#include <cassert>
#include <boost/log/trivial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost;
using boost::numeric::ublas::detail::equals;
using boost::numeric::ublas::matrix; 

namespace problem {

RCCProblem::RCCProblem() {
	// TODO Auto-generated constructor stub

}

RCCProblem::~RCCProblem() {
	// TODO Auto-generated destructor stub
}

int RCCProblem::getType() const {
	return ClusteringProblem::RCC_PROBLEM;
}

/**
 * Returns he Relaxed Imbalance of a partition P (RI(P)).
 */
Imbalance RCCProblem::objectiveFunction(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	std::vector<unsigned long> cluster = c.getVertexInClusterVector();
	std::vector<unsigned long> clusterIdList = c.getClusterIdList();
	std::vector<unsigned long>::iterator it;
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	const matrix<double> posSum(c.positiveSum);
	const matrix<double> negSum(c.negativeSum);
	c.positiveSum.resize(nc, nc, false);
	c.negativeSum.resize(nc, nc, false);
	c.positiveSum.assign(zero_matrix<double>(nc,nc));
	c.negativeSum.assign(zero_matrix<double>(nc,nc));

	// BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting full obj function calculation." << endl;
	// c.printClustering();

	// For each vertex i
	for(int i = 0; i < n; i++) {
		DirectedGraph::out_edge_iterator f, l;
		DirectedGraph::edge_descriptor e;
		it = find(clusterIdList.begin(),clusterIdList.end(), cluster[i]);
		unsigned long ki = std::distance(clusterIdList.begin(), it);
		assert(ki < nc);
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			double weight = ((Edge*)f->get_property())->weight;
			e = *f;
			Vertex dest = target(e, g.graph);
			int j = dest.id;
			// ignores edge loops
            if(i == j)  continue;
			bool sameCluster = cluster[i] == cluster[j];
			if(sameCluster) {
				if(weight > 0) { // positive edge
					c.positiveSum(ki, ki) += weight;
				} else { // negative edge
					c.negativeSum(ki, ki) += fabs(weight);
				}
			} else {
			    // Finds out to which cluster vertex j belongs
				it = find(clusterIdList.begin(),clusterIdList.end(), cluster[j]);
				unsigned long kj = std::distance(clusterIdList.begin(), it);
				assert(kj < nc);
				if(weight > 0) { // positive edge
					c.positiveSum(ki, kj) += weight;
				} else { // negative edge
					c.negativeSum(ki, kj) += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < nc; k2++) {
			if(k1 < k2) {
				externalSum += min(c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
						c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}
	
	// BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;
	/*
	const double epsilon = std::numeric_limits<double>::epsilon();
	if(posSum.size1() > 0) {
		if(equals(c.positiveSum, posSum, epsilon, epsilon) and equals(c.negativeSum, negSum, epsilon, epsilon)) {
			BOOST_LOG_TRIVIAL(trace) << "Obj function calculation MATCH.";
		} else {
			BOOST_LOG_TRIVIAL(trace) << "Obj function calculation DOES NOT MATCH.";
		}
	} */

	return Imbalance(internalSum, externalSum);
}

// TODO consertar aqui!
Imbalance RCCProblem::calculateDeltaPlusObjectiveFunction(SignedGraph& g, Clustering& c,
		const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	std::vector<unsigned long> cluster = c.getVertexInClusterVector();
	std::vector<unsigned long> clusterIdList = c.getClusterIdList();
	std::vector<unsigned long>::iterator it;

	// matrix constains summation of values between pairs of clusters
	// matrix diagonal contains internal sum values
	assert(nc >= c.positiveSum.size1());
	// If a new cluster has been created...
	if(nc > c.positiveSum.size1()) {
		int previousSize = c.positiveSum.size1();
		//BOOST_LOG_TRIVIAL(trace) << "Resizing matrix from dim " << previousSize << " to dim " << nc;
		//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum before resize: " << c.positiveSum << endl;

		c.positiveSum.resize(nc, nc);
		c.negativeSum.resize(nc, nc);
		// assign zero to each new line, row element
		for(int i = 0; i < nc; i++) {
			for(int j = previousSize; j < nc; j++) {
				c.positiveSum(i, j) = 0.0;
				c.positiveSum(j, i) = 0.0;
				c.negativeSum(i, j) = 0.0;
				c.negativeSum(j, i) = 0.0;
			}
		}
		//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum after resize: " << c.positiveSum << endl;
	}
	//BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting delta plus obj function calculation. k = " << k;
	// c.printClustering();

	// Finds out to which cluster vertex i belongs (currentCluster)
	it = find(clusterIdList.begin(),clusterIdList.end(), cluster[i]);
	unsigned long ki = std::distance(clusterIdList.begin(), it);
	assert(ki < nc);

	DirectedGraph::out_edge_iterator f, l;
	DirectedGraph::edge_descriptor e;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
		double weight = ((Edge*)f->get_property())->weight;
		e = *f;
        Vertex dest = target(e, g.graph);
        int j = dest.id;
		// ignores edge loops
        if(i == j)  continue;

		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		bool sameCluster = cluster[i] == cluster[j];
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			it = find(clusterIdList.begin(),clusterIdList.end(), cluster[j]);
			unsigned long kj = std::distance(clusterIdList.begin(), it);
			assert(kj < nc);
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) += weight;
			} else { // negative edge
				c.negativeSum(ki, kj) += fabs(weight);
			}
		}
	}
	DirectedGraph::in_edge_iterator f2, l2;
	// For each in edge of i => edge (j, i)
	for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {
		double weight = ((Edge*)f2->get_property())->weight;
		Vertex src = source(*f2, g.graph), targ = target(*f2, g.graph);
		int j = src.id;
		// ignores edge loops
		if(i == j)  continue;
		assert(i == targ.id);

		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		bool sameCluster = cluster[i] == cluster[j];
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			it = find(clusterIdList.begin(),clusterIdList.end(), cluster[j]);
			unsigned long kj = std::distance(clusterIdList.begin(), it);
			assert(kj < nc);
			if(weight > 0) { // positive edge
				c.positiveSum(kj, ki) += weight;
			} else { // negative edge
				c.negativeSum(kj, ki) += fabs(weight);
			}
		}
	}

	// recalculates the imbalance based on the matrices
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < nc; k2++) {
			if(k1 < k2) {
				externalSum += min(c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
										c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}

	//BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	//BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;

	return Imbalance(internalSum, externalSum);
}

Imbalance RCCProblem::calculateDeltaMinusObjectiveFunction(SignedGraph& g, Clustering& c,
		const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	std::vector<unsigned long> cluster = c.getVertexInClusterVector();
	std::vector<unsigned long> clusterIdList = c.getClusterIdList();
	std::vector<unsigned long>::iterator it;
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	assert(nc <= c.positiveSum.size1());
	//BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting delta minus obj function calculation. k = " << k;
	// c.printClustering();

	// Find out to which cluster vertex i belongs (currentCluster)
	it = find(clusterIdList.begin(),clusterIdList.end(), cluster[i]);
	unsigned long ki = std::distance(clusterIdList.begin(), it);
	assert(ki < nc);

	DirectedGraph::out_edge_iterator f, l;
	DirectedGraph::edge_descriptor e;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
		double weight = ((Edge*)f->get_property())->weight;
		e = *f;
        Vertex dest = target(e, g.graph);
        int j = dest.id;
		// ignores edge loops
        if(i == j)  continue;
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = cluster[i] == cluster[j];
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			it = find(clusterIdList.begin(),clusterIdList.end(), cluster[j]);
			unsigned long kj = std::distance(clusterIdList.begin(), it);
			assert(kj < nc);
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) -= weight;
			} else { // negative edge
				c.negativeSum(ki, kj) -= fabs(weight);
			}
		}
	}
	// For each in edge of i => edge (j, i)
	DirectedGraph::in_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {
		double weight = ((Edge*)f2->get_property())->weight;
		Vertex src = source(*f2, g.graph), targ = target(*f2, g.graph);
        int j = src.id;
		// ignores edge loops
		if(i == j)  continue;
        assert(i == targ.id);
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = cluster[i] == cluster[j];
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			it = find(clusterIdList.begin(),clusterIdList.end(), cluster[j]);
			unsigned long kj = std::distance(clusterIdList.begin(), it);
			assert(kj < nc);
			if(weight > 0) { // positive edge
				c.positiveSum(kj, ki) -= weight;
			} else { // negative edge
				c.negativeSum(kj, ki) -= fabs(weight);
			}
		}
	}

    // If cluster k has only vertex i, it will be deleted...
    if(c.getClusterSize(k) == 1) {
		// remove line and column corresponding to the cluster being removed (cluster k)
		matrix<double> tempPos(c.positiveSum);
		matrix<double> tempNeg(c.negativeSum);
		//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum before shrink: " << c.positiveSum << endl;
		nc--;

		c.positiveSum.resize(nc, nc, false);
		for(int i = 0, aux_i = 0; i < tempPos.size1(); i++) {
				if(i != k) {
						for(int j = 0, aux_j = 0; j < tempPos.size2(); j++) {
								if(j != k) {
										c.positiveSum(aux_i, aux_j) = tempPos(i, j);
										aux_j++;
								}
						}
						aux_i++;
				}
		}
		c.negativeSum.resize(nc, nc, false);
		for(int i = 0, aux_i = 0; i < tempNeg.size1(); i++) {
				if(i != k) {
						for(int j = 0, aux_j = 0; j < tempNeg.size2(); j++) {
								if(j != k) {
										c.negativeSum(aux_i, aux_j) = tempNeg(i, j);
										aux_j++;
								}
						}
						aux_i++;
				}
		}
		//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum after shrink: " << c.positiveSum << endl;
    }
	// recalculates the imbalance based on the positive and negative matrices
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < c.positiveSum.size1(); k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < c.positiveSum.size1(); k2++) {
			if(k1 < k2) {
				externalSum += min(c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
										c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}

	//BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	//BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	//BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;

	return Imbalance(internalSum, externalSum);
}

} /* namespace problem */
