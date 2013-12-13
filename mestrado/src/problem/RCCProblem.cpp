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

// TODO: implement this
/**
 * Returns he Relaxed Imbalance of a partition P (RI(P)).
 */
Imbalance RCCProblem::objectiveFunction(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	c.positiveSum.resize(nc, nc, false);
	c.negativeSum.resize(nc, nc, false);
	c.positiveSum.assign(zero_matrix<double>(nc,nc));
	c.negativeSum.assign(zero_matrix<double>(nc,nc));

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo da funcao objetivo." << endl;

	// For each vertex i
	for(int i = 0; i < n; i++) {
		BoolArray cluster;
		unsigned long ki;
		// Find out to which cluster vertex i belongs
		for(ki = 0; ki < nc; ki++) {
			cluster = c.getCluster(ki);
			if(cluster[i]) {
				break;
			}
		}
		DirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			double weight = ((Edge*)f->get_property())->weight;
			int j = target(*f, g.graph);
			bool sameCluster = cluster[j];
			if(sameCluster) {
				if(weight > 0) { // positive edge
					c.positiveSum(ki, ki) += weight;
				} else { // negative edge
					c.negativeSum(ki, ki) += fabs(weight);
				}
			} else {
			    // Find out to which cluster vertex j belongs
				unsigned long kj = 0;
				for(kj = 0; kj < nc; kj++) {
					cluster = c.getCluster(kj);
					if(cluster[j]) {
						break;
					}
				}
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
		for(unsigned long k2 = 0; k2 < k1; k2++) {
			externalSum += min(c.positiveSum(k1, k2), c.negativeSum(k1, k2));
		}
	}
	
	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum: " << c.positiveSum << endl;

	return Imbalance(internalSum, externalSum);
}

Imbalance RCCProblem::calculateDeltaPlusObjectiveFunction(SignedGraph& g, Clustering& c,
		const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	assert(nc >= c.positiveSum.size1());
	if(nc > c.positiveSum.size1()) {
		int previousSize = c.positiveSum.size1();
		BOOST_LOG_TRIVIAL(trace) << "Redimensionando a matriz da dim " << previousSize << "para a dim " << nc;
		BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum before resize: " << c.positiveSum << endl;

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
		BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum after resize: " << c.positiveSum << endl;
	}
	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo do delta mais da funcao objetivo. k = " << k;

	BoolArray currentCluster;
	unsigned long ki;
	// Find out to which cluster vertex i belongs (currentCluster)
	for(ki = 0; ki < nc; ki++) {
		currentCluster = c.getCluster(ki);
		if(currentCluster[i]) {
			break;
		}
	}
	// gets vertex i's new cluster
	BoolArray newCluster = c.getCluster(k);

	DirectedGraph::out_edge_iterator f, l;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
		double weight = ((Edge*)f->get_property())->weight;
		int j = target(*f, g.graph);

		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		bool sameCluster = newCluster[j];
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) += weight;
				//c.positiveSum(kj, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, kj) += fabs(weight);
				//c.negativeSum(kj, ki) += fabs(weight);
			}
		}
	}
	DirectedGraph::in_edge_iterator f2, l2;
	// For each in edge of i => edge (j, i)
	for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {
		double weight = ((Edge*)f2->get_property())->weight;
		int j = target(*f2, g.graph);

		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		bool sameCluster = newCluster[j];
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				//c.positiveSum(kj, ki) += weight;
				c.positiveSum(ki, kj) += weight;
			} else { // negative edge
				//c.negativeSum(kj, ki) += fabs(weight);
				c.negativeSum(ki, kj) += fabs(weight);
			}
		}
	}
	// recalculates the imbalance based on the matrices
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < k1; k2++) {
			externalSum += min(c.positiveSum(k1, k2), c.negativeSum(k1, k2));
		}
	}

	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum: " << c.positiveSum << endl;

	return Imbalance(internalSum, externalSum);
	// TODO marcar diferenças entre antes e depois para calcular apenas o que mudou nas matrizes
}

Imbalance RCCProblem::calculateDeltaMinusObjectiveFunction(SignedGraph& g, Clustering& c,
		const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	assert(nc <= c.positiveSum.size1());
	if(nc < c.positiveSum.size1()) {
		// remove line and column corresponding to the cluster being removed (cluster k)
		// assign zero to each new line, row element
		matrix<double> tempPos(c.positiveSum);
		matrix<double> tempNeg(c.negativeSum);
		BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum before shrink: " << c.positiveSum << endl;

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
		BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum after shrink: " << c.positiveSum << endl;
	}

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Disparando calculo do delta menos da funcao objetivo. k = " << k;

	BoolArray currentCluster;
	unsigned long ki;
	// Find out to which cluster vertex i belongs (currentCluster)
	for(ki = 0; ki < nc; ki++) {
		currentCluster = c.getCluster(ki);
		if(currentCluster[i]) {
			break;
		}
	}
	// gets vertex i's new cluster
	BoolArray newCluster = c.getCluster(k);

	DirectedGraph::out_edge_iterator f, l;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
		double weight = ((Edge*)f->get_property())->weight;
		int j = target(*f, g.graph);
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = currentCluster[j];
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				c.positiveSum(ki, kj) -= weight;
				//c.positiveSum(kj, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, kj) -= fabs(weight);
				//c.negativeSum(kj, ki) -= fabs(weight);
			}
		}
	}
	// For each in edge of i => edge (j, i)
	DirectedGraph::in_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {
		double weight = ((Edge*)f2->get_property())->weight;
		int j = target(*f2, g.graph);
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = currentCluster[j];
		if(sameCluster) {
			if(weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = 0;
			for(kj = 0; kj < nc; kj++) {
				BoolArray cluster = c.getCluster(kj);
				if(cluster[j]) {
					break;
				}
			}
			if(weight > 0) { // positive edge
				//c.positiveSum(kj, ki) -= weight;
				c.positiveSum(ki, kj) -= weight;
			} else { // negative edge
				//c.negativeSum(kj, ki) -= fabs(weight);
				c.negativeSum(ki, kj) -= fabs(weight);
			}
		}
	}
	// recalculates the imbalance based on the matrices
	double internalSum = 0.0, externalSum = 0.0;
	for(unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for(unsigned long k2 = 0; k2 < k1; k2++) {
			externalSum += min(c.positiveSum(k1, k2), c.negativeSum(k1, k2));
		}
	}

	BOOST_LOG_TRIVIAL(trace) << "Valor calculado: " << (internalSum + externalSum) << endl;
	BOOST_LOG_TRIVIAL(trace) << "Matriz positiveSum: " << c.positiveSum << endl;

	return Imbalance(internalSum, externalSum);
	// TODO marcar diferenças entre antes e depois para calcular apenas o que mudou nas matrizes
}

} /* namespace problem */
