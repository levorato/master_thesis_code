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
#include <list>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost;
using boost::numeric::ublas::detail::equals;
using boost::numeric::ublas::matrix;
using boost::unordered_set;
using namespace clusteringgraph;

namespace problem {

RCCProblem::RCCProblem() :
		k(0) {
	// TODO Auto-generated constructor stub

}

RCCProblem::~RCCProblem() {
	// TODO Auto-generated destructor stub
}

int RCCProblem::getType() const {
	return ClusteringProblem::RCC_PROBLEM;
}

string RCCProblem::getName() {
	return "RCC";
}

/**
 * Returns he Relaxed Imbalance of a partition P (RI(P)).
 */
Imbalance RCCProblem::objectiveFunction(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	ClusterArray myCluster = c.getClusterArray();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	const matrix<double> posSum(c.positiveSum);
	const matrix<double> negSum(c.negativeSum);
	c.positiveSum.resize(nc, nc, false);
	c.negativeSum.resize(nc, nc, false);
	c.positiveSum.assign(zero_matrix<double>(nc, nc));
	c.negativeSum.assign(zero_matrix<double>(nc, nc));

	BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting full obj function calculation.";
	// c.printClustering();
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i
	for (int i = 0; i < n; i++) {
		long ki = myCluster[i];
		assert(ki < nc);
		ParallelGraph::out_edge_iterator f, l;
		ParallelGraph::edge_descriptor e;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(vertex(i, g.graph), g.graph); f != l; ++f) {
			e = *f;
			double weight = ew[e].weight;
			Vertex dest = target(e, g.graph).local;
			int j = dest.id;
			long kj = myCluster[j];
			// ignores edge loops
			if (i == j)
				continue;
			bool sameCluster = (ki == kj);
			if (sameCluster) {
				if (weight > 0) { // positive edge
					c.positiveSum(ki, ki) += weight;
				} else { // negative edge
					c.negativeSum(ki, ki) += fabs(weight);
				}
			} else {
				assert(kj < nc);
				if (weight > 0) { // positive edge
					c.positiveSum(ki, kj) += weight;
				} else { // negative edge
					c.negativeSum(ki, kj) += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	for (unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for (unsigned long k2 = 0; k2 < nc; k2++) {
			if (k1 < k2) {
				externalSum += min(
						c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
						c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}

	// BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;

	//const double epsilon = std::numeric_limits<double>::epsilon();
	//if(posSum.size1() > 0) {
	//	if(equals(c.positiveSum, posSum, epsilon, epsilon) and equals(c.negativeSum, negSum, epsilon, epsilon)) {
	// BOOST_LOG_TRIVIAL(trace) << "Obj function calculation MATCH.";
	//	} else {
	// BOOST_LOG_TRIVIAL(trace) << "Obj function calculation DOES NOT MATCH.";
	//	}
	//}

	return Imbalance(internalSum, externalSum);
}

Imbalance RCCProblem::calculateDeltaPlusObjectiveFunction(SignedGraph& g,
		Clustering& c, const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	ClusterArray myCluster = c.getClusterArray();

	// matrix summation of values between pairs of clusters
	// matrix diagonal contains internal sum values
	assert(nc >= c.positiveSum.size1());
	// If a new cluster has been created...
	if (nc > c.positiveSum.size1()) {
		int previousSize = c.positiveSum.size1();
		 // BOOST_LOG_TRIVIAL(trace) << "Resizing matrix from dim " << previousSize << " to dim " << nc;
		 // BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum before resize: " << c.positiveSum << endl;

		c.positiveSum.resize(nc, nc);
		c.negativeSum.resize(nc, nc);
		// assign zero to each new line, row element
		for (int i = 0; i < nc; i++) {
			for (int j = previousSize; j < nc; j++) {
				c.positiveSum(i, j) = 0.0;
				c.positiveSum(j, i) = 0.0;
				c.negativeSum(i, j) = 0.0;
				c.negativeSum(j, i) = 0.0;
			}
		}
		 // BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum after resize: " << c.positiveSum << endl;
	}
	 // BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting delta plus obj function calculation. k = " << k;
	// c.printClustering();

	unsigned long ki = myCluster[i];
	assert(ki < nc);
	// gets vertex i's new cluster is k

	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	ParallelGraph::out_edge_iterator f, l;
	ParallelGraph::edge_descriptor e;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(vertex(i, g.graph), g.graph); f != l; ++f) {
		e = *f;
		double weight = ew[e].weight;
		Vertex dest = target(e, g.graph).local;
		int j = dest.id;
		long kj = myCluster[j];
		// ignores edge loops
		if (i == j)
			continue;

		// vertex i is being moved to cluster k
		// if applicable, adds the values corresponding to vertex i in its new currentCluster (k)
		bool sameCluster = (kj == k);
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 2 (ADDITION) -----------------
		if (sameCluster) {
			if (weight > 0) { // positive edge
				c.positiveSum(ki, ki) += weight;
			} else { // negative edge
				c.negativeSum(ki, ki) += fabs(weight);
			}
		} else {
			// treats the cases where not all vertices are in the clusters (e.g. construct clustering)
			if (kj < nc and kj != Clustering::NO_CLUSTER) {
				if (weight > 0) { // positive edge
					c.positiveSum(ki, kj) += weight;
				} else { // negative edge
					c.negativeSum(ki, kj) += fabs(weight);
				}
			}
		}
	}
	// recalculates the imbalance based on the matrices
	double internalSum = 0.0, externalSum = 0.0;
	for (unsigned long k1 = 0; k1 < nc; k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for (unsigned long k2 = 0; k2 < nc; k2++) {
			if (k1 < k2) {
				externalSum += min(
						c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
						c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}

	// BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;

	return Imbalance(internalSum, externalSum);
}

Imbalance RCCProblem::calculateDeltaMinusObjectiveFunction(SignedGraph& g,
		Clustering& c, const unsigned long& k, const unsigned long& i) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	ClusterArray myCluster = c.getClusterArray();
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	assert(nc <= c.positiveSum.size1());
	// BOOST_LOG_TRIVIAL(trace) << "[RCCProblem] Starting delta minus obj function calculation. k = " << k;
	// c.printClustering();

	unsigned long ki = myCluster[i];
	assert(ki < nc);

	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	ParallelGraph::out_edge_iterator f, l;
	ParallelGraph::edge_descriptor e;
	// For each out edge of i => edge (i, j)
	for (boost::tie(f, l) = out_edges(vertex(i, g.graph), g.graph); f != l; ++f) {
		e = *f;
		double weight = ew[e].weight;
		Vertex dest = target(e, g.graph).local;
		int j = dest.id;
		long kj = myCluster[j];
		// ignores edge loops
		if (i == j)
			continue;
		// INTERNAL AND EXTERNAL SUM RECALCULATION - STEP 1 (SUBTRACTION) -----------------
		// if applicable, subtracts the values corresponding to vertex i in its current currentCluster (ki)
		bool sameCluster = (kj == k);
		if (sameCluster) {
			if (weight > 0) { // positive edge
				c.positiveSum(ki, ki) -= weight;
			} else { // negative edge
				c.negativeSum(ki, ki) -= fabs(weight);
			}
		} else {
			// Find out to which cluster vertex j belongs
			unsigned long kj = myCluster[j];
			assert(kj < nc and kj != Clustering::NO_CLUSTER);
			if (weight > 0) { // positive edge
				c.positiveSum(ki, kj) -= weight;
			} else { // negative edge
				c.negativeSum(ki, kj) -= fabs(weight);
			}
		}
	}

	// If cluster k has only vertex i, it will be deleted...
	if (c.getClusterSize(k) == 1) {
		// remove line and column corresponding to the cluster being removed (cluster k)
		matrix<double> tempPos(c.positiveSum);
		matrix<double> tempNeg(c.negativeSum);
		// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum before shrink: " << c.positiveSum << endl;
		nc--;

		c.positiveSum.resize(nc, nc, false);
		for (int i = 0, aux_i = 0; i < tempPos.size1(); i++) {
			if (i != k) {
				for (int j = 0, aux_j = 0; j < tempPos.size2(); j++) {
					if (j != k) {
						c.positiveSum(aux_i, aux_j) = tempPos(i, j);
						aux_j++;
					}
				}
				aux_i++;
			}
		}
		c.negativeSum.resize(nc, nc, false);
		for (int i = 0, aux_i = 0; i < tempNeg.size1(); i++) {
			if (i != k) {
				for (int j = 0, aux_j = 0; j < tempNeg.size2(); j++) {
					if (j != k) {
						c.negativeSum(aux_i, aux_j) = tempNeg(i, j);
						aux_j++;
					}
				}
				aux_i++;
			}
		}
		// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum after shrink: " << c.positiveSum << endl;
	}
	// recalculates the imbalance based on the positive and negative matrices
	double internalSum = 0.0, externalSum = 0.0;
	for (unsigned long k1 = 0; k1 < c.positiveSum.size1(); k1++) {
		internalSum += min(c.positiveSum(k1, k1), c.negativeSum(k1, k1));
		for (unsigned long k2 = 0; k2 < c.positiveSum.size1(); k2++) {
			if (k1 < k2) {
				externalSum += min(
						c.positiveSum(k1, k2) + c.positiveSum(k2, k1),
						c.negativeSum(k1, k2) + c.negativeSum(k2, k1));
			}
		}
	}

	// BOOST_LOG_TRIVIAL(trace) << "Calculated value: " << (internalSum + externalSum) << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix positiveSum: " << c.positiveSum << endl;
	// BOOST_LOG_TRIVIAL(trace) << "Matrix negativeSum: " << c.negativeSum << endl;

	return Imbalance(internalSum, externalSum);
}

/**
  * Contabilizar, para cada vertice, o total de arestas que estao contribuindo para o
  * desequilibrio (positivo + negativo). Assim a gente vai poder comparar quem tem
  * mais relacoes em desequilibrio e quem sabe tirar alguma conclusao disso.
  * Esse script será executado sobre os resultados de todos os grafos de intâncias reais.
*/
string RCCProblem::analyzeImbalance(SignedGraph& g, Clustering& c) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	ClusterArray myCluster = c.getClusterArray();
	stringstream ss1, ss2, ss3;
	// os valores de soma entre clusters devem compor uma matriz
	// as diagnonais da matriz contem os valores das somas internas
	const matrix<double> posSum(c.positiveSum);
	const matrix<double> negSum(c.negativeSum);
	c.positiveSum.resize(nc, nc, false);
	c.negativeSum.resize(nc, nc, false);
	c.positiveSum.assign(zero_matrix<double>(nc, nc));
	c.negativeSum.assign(zero_matrix<double>(nc, nc));
	// Cluster to cluster matrix containing positive / negative contribution to imbalance
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	BOOST_LOG_TRIVIAL(info) << "[RCCProblem] Starting imbalance analysis.";
	ss1 << endl << "Imbalance analysis (out edges contribution):" << endl;
	ss1 << "Vertex,PositiveSum,NegativeSum" << endl;
	ss2 << "Imbalance analysis (in edges contribution):" << endl;
	ss2 << "Vertex,PositiveSum,NegativeSum" << endl;

	// For each vertex i
	for (int i = 0; i < n; i++) {
		long ki = myCluster[i];
		assert(ki < nc);
		ParallelGraph::out_edge_iterator f, l;
		ParallelGraph::edge_descriptor e;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(vertex(i, g.graph), g.graph); f != l; ++f) {
			e = *f;
			double weight = ew[e].weight;
			Vertex dest = target(e, g.graph).local;
			int j = dest.id;
			// ignores edge loops
			if (i == j)
				continue;
			bool sameCluster = (myCluster[i] == myCluster[j]);
			if (sameCluster) {
				if (weight > 0) { // positive edge
					c.positiveSum(ki, ki) += weight;
				} else { // negative edge
					c.negativeSum(ki, ki) += fabs(weight);
				}
			} else {
				// Find out to which cluster vertex j belongs
				long kj = myCluster[j];
				assert(kj < nc);
				if (weight > 0) { // positive edge
					c.positiveSum(ki, kj) += weight;
				} else { // negative edge
					c.negativeSum(ki, kj) += fabs(weight);
				}
			}
		}
	}
	double internalSum = 0.0, externalSum = 0.0;
	list<EdgeContribution> edgeList;
	for (unsigned long k1 = 0; k1 < nc; k1++) {
		if(c.positiveSum(k1, k1) < c.negativeSum(k1, k1)) {
			list<EdgeContribution> l = computeEdges(g, c, k1, k1, RCCProblem::POSITIVE_EDGE);
			edgeList.splice(edgeList.end(), l);
			internalSum += c.positiveSum(k1, k1);
		} else {
			list<EdgeContribution> l = computeEdges(g, c, k1, k1, RCCProblem::NEGATIVE_EDGE);
			edgeList.splice(edgeList.end(), l);
			internalSum += c.negativeSum(k1, k1);
		}
		for (unsigned long k2 = 0; k2 < nc; k2++) {
			if (k1 < k2) {
				double posSum = c.positiveSum(k1, k2) + c.positiveSum(k2, k1);
				double negSum = c.negativeSum(k1, k2) + c.negativeSum(k2, k1);
				if(posSum < negSum) {
					list<EdgeContribution> l = computeEdges(g, c, k1, k2, RCCProblem::POSITIVE_EDGE);
					edgeList.splice(edgeList.end(), l);
					list<EdgeContribution> l2 = computeEdges(g, c, k2, k1, RCCProblem::POSITIVE_EDGE);
					edgeList.splice(edgeList.end(), l2);
					externalSum += posSum;
				} else {
					list<EdgeContribution> l = computeEdges(g, c, k1, k2, RCCProblem::NEGATIVE_EDGE);
					edgeList.splice(edgeList.end(), l);
					list<EdgeContribution> l2 = computeEdges(g, c, k2, k1, RCCProblem::NEGATIVE_EDGE);
					edgeList.splice(edgeList.end(), l2);
					externalSum += negSum;
				}
			}
		}
	}
	// computes edges' positive and negative contribution to relaxed imbalance
	std::vector<double> positiveContribution(n, 0.0), negativeContribution(n, 0.0);
	// Out edges contribution for each vertex
	for(list<EdgeContribution>::iterator it = edgeList.begin(); it != edgeList.end(); it++) {
		EdgeContribution edge = *it;
		if(edge.value > 0) {
			positiveContribution[edge.i] += edge.value;
		} else {
			negativeContribution[edge.i] += (-1) * edge.value;
		}
	}
	double sum1 = 0.0, sum2 = 0.0;
	for(int i = 0; i < n; i++) {
		ss1 << i << "," << positiveContribution[i] << "," << negativeContribution[i] << endl;
		sum1 += positiveContribution[i] + negativeContribution[i];
		positiveContribution[i] = 0.0;
		negativeContribution[i] = 0.0;
	}
	// In edges contribution for each vertex
	for(list<EdgeContribution>::iterator it = edgeList.begin(); it != edgeList.end(); it++) {
		EdgeContribution edge = *it;
		if(edge.value > 0) {
			positiveContribution[edge.j] += edge.value;
		} else {
			negativeContribution[edge.j] += (-1) * edge.value;
		}
	}
	for(int i = 0; i < n; i++) {
		ss2 << i << "," << positiveContribution[i] << "," << negativeContribution[i] << endl;
		sum2 += positiveContribution[i] + negativeContribution[i];
	}

	// Matrix containing cluster to cluster contribution to relaxed imbalance
	for (unsigned long k1 = 0; k1 < nc; k1++) {
		if(c.positiveSum(k1, k1) < c.negativeSum(k1, k1)) {
			clusterImbMatrix(k1, k1) = c.positiveSum(k1, k1);
		} else {
			clusterImbMatrix(k1, k1) = - c.negativeSum(k1, k1);
		}
		for (unsigned long k2 = 0; k2 < nc; k2++) {
			if (k1 < k2) {
				if(c.positiveSum(k1, k2) + c.positiveSum(k2, k1) <
						c.negativeSum(k1, k2) + c.negativeSum(k2, k1)) {
					clusterImbMatrix(k1, k2) = c.positiveSum(k1, k2);
					clusterImbMatrix(k2, k1) = c.positiveSum(k2, k1);
				} else {
					clusterImbMatrix(k1, k2) = - (c.negativeSum(k1, k2));
					clusterImbMatrix(k2, k1) = - (c.negativeSum(k2, k1));
				}
			}
		}
	}
	ss3 << endl << "Cluster contribution to imbalance analysis (cluster-cluster matrix):" << endl;
	for(int i = 0; i < nc; i++) {
		for(int j = 0; j < nc; j++) {
			ss3 << clusterImbMatrix(i, j) << ", ";
		}
		ss3 << endl;
	}

	BOOST_LOG_TRIVIAL(info) << "[RCCProblem] Graph analysis done. Obj = " << (internalSum + externalSum);
	return ss1.str() + ss2.str() + ss3.str();
}

list<EdgeContribution> RCCProblem::computeEdges(SignedGraph& g, Clustering& c, int c1, int c2, int edgeType) {
	int nc = c.getNumberOfClusters();
	int n = g.getN();
	ClusterArray myCluster = c.getClusterArray();
	list<EdgeContribution> edgeList;
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i in cluster c1
	for (int i = 0; i < n; i++) {
		if(myCluster[i] == c1) {
			ParallelGraph::out_edge_iterator f, l;
			ParallelGraph::edge_descriptor e;
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(vertex(i, g.graph), g.graph); f != l; ++f) {
				e = *f;
				double weight = ew[e].weight;
				e = *f;
				Vertex dest = target(e, g.graph).local;
				int j = dest.id;
				// ignores edge loops
				if (i == j)
					continue;
				// vertex j must be in cluster c2
				if(myCluster[j] == c2) {
					if (weight > 0 && edgeType == RCCProblem::POSITIVE_EDGE) { // we are looking for pos. edges
						edgeList.push_back(EdgeContribution(i, j, weight));
					} else if(weight < 0 && edgeType == RCCProblem::NEGATIVE_EDGE) { // we are looking for neg. edges
						edgeList.push_back(EdgeContribution(i, j, weight));
					}
				}
			}
		}
	}
	return edgeList;
}

} /* namespace problem */
