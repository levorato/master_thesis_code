/*
 * ImbalanceGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ImbalanceGainFunction.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"

#include <boost/log/trivial.hpp>

namespace resolution {
namespace construction {

ImbalanceGainFunction::ImbalanceGainFunction(SignedGraph* g) :
		GainFunction::GainFunction(g) {
	// TODO Auto-generated constructor stub

}

ImbalanceGainFunction::~ImbalanceGainFunction() {
	// TODO Auto-generated destructor stub
}

GainCalculation ImbalanceGainFunction::calculateIndividualGain(
		ClusteringProblem& p, Clustering& c, int v) {

	GainCalculation newCalc;
	int n = graph->getN();
	int nc = c.getNumberOfClusters();
	 BOOST_LOG_TRIVIAL(trace)<< "New gain function for vertex " << v << "...";

	if(p.getType() == ClusteringProblem::CC_PROBLEM) {
		return this->calculateIndividualGainCCProblem(p, c, v);
	}
	// BOOST_LOG_TRIVIAL(trace)<< "Traditional gain function...";
	Clustering cTemp = c;
	GainCalculation gainCalculation;
	double min = std::numeric_limits<double>::max();
	gainCalculation.vertex = v;
	gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
	long k = 0;
	if (p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}

	// For each cluster ci...
	// int nc = c.getNumberOfClusters();
	for (unsigned long ci = 0; ci < nc; ci++) {
		// cout << "Cluster " << k << endl;
		cTemp = c;
		cTemp.addNodeToCluster(*graph, p, v, ci);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = ci;
		}
	}
	// For a new cluster ci+1
	// cout << "New cluster" << endl;
	if ((p.getType() == ClusteringProblem::CC_PROBLEM)
					or ((p.getType() == ClusteringProblem::RCC_PROBLEM)
							and (c.getNumberOfClusters() < k))) {
		cTemp = c;
		cTemp.addCluster(*graph, p, v);
		Imbalance imb = cTemp.getImbalance();
		if (imb.getValue() < min) {
			min = imb.getValue();
			gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
		}
	}
	gainCalculation.gainValue = min;
	// cout << "gain(a) = " << min << endl;
	// BOOST_LOG_TRIVIAL(trace) << gainCalculation.clusterNumber << " vs " << newCalc.clusterNumber;
	// BOOST_LOG_TRIVIAL(trace) << gainCalculation.gainValue  << " vs " << newCalc.gainValue;

	return gainCalculation;
}

GainCalculation ImbalanceGainFunction::calculateIndividualGainCCProblem(
		ClusteringProblem& p, Clustering& c, int v) {

	Clustering cTemp = c;
	int n = graph->getN();
	int nc = c.getNumberOfClusters();
	double currentImbalance = c.getImbalance().getValue();
	ClusterArray myCluster = c.getClusterArray();
	for(int e = 0; e < myCluster.size(); e++) {
		if (myCluster[e] == Clustering::NO_CLUSTER) {
			myCluster[e] = nc;
		}
	}
	// Array that stores the sum of edge weights between vertex i and all clusters
	if( (h_VertexClusterPosSum.size1() == 0) or (nc == 0) ) {
		BOOST_LOG_TRIVIAL(info)<< "Initializing sum arrays for the first time in construction phase...";
		BOOST_LOG_TRIVIAL(debug)<< "n = " << n << ", nc = " << nc;
		h_VertexClusterPosSum = zero_matrix<double>(n, (nc + 1));
		h_VertexClusterNegSum = zero_matrix<double>(n, (nc + 1));
		boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(graph->graph));
		ParallelGraph::edge_descriptor e;

		for(int i = 0; i < n; i++) {  // for each vertex i
			ParallelGraph::out_edge_iterator f, l;
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(vertex(i, *(graph->graph)), *(graph->graph)); f != l; ++f) {
				int j = target(*f, *(graph->graph)).local;
				e = *f;
				double weight = ew[e].weight;
				if(weight > 0) {
					h_VertexClusterPosSum(i, myCluster[j]) += fabs(weight);
				} else {
					h_VertexClusterNegSum(i, myCluster[j]) += fabs(weight);
				}
			}
		}
		BOOST_LOG_TRIVIAL(trace)<< "Sum arrays calculated. Calculating destImbArrays...";
	} else {
		BOOST_LOG_TRIVIAL(trace)<< "Reusing sum arrays. Calculating destImbArrays...";
	}

	/*
	 * Na inserção em cluster novo, contar apenas as relações externas entre o vértice i e os vértices
	 * que estão sem cluster. Não incluir as relações internas quando k = nc.
	 * Quando o vértice i for inserido em cluster existente, contar as relações internas a k2 negativas,
	 * bem como as relações externas a k2 positivas com i.
	 */
	std::vector<double> destImbArray(nc + 1, double(0.0));
	for(int k2 = 0; k2 <= nc; k2++) {
		destImbArray[k2] = currentImbalance;
		if(k2 < nc) {
			destImbArray[k2] += h_VertexClusterNegSum(v, k2);
		}
		for(int k = 0; k < nc; k++) {
			if(k != k2) {
				destImbArray[k2] += h_VertexClusterPosSum(v, k);
			}
		}
		destImbArray[k2] += h_VertexClusterPosSum(v, nc);
		// BOOST_LOG_TRIVIAL(trace)<< "destImbArray[" << k2 << "] = " << destImbArray[k2];
	}
	BOOST_LOG_TRIVIAL(trace)<< "Calculating min element...";
	std::vector<double>::iterator iter = std::min_element(destImbArray.begin(), destImbArray.end());
	double min_val = *iter;
	int position = iter - destImbArray.begin();
	 BOOST_LOG_TRIVIAL(trace)<< "Min element is " << min_val << " and destCluster is " << position;
	GainCalculation newCalc;
	newCalc.vertex = v;
	newCalc.clusterNumber = position;
	if(newCalc.clusterNumber == nc) {
		newCalc.clusterNumber = Clustering::NEW_CLUSTER;
	}
	newCalc.gainValue = min_val;

	// Updates the vertex-cluster sum arrays
	int k1 = nc;
	int k2 = position;
	if(newCalc.clusterNumber == Clustering::NEW_CLUSTER) {
		int new_nc = nc + 1;
		h_VertexClusterPosSum.resize(n, (new_nc + 1), true);
		h_VertexClusterNegSum.resize(n, (new_nc + 1), true);

		// vertex v is being moved to a new cluster
		// move a fileira correspondente ao cluster k = nc na matriz de soma, shiftando os dados para a direita (nc + 1)
		for(int v = 0; v < n; v++) {
			h_VertexClusterPosSum(v, (new_nc)) = h_VertexClusterPosSum(v, (nc));
			h_VertexClusterNegSum(v, (new_nc)) = h_VertexClusterNegSum(v, (nc));
		}
		// zera a fileira movida anteriormente
		for(int v = 0; v < n; v++) {
			h_VertexClusterPosSum(v, (nc)) = 0.0;
			h_VertexClusterNegSum(v, (nc)) = 0.0;
		}
		k1 = new_nc;
	}
	ParallelGraph::out_edge_iterator f, l;
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(graph->graph));
	ParallelGraph::edge_descriptor e;
	// For each out edge of v
	for (boost::tie(f, l) = out_edges(vertex(v, *(graph->graph)), *(graph->graph)); f != l; ++f) {
		int j = target(*f, *(graph->graph)).local;
		e = *f;
		double weight = ew[e].weight;
		if(weight > 0) {
			h_VertexClusterPosSum(j, k1) -= fabs(weight);
			h_VertexClusterPosSum(j, k2) += fabs(weight);
		} else {
			h_VertexClusterNegSum(j, k1) -= fabs(weight);
			h_VertexClusterNegSum(j, k2) += fabs(weight);
		}
	}
	return newCalc;
}

void ImbalanceGainFunction::calculateGainList(ClusteringProblem &p,
		Clustering &c, list<GainCalculation>& nodeList) {
	// cout << "Calculating gain list..." << endl;
	list<GainCalculation, allocator<GainCalculation> >::iterator pos;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if (p.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(p);
		k = rp.getK();
	}
	unsigned int i = 0;
	Clustering cTemp = c;
	for (i = 0, pos = nodeList.begin(); i < nodeList.size(); ++pos, ++i) {
		GainCalculation& gainCalculation = *pos;
		unsigned long a = gainCalculation.vertex;
		// cout << "Vertex " << a << endl;
		double min = std::numeric_limits<double>::max();
		gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;

		// For each cluster cl...
		int nc = c.getNumberOfClusters();
		for (unsigned long cl = 0; cl < nc; cl++) {
			// cout << "Cluster " << cl << endl;
			cTemp = c;
			cTemp.addNodeToCluster(*graph, p, a, cl);
			Imbalance imb = cTemp.getImbalance();
			if (imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = cl;
			}
		}
		// For a new cluster cl+1
		// cout << "New cluster" << endl;
		if ((p.getType() == ClusteringProblem::CC_PROBLEM)
				or ((p.getType() == ClusteringProblem::RCC_PROBLEM)
						and (c.getNumberOfClusters() < k))) {
			cTemp = c;
			cTemp.addCluster(*graph, p, a);
			Imbalance imb = cTemp.getImbalance();
			if (imb.getValue() < min) {
				min = imb.getValue();
				gainCalculation.clusterNumber = Clustering::NEW_CLUSTER;
			}
		}
		gainCalculation.gainValue = min;
		// cout << "gain(a) = " << min << endl;
	}
}

int ImbalanceGainFunction::getType() {
	return GainFunction::IMBALANCE;
}

GainFunction::GainFunctionComparison ImbalanceGainFunction::getComparator() {
	return GainFunctionComparison(true);
}

} /* namespace grasp */
} /* namespace resolution */
