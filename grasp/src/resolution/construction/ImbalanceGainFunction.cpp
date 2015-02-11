/*
 * ImbalanceGainFunction.cpp
 *
 *  Created on: Jun 18, 2013
 *      Author: mario
 */

#include "include/ImbalanceGainFunction.h"
#include "problem/include/CCProblem.h"
#include "problem/include/RCCProblem.h"
#include "graph/include/Graph.h"

#include <algorithm>
#include <cstdio>
#include <vector>

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

	int nc = c.getNumberOfClusters();
	int n = this->graph->getN();
	ClusterArray myCluster = c.getClusterArray();
	int count = 0;
	for (int e = 0; e < myCluster.size(); e++) {
		if (myCluster[e] == Clustering::NO_CLUSTER) {
			myCluster[e] = nc;
			count++;
		}
	}
	/*
	cout << "Converted " << count << " vertices to nc cluster." << endl;
	cout << "The current number of clusters is " << nc << " and imbalance is "
			<< c.getImbalance().getValue() << endl;
	cout << "Finding the best cluster to move vertex " << v << " to..." << endl;
	*/

	// Array that stores the sum of edge weights between vertex i and all clusters
	std::vector<double> h_VertexClusterPosSum(n * (nc + 1));
	std::vector<double> h_VertexClusterNegSum(n * (nc + 1));
	for (int i = 0; i < n * (nc + 1); i++) {
		h_VertexClusterPosSum[i] = 0.0;
		h_VertexClusterNegSum[i] = 0.0;
	}
	// For each vertex, creates a list of in and out edges
	for (int edge = 0, i = 0; i < n; i++) {  // For each vertex i
		DirectedGraph::out_edge_iterator f, l;  // For each out edge of i
		int count = 0;
		for (boost::tie(f, l) = out_edges(i, graph->graph); f != l; ++f) { // out edges of i
			double weight = ((Edge*) f->get_property())->weight;
			int j = target(*f, graph->graph);
			count++;
			edge++;
			if (weight > 0) {
				h_VertexClusterPosSum[myCluster[j] * n + i] += fabs(weight);
			} else {
				h_VertexClusterNegSum[myCluster[j] * n + i] += fabs(weight);
			}
		}
		DirectedGraph::in_edge_iterator f2, l2;  // For each in edge of i
		for (boost::tie(f2, l2) = in_edges(i, graph->graph); f2 != l2; ++f2) { // in edges of i
			double weight = ((Edge*) f2->get_property())->weight;
			int j = source(*f2, graph->graph);
			count++;
			edge++;
			if (weight > 0) {
				h_VertexClusterPosSum[myCluster[j] * n + i] += fabs(weight);
			} else {
				h_VertexClusterNegSum[myCluster[j] * n + i] += fabs(weight);
			}
		}
	}

	std::vector<float> d_destFunctionValue;
	// result / destination vector
	int numberOfChunks = nc + 1;
	d_destFunctionValue.resize(numberOfChunks);

	for (int k2 = 0; k2 <= nc; k2++) {
		// executes local search for vertex i, moving it to cluster k2
		// Option 1: vertex i is moved from k1 to another existing cluster k2 != k1
		// Option 2: vertex i is moved from k1 to a new cluster (k2 = nc)
		// only calculates the cost of inserting vertex i into cluster k2
		// vertex i is in no cluster
		/*
		 * Na inserção em cluster novo, contar apenas as relações externas entre o vértice i e os vértices
		 * que estão sem cluster. Não incluir as relações internas quando k = nc.
		 * Quando o vértice i for inserido em cluster existente, contar as relações internas a k2 negativas,
		 * bem como as relações externas a k2 positivas com i.
		 */
		// updates thread idx / vertex i to cluster k2 imbalance result
		d_destFunctionValue[k2] = c.getImbalance().getValue();
		if(k2 < nc) {
			d_destFunctionValue[k2] += h_VertexClusterNegSum[v + k2 * n];
			// d_destFunctionValue[k2] -= h_VertexClusterPosSum[v + k2 * n];
		}
		for(int k = 0; k < nc; k++) {
			// printf("%d   %.2f\n", v+k*n, h_VertexClusterPosSum[v + k * n]);
			if(k != k2) {
				d_destFunctionValue[k2] += h_VertexClusterPosSum[v + k * n];
			}
		}

		//printf("%d   %.2f\n", v+nc*n, h_VertexClusterPosSum[v + nc * n]);
		d_destFunctionValue[k2] += h_VertexClusterPosSum[v + nc * n];

		// printf("destFunctionValue[%d] = %.2f\n", k2, d_destFunctionValue[k2]);
	}

	// printf("Begin reduce / post-process...\n");
	std::vector<float>::iterator iter = std::min_element(
			d_destFunctionValue.begin(),
			d_destFunctionValue.begin() + numberOfChunks);
	std::vector<float>::iterator it;
	/*
	for (it = d_destFunctionValue.begin();
			it != d_destFunctionValue.begin() + numberOfChunks; ++it) {
		std::cout << (*it) << ' ';
	}
	std::cout << endl; */

	for (int e = 0; e < myCluster.size(); e++) {
		if (myCluster[e] == nc) {
			myCluster[e] = Clustering::NO_CLUSTER;
		}
	}

	float min_val = *iter;
	// determines the position of the best improvement found in the result vector
	uint position = iter - d_destFunctionValue.begin();
	int resultIdx = position;
	int clusterNumber = resultIdx;
	// printf("Idx = %d: The best src vertex is %d to cluster %d with I(P) = %.2f\n", resultIdx, bestSrcVertex, destcluster, destFunctionValue);
	if (min_val < 0) {
		printf("WARNING: I(P) < 0 !!!\n");
	}
	if (clusterNumber == nc) {
		clusterNumber = Clustering::NEW_CLUSTER;
	}
	// printf("CUDA I(P) = %.2f : vertex %d goes to cluster %d\n", min_val, v, clusterNumber);

	// For each cluster ci...
	nc = c.getNumberOfClusters();
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

	// printf("CPU I(P) = %.2f : vertex %d goes to cluster %d\n", gainCalculation.gainValue, v, gainCalculation.clusterNumber);

	assert(min_val == gainCalculation.gainValue);
	gainCalculation.gainValue = min_val;
	gainCalculation.clusterNumber = clusterNumber;

	// cout << "gain(a) = " << min << endl;
	return gainCalculation;
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
