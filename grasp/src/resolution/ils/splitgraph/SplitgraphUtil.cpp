/*
 * SplitgraphUtil.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: mlevorato
 */

#include <cmath>

#include <boost/math/special_functions/round.hpp>
#include <boost/log/trivial.hpp>
#include <boost/lexical_cast.hpp>

#include "include/SplitgraphUtil.h"
#include "util/include/RandomUtil.h"

using namespace util;

namespace resolution {
namespace ils {

SplitgraphUtil::SplitgraphUtil() {
	// TODO Auto-generated constructor stub

}

SplitgraphUtil::~SplitgraphUtil() {
	// TODO Auto-generated destructor stub
}

std::vector<long> SplitgraphUtil::getListOfVeticesInCluster(SignedGraph& g, const Clustering& globalClustering,
		long clusterNumber) {

	long n = g.getN();
	ClusterArray globalClusterArray = globalClustering.getClusterArray();
	std::vector<long> vertexList;
	for(long vx = 0; vx < n; vx++) {
		if(globalClusterArray[vx] == clusterNumber) {
			//BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Choosing vertex " << vx << " to be moved, it belongs to cluster" << clusterToMove << ".";
			vertexList.push_back(vx);
		}
	}
	return vertexList;
}

std::vector<Coordinate> SplitgraphUtil::obtainListOfClustersFromProcess(SignedGraph& g,
		const Clustering& globalClustering, int processNumber) {

	long nc = globalClustering.getNumberOfClusters();
	ClusterArray globalCluster = globalClustering.getClusterArray();
	// to which process a global cluster belongs to?
	const std::vector<unsigned int> clusterProcessOrigin = globalClustering.getClusterProcessOrigin();

	std::vector<Coordinate> clusterList;
	for(long clusterNum = 0; clusterNum < nc; clusterNum++) {
		if(clusterProcessOrigin[clusterNum] == processNumber) {
			clusterList.push_back(Coordinate(clusterNum, processNumber, globalClustering.getClusterSize(clusterNum)));
		}
	}

	return clusterList;
}

Imbalance SplitgraphUtil::calculateExternalImbalanceSumBetweenProcesses(const ImbalanceMatrix& processClusterImbMatrix) {
	double externalImbalancePosSum = 0.0;
	double externalImbalanceNegSum = 0.0;
	for (unsigned i = 0; i < processClusterImbMatrix.pos.size1(); ++i) {
		for (unsigned j = 0; j < processClusterImbMatrix.pos.size2(); ++j) {
			if(i != j) {
				externalImbalancePosSum += processClusterImbMatrix.pos(i, j);
				externalImbalanceNegSum += processClusterImbMatrix.neg(i, j);
			}
		}
	}
	return Imbalance(externalImbalancePosSum, externalImbalanceNegSum);
}

Imbalance SplitgraphUtil::calculateInternalImbalanceSumOfAllProcesses(std::vector<Imbalance>& internalProcessImbalance) {
	double internalImbalancePosSum = 0.0;
	double internalImbalanceNegSum = 0.0;
	for(unsigned i = 0; i < internalProcessImbalance.size(); i++) {
		internalImbalancePosSum += internalProcessImbalance[i].getPositiveValue();
		internalImbalanceNegSum += internalProcessImbalance[i].getNegativeValue();
	}
	return Imbalance(internalImbalancePosSum, internalImbalanceNegSum);
}

Imbalance SplitgraphUtil::calculateProcessInternalImbalance(SignedGraph *g, Clustering& c,
		unsigned int processNumber) {

	double positiveSum = 0.0, negativeSum = 0.0;
	int nc = c.getNumberOfClusters();
	int n = g->getN();
	ClusterArray myCluster = c.getClusterArray();
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
	UndirectedGraph::edge_descriptor e;

	// creates a bool array where every cluster in process 'processNumber' is marked with 1
	const std::vector<unsigned int> clusterProcessOrigin = c.getClusterProcessOrigin();
	std::vector<bool> processContainsCluster(nc, false);
	for(long k = 0; k < nc; k++) {
		if(clusterProcessOrigin[k] == processNumber) {
			processContainsCluster[k] = true;
		}
	}
	// creates a bool array where every vertex in process 'processNumber' is marked with 1
	std::vector<bool> processContainsVertex(n, false);
	for(long i = 0; i < n; i++) {
		if(processContainsCluster[myCluster[i]]) {  // the vertex belongs to this process
			processContainsVertex[i] = true;
		}
	}

	// For each vertex i that belongs to process 'processNumber'
	for(long i = 0; i < n; i++) {
		if(processContainsVertex[i]) {
			long ki = myCluster[i];
			UndirectedGraph::out_edge_iterator f, l;
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(i, g->graph); f != l; ++f) {
				e = *f;
				double weight = ew[e].weight;
				long j = target(*f, g->graph);
				if(processContainsVertex[j]) {  // only takes into account vertices inside the process
					long kj = myCluster[j];
					bool sameCluster = (ki == kj);
					if(weight < 0 and sameCluster) {  // negative edge
						// i and j are in the same cluster
						negativeSum += fabs(weight);
					} else if(weight > 0 and (not sameCluster)) {  // positive edge
						// i and j are NOT in the same cluster
						positiveSum += weight;
					}
				}
			}
		}
	}
	return Imbalance(positiveSum, negativeSum);
}

/**
  *  Calculates the positive and negative degrees of each vertex v in clusterX *relative to clusterX only*.
*/
std::list<VertexDegree> SplitgraphUtil::calculateDegreesInsideCluster(SignedGraph *g,
		Clustering& globalClustering, long clusterX) {

	std::vector<long> verticesInClusterX = getListOfVeticesInCluster(*g, globalClustering, clusterX);
	ClusterArray myCluster(g->getN(), Clustering::NO_CLUSTER);
	for(std::vector<long>::iterator it = verticesInClusterX.begin(); it != verticesInClusterX.end(); it++) {
		long v = *it;
		myCluster[v] = 1;
	}
	std::list<VertexDegree> degreesInsideClusterX;
	for(std::vector<long>::iterator it = verticesInClusterX.begin(); it != verticesInClusterX.end(); it++) {
		long v = *it;
		// calculates the positive and negative degrees of vertex v *relative to clusterX only*
		double negDegree = fabs(g->getNegativeEdgeSumBetweenVertexAndClustering(v, myCluster));
		double posDegree = g->getPositiveEdgeSumBetweenVertexAndClustering(v, myCluster);
		degreesInsideClusterX.push_back(VertexDegree(v, posDegree, negDegree));
	}
	return degreesInsideClusterX;
}

long SplitgraphUtil::chooseRandomVertex(std::list<VertexDegree>& vertexList, long x) {
	// Generates a random number between 1 and x
	// distribution that maps to 1..x
	if(x - 1 < 0) {
		x++;
	}
	// random number generators used in loop randomization
	RandomUtil randomUtil;
	unsigned int selectedVertexSetIndex = randomUtil.next(0, x - 1);
	VertexDegree selectedVertex;

	// Finds the Vertex
	std::list<VertexDegree>::iterator pos = vertexList.begin();
	std::advance(pos, selectedVertexSetIndex);
	selectedVertex = *pos;
	vertexList.erase(pos);
	return selectedVertex.id;
}

std::vector<Coordinate> SplitgraphUtil::obtainListOfImbalancedClusters(SignedGraph& g,
		Clustering& globalClustering) {

	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	// Cluster to cluster matrix containing positive / negative contribution to imbalance
	long nc = globalClustering.getNumberOfClusters();
	ClusterArray globalCluster = globalClustering.getClusterArray();
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);
	matrix<double> clusterEdgeSumMatrix = zero_matrix<double>(nc, nc);
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = globalCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			bool sameCluster = (globalCluster[targ.id] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				clusterImbMatrix(ki, globalCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				clusterImbMatrix(ki, globalCluster[targ.id]) += fabs(weight);
			}
			clusterEdgeSumMatrix(ki, globalCluster[targ.id]) += fabs(weight);
		}
	}

	std::vector<Coordinate> imbalancedClusterList;
	// For each cluster ki
	for(long ki = 0; ki < nc; ki++) {
		double percentageOfInternalNegativeEdges = 0.0;
		double percentageOfExternalPositiveEdges = 0.0;

		double sumOfInternalEdges = clusterEdgeSumMatrix(ki, ki);
		double sumOfInternalNegativeEdges = clusterImbMatrix(ki, ki);
		double sumOfExternalEdges = 0.0;
		double sumOfPositiveExternalEdges = 0.0;
		for(long kj = 0; kj < nc; kj++) {
			if(ki != kj) {
				sumOfExternalEdges += clusterEdgeSumMatrix(ki, kj) + clusterEdgeSumMatrix(kj, ki);
				sumOfPositiveExternalEdges += clusterImbMatrix(ki, kj) + clusterImbMatrix(kj, ki);
			}
		}
		if(sumOfInternalEdges > 0) {
			percentageOfInternalNegativeEdges = (100 * sumOfInternalNegativeEdges) / sumOfInternalEdges;  // item A
		}
		if(sumOfExternalEdges > 0) {
			percentageOfExternalPositiveEdges = (100 * sumOfPositiveExternalEdges) / sumOfExternalEdges;  // item B
		}
		double imbalancePercentage = max(percentageOfInternalNegativeEdges, percentageOfExternalPositiveEdges);
		double imbalance = sumOfInternalNegativeEdges + sumOfPositiveExternalEdges;
		//double imbalancePercentage = percentageOfExternalPositiveEdges;
		//imbalancedClusterList.push_back(Coordinate(ki, globalClustering.getClusterSize(ki), imbalancePercentage));
		// WHEN MOVING CLUSTERS BETWEEN PROCESSES, ORDERING BY CLUSTER IMBALANCE,
		// THE USE OF ABSOLUTE IMBALANCE VALUES (BELOW) YIELDS BETTER RESULTS THAN PERCENTUAL IMBALANCE (ABOVE)
		imbalancedClusterList.push_back(Coordinate(ki, globalClustering.getClusterSize(ki), imbalance));
	}
	return imbalancedClusterList;
}

std::vector<Coordinate> SplitgraphUtil::obtainListOfOverloadedProcesses(SignedGraph& g,
		const ProcessClustering& processClustering, long maximumNumberOfVertices) {
	long n = g.getN();
	const Clustering& splitgraphClustering = processClustering.getSplitgraphClustering();
	int numberOfProcesses = splitgraphClustering.getNumberOfClusters();
	// A vertex-overloaded process is a process with more than (maximumNumberOfVertices) vertices.
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] List of overloaded processes (with more than "
			<< maximumNumberOfVertices << " vertices): ";

	std::vector<Coordinate> overloadedProcessList;
	stringstream processListStr;
	// For each process px
	for(long px = 0; px < numberOfProcesses; px++) {
		if(splitgraphClustering.getClusterSize(px) > maximumNumberOfVertices) {
			overloadedProcessList.push_back(Coordinate(px, 0, splitgraphClustering.getClusterSize(px)));
			processListStr << px << " ";
		}
	}
	if(processListStr.str() == "") {
		BOOST_LOG_TRIVIAL(debug) << "None.";
	} else {
		BOOST_LOG_TRIVIAL(debug) << processListStr.str() << " ";
	}
	return overloadedProcessList;
}

std::vector<Coordinate> SplitgraphUtil::obtainListOfOverloadedProcesses(SignedGraph& g,
		const ProcessClustering& processClustering) {

	long n = g.getN();
	const Clustering& splitgraphClustering = processClustering.getSplitgraphClustering();
	int numberOfProcesses = splitgraphClustering.getNumberOfClusters();
	// A vertex-overloaded process is a process with more than (n / numberOfProcesses) vertices.
	long numberOfEquallyDividedVertices = (long)ceil(n / (double)numberOfProcesses);
	return obtainListOfOverloadedProcesses(g, processClustering, numberOfEquallyDividedVertices);
}

ImbalanceMatrix SplitgraphUtil::calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster,
		std::vector< pair<long, double> >& vertexImbalance, const int& numberOfProcesses) {

	long nc = numberOfProcesses;
	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	// Process to process matrix containing positive / negative contribution to imbalance
	ImbalanceMatrix clusterImbMatrix(nc);
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	// Vector containing each vertex contribution to imbalance: vertexImbalance

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			bool sameCluster = (myCluster[targ.id] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				negativeSum += fabs(weight);
				clusterImbMatrix.neg(ki, myCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
				clusterImbMatrix.pos(ki, myCluster[targ.id]) += fabs(weight);
			}
		}
		vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
	return clusterImbMatrix;
}

void SplitgraphUtil::updateProcessToProcessImbalanceMatrix(SignedGraph& g,
			const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix, const int& numberOfProcesses) {

	long nc = numberOfProcesses;
	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i in listOfModifiedVertices
	for(long item = 0; item < listOfModifiedVertices.size(); item++) {
		long i = listOfModifiedVertices[item];
		long old_ki = previousSplitgraphClusterArray[i];
		long new_ki = newSplitgraphClusterArray[i];
		UndirectedGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
			e = *f;
			Vertex src = source(e, g.graph), targ = target(e, g.graph);
			double weight = ew[e].weight;
			// Etapa 1: subtracao dos imbalances antigos
			long old_kj = previousSplitgraphClusterArray[targ.id];
			bool sameCluster = (old_kj == old_ki);
			// TODO VERIFICAR SE DEVE SER RECALCULADO TAMBEM O PAR OPOSTO DA MATRIZ: (KJ, KI)
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j were in the same cluster
				processClusterImbMatrix.neg(old_ki, old_kj) -= fabs(weight);
				processClusterImbMatrix.neg(old_kj, old_ki) -= fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j were NOT in the same cluster
				processClusterImbMatrix.pos(old_ki, old_kj) -= fabs(weight);
				processClusterImbMatrix.pos(old_kj, old_ki) -= fabs(weight);
			}
			// Etapa 2: acrescimo dos novos imbalances
			long new_kj = newSplitgraphClusterArray[targ.id];
			sameCluster = (new_kj == new_ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are now in the same cluster
				processClusterImbMatrix.neg(new_ki, new_kj) += fabs(weight);
				processClusterImbMatrix.neg(new_kj, new_ki) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				processClusterImbMatrix.pos(new_ki, new_kj) += fabs(weight);
				processClusterImbMatrix.pos(new_kj, new_ki) += fabs(weight);
			}
		}
		// NAO EH NECESSARIO ATUALIZAR O VERTEX IMBALANCE ABAIXO, POIS A BUSCA LOCAL (VND) EH FIRST IMPROVEMENT
		// vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
}

Coordinate SplitgraphUtil::findMaximumElementInMatrix(ImbalanceMatrix &mat) {
	int x = 0, y = 0;
	double maxValue = -DBL_MAX;

	// TODO INCLUDE HANDLING FOR MATH PRECISION DOUBLE (EPSILON)
	for(int i = 0; i < mat.pos.size1(); i++) {
		for(int j = 0; j < mat.pos.size2(); j++) {
			if((mat.pos(i, j) + mat.neg(i, j)) > maxValue) {
				maxValue = mat.pos(i, j) + mat.neg(i, j);
				x = i;
				y = j;
			}
		}
	}
	return Coordinate(x, y, maxValue);
}

std::vector< Coordinate > SplitgraphUtil::getMatrixElementsAsList(ImbalanceMatrix &mat) {
	std::vector< Coordinate > returnList;

	for(int i = 0; i < mat.pos.size1(); i++) {
		for(int j = 0; j < mat.pos.size2(); j++) {
			returnList.push_back(Coordinate(i, j, (mat.pos(i, j) + mat.neg(i, j))));
		}
	}
	return returnList;
}

long SplitgraphUtil::findMostImbalancedVertexInProcessPair(SignedGraph& g,
		const ClusterArray& splitGraphCluster, const ClusterArray& globalCluster, Coordinate processPair) const {

	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	double maxImbalance = -DBL_MAX;
	long maxVertex = 0;

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = globalCluster[i];
		UndirectedGraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// only processes vertexes belonging to the process pair
		if((splitGraphCluster[i] == processPair.x) or (splitGraphCluster[i] == processPair.y)) {
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
				e = *f;
				Vertex src = source(e, g.graph), targ = target(e, g.graph);
				double weight = ew[e].weight;
				long kj = globalCluster[targ.id];
				bool sameCluster = (kj == ki);
				if(weight < 0 and sameCluster) {  // negative edge
					// i and j are in the same cluster
					negativeSum += fabs(weight);
				} else if(weight > 0 and (not sameCluster)){ // and ((kj == processPair.x) or (kj == processPair.y))) {  // positive edge
					// only counts the imbalance between the clusters / processes in the processPair (x and y)
					// i and j are NOT in the same cluster
					positiveSum += weight;
				}
			}
			// checks if the imbalance caused by vertex i is the biggest one
			// TODO INCLUDE HANDLING FOR MATH PRECISION DOUBLE (EPSILON)
			if(positiveSum + negativeSum > maxImbalance) {
				maxImbalance = positiveSum + negativeSum;
				maxVertex = i;
			}
		}
	}
	// BOOST_LOG_TRIVIAL(info) << "Vertex " << maxVertex << " has imbalance sum of " << maxImbalance;
	return maxVertex;
}

/**
  *  Returns a list containing the vertices that belong to the positive clique C+,
  *  found inside global cluster X. This is a greedy heuristic.
*/
std::vector<long> SplitgraphUtil::findPseudoCliqueC(SignedGraph *g, Clustering& globalClustering,
		long clusterX) {

	string strCluster = boost::lexical_cast<std::string>(clusterX);
	BOOST_LOG_TRIVIAL(info) << "Invoking findPseudoCliqueC for global cluster number " << strCluster << "...";
	long n = g->getN();
	long m = g->getM();

	// *** A - Constructive phase: building a maximal clique
	//  Finding a maximal clique is straightforward: Starting with an arbitrary clique
	//  (for instance, a single vertex), grow the current clique one vertex at a time
	//  by iterating over the graph’s remaining vertices, adding a vertex if it is connected
	//  to each vertex in the current clique, and discarding it otherwise.
	//  This algorithm runs in linear time.
	// 1. For every vertex v in clusterX, let Dx_+(v) = number of positive neighbors inside clusterX.
	//   1.1. Obtains a list containing vertex degree info (vertex number and its positive and negative degrees)
	std::list<VertexDegree> degreesInsideClusterX = calculateDegreesInsideCluster(g, globalClustering, clusterX);
	//   1.2. Sorts the vertex list according to a weighted formula of positive and negative degrees (descending order)
	degreesInsideClusterX.sort(weighted_pos_neg_degree_ordering_desc());

	// 2. The positive clique C+ begins with the vertex v that has the biggest Dx_+(v) in clusterX
	//  Two data structures are used to control the vertices inside cliqueC: a vector and a binary cluster array
	std::list<long> cliqueC;
	ClusterArray cliqueCClusterArray(g->getN(), Clustering::NO_CLUSTER);
	long u = chooseRandomVertex(degreesInsideClusterX, boost::math::lround(CLUSTERING_ALPHA * degreesInsideClusterX.size()));
	cliqueC.push_back(u);
	cliqueCClusterArray[u] = 1;
	BOOST_LOG_TRIVIAL(info) << "Global cluster number " << strCluster << " currently has size " << degreesInsideClusterX.size();
	// BOOST_LOG_TRIVIAL(info) << "Building cliqueC with initial vertex " << u << "...";

	// 3. At each step, add to the clique C+ the vertex u (that also belongs to cluster X), where
	//    C+ \union {u} is also a positive clique, and the vertex u is a random vertex between
	//    the first (alpha x ||Dx_+(v)||) in clusterX => randomness factor
	long cliqueSize = 1;
	while (degreesInsideClusterX.size() > 0) { // list != empty
		// Choose vertex u randomly among the first (alpha x |lc|) elements of lc
		// (alpha x |lc|) is a rounded number
		u = chooseRandomVertex(degreesInsideClusterX, boost::math::lround(CLUSTERING_ALPHA * degreesInsideClusterX.size()));
		// only adds u to the cliqueC if C+ \union {u} is also a positive clique
		if(isPseudoClique(g, cliqueC, cliqueCClusterArray, u)) {
			cliqueC.push_back(u);
			cliqueCClusterArray[u] = 1;
			cliqueSize++;
			// BOOST_LOG_TRIVIAL(info) << "Adding vertex " << it->id << " to the cliqueC.";
		}
		// limits the clique size to half of the source cluster size
		if(cliqueSize >= boost::math::lround(degreesInsideClusterX.size() / (double)2.0))  break;
	}
	BOOST_LOG_TRIVIAL(info) << "Built initial cliqueC of size " << cliqueC.size() << ".";
	// 4. The constructive procedure stops when the clique C+ is maximal.
	// That is, the list of vertices in clusterX has been fully traversed.

	// *** B - Unique local search (executed only once)
	// 1. For every vertex v in clique C+, test a 1-2-swap neighbourhood, that is, try to
	//    remove 1 vertex from C+ and insert other 2 vertices that belong to clusterX (but still not to C+),
	//    preserving the property of maximal positive clique.
	/*
	BOOST_LOG_TRIVIAL(info) << "CliqueC local search...";
	bool improvement = false;
	int iter = 0;
	std::vector<long> verticesInClusterX = util.getListOfVeticesInCluster(*g, globalClustering, clusterX);
	do {
		improvement = false;
		iter++;
		// BOOST_LOG_TRIVIAL(info) << "Clique local search iteration " << iter << "...";
		for(std::list<long>::iterator it = cliqueC.begin(); it != cliqueC.end(); it++) {
			long v = *it;
			// try to remove vertex v from tempCliqueC
			// BOOST_LOG_TRIVIAL(info) << "Trying to remove vertex " << v << " from the clique.";
			cliqueCClusterArray[v] = Clustering::NO_CLUSTER;

			// 1-2-swap neighbourhood:
			// try to insert 2 vertices that belong to clusterX but are still not in cliqueC (and are not the vertex v removed)
			//  for every combination of 2 vertices (a and b) in clusterX but not in cliqueC:
			for(std::vector<long>::iterator it_a = verticesInClusterX.begin(); it_a != verticesInClusterX.end(); it_a++) {
				long a = *it_a;
				if((cliqueCClusterArray[a] < 0) and (a != v)) {  // vertex a belongs to clusterX but is not in cliqueC
					for(std::vector<long>::iterator it_b = verticesInClusterX.begin(); it_b != verticesInClusterX.end(); it_b++) {
						long b = *it_b;
						if((cliqueCClusterArray[b] < 0) and (b != v) and (b < a)) {  // vertex b belongs to clusterX but is not in cliqueC
							// check if cliqueC + {a, b} is still a positive clique
							if(isPseudoClique(g, cliqueCClusterArray, cliqueC.size() - 1, a, b)) {
								// finally insert vertices a and b in cliqueC
								// BOOST_LOG_TRIVIAL(info) << "Removed vertex " << v << " from the clique.";
								// BOOST_LOG_TRIVIAL(info) << "Added vertices " << a << " and " << b  << " to the clique.";

								cliqueCClusterArray[a] = 1;
								cliqueCClusterArray[b] = 1;
								cliqueC.push_back(a);
								cliqueC.push_back(b);
								cliqueC.erase(std::remove(cliqueC.begin(), cliqueC.end(), v), cliqueC.end());
								// BOOST_LOG_TRIVIAL(info) << "CliqueC is now of size " << cliqueC.size() << ".";
								improvement = true;
								break;
							}
						}
						if(improvement)  break;
					}
				}
			}
			if(improvement){
				break;  // restarts the vertex v in cliqueC loop since cliqueC was modified
				// TODO replace this line with a time limit verification (10s?)
			} else {  // undo the movement
				cliqueCClusterArray[v] = 1;
			}
		}
	} while(improvement);
	// 2. Repeat this process until C+ has reached a local optimum related to the 1-2-swap neighbourhood.
	BOOST_LOG_TRIVIAL(info) << "Found a pseudoCliqueC of size " << cliqueC.size() << ".";
*/
	// converts the list to a vector of long
	std::vector<long> cliqueCvec(cliqueC.size());
	copy(cliqueC.begin(), cliqueC.end(), cliqueCvec.begin());
	return cliqueCvec;
}

/**
  *  Determines if the set cliqueC \union {u} is a pseudo clique in the graph g.
  *  By definition, cliqueC is already a pseudo clique.
*/
bool SplitgraphUtil::isPseudoClique(SignedGraph *g, std::list<long>& cliqueC,
		ClusterArray& cliqueCClusterArray, long u) {

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
	UndirectedGraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCount = 0, positiveEdgeCount = 0, negativeEdgeCount = 0;

	// validates the edges from vertex u to cliqueC
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(u, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if(cliqueCClusterArray[j] >= 0) {
			totalEdgeCount++;
			if(ew[e].weight > 0) {
				positiveEdgeCount++;
			} else {
				negativeEdgeCount++;
			}
		}
	}
	// Must allow at least POS_EDGE_PERC_RELAX % of positive connections
	// Must allow no more than NEG_EDGE_PERC_RELAX % of negative connections
	double percOfPositiveConnections = positiveEdgeCount / (double)totalEdgeCount;
	double percOfNegativeConnections = negativeEdgeCount / (double)totalEdgeCount;
	bool isPseudoClique = ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
	return isPseudoClique;
}

/**
  *  Determines if the set cliqueC \union {u, v} is a pseudo clique in the graph g.
  *  By definition, cliqueC is already a pseudo clique.
*/
bool SplitgraphUtil::isPseudoClique(SignedGraph *g, ClusterArray& cliqueCClusterArray,
		long cliqueSize, long u, long v) {

	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g->graph);
	UndirectedGraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCountU = 0, positiveEdgeCountU = 0, negativeEdgeCountU = 0;
	long totalEdgeCountV = 0, positiveEdgeCountV = 0, negativeEdgeCountV = 0;

	// validates the edges from vertex u to cliqueC
	UndirectedGraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(u, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if (cliqueCClusterArray[j] >= 0) {
			totalEdgeCountU++;
			if(ew[e].weight < 0) {
				negativeEdgeCountU++;
			} else {
				positiveEdgeCountU++;
			}
		}
	}
	double percOfPositiveConnections = positiveEdgeCountU / (double)totalEdgeCountU;
	double percOfNegativeConnections = negativeEdgeCountU / (double)totalEdgeCountU;
	bool isPseudoClique = ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
	if( not isPseudoClique ) {
		return false;
	}

	// validates the edges from vertex v to cliqueC
	for (boost::tie(f2, l2) = out_edges(v, g->graph); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, g->graph);
		if (cliqueCClusterArray[j] >= 0) {
			totalEdgeCountV++;
			if(ew[e].weight < 0) {
				negativeEdgeCountV++;
			} else {
				positiveEdgeCountV++;
			}
		}
	}
	percOfPositiveConnections = positiveEdgeCountV / (double)totalEdgeCountV;
	percOfNegativeConnections = negativeEdgeCountV / (double)totalEdgeCountV;
	return ( percOfPositiveConnections >= POS_EDGE_PERC_RELAX )
				and ( percOfNegativeConnections <= NEG_EDGE_PERC_RELAX );
}

std::vector<Imbalance> SplitgraphUtil::calculateProcessInternalImbalance(SignedGraph& g,
		ClusterArray& splitGraphCluster, ClusterArray& globalCluster, int numberOfProcesses) {

	long n = g.getN();
	UndirectedGraph::edge_descriptor e;
	boost::property_map<UndirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	std::vector<Imbalance> processInternalImbalance;

	for(int px = 0; px < numberOfProcesses; px++) {
		// For each vertex i
		double positiveSum = double(0.0), negativeSum = double(0.0);
		for(long i = 0; i < n; i++) {
			if(splitGraphCluster[i] == px) {
				long ki = globalCluster[i];
				UndirectedGraph::out_edge_iterator f, l;
				// For each out edge of i
				for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
					e = *f;
					Vertex src = source(e, g.graph), targ = target(e, g.graph);
					long j = targ.id;
					if(splitGraphCluster[j] == px) {
						long kj = globalCluster[targ.id];
						double weight = ew[e].weight;
						bool sameCluster = (kj == ki);
						if(weight < 0 and sameCluster) {  // negative edge
							// i and j are in the same cluster
							negativeSum += fabs(weight);
						} else if(weight > 0 and (not sameCluster)){ // and ((kj == processPair.x) or (kj == processPair.y))) {  // positive edge
							// only counts the imbalance between the clusters / processes in the processPair (x and y)
							// i and j are NOT in the same cluster
							positiveSum += weight;
						}
					}
				}
			}
		}
		processInternalImbalance.push_back(Imbalance(positiveSum, negativeSum));
	}
	stringstream ss;
	ss << "Calulated InternalImbalanceVector: ";
	for(int x = 0; x < numberOfProcesses; x++) {
		ss << processInternalImbalance[x].getValue() << " ";
	}
	BOOST_LOG_TRIVIAL(info) << ss.str();
	return processInternalImbalance;
}

} /* namespace ils */
} /* namespace resolution */
