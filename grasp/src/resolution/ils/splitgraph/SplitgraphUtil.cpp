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

// REIMPLEMENTADA DE ACORDO COM A BOOST PBGL
ImbalanceMatrix SplitgraphUtil::calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster,
		std::vector< pair<long, double> >& vertexImbalance, const int& numberOfProcesses) {
	boost::mpi::communicator world;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] calculateProcessToProcessImbalanceMatrix invoked.";
	for(int px = 1; px < numberOfProcesses; px++) {
		// worker processes, obtain answer via MPI
		InputMessageSplitUtil imsg;
		imsg.myCluster = myCluster;
		imsg.vertexImbalance = vertexImbalance;
		imsg.numberOfProcesses = numberOfProcesses;
		imsg.functionRequested = InputMessageSplitUtil::calculateProcessToProcessImbalanceMatrix;
		world.send(px, MPIMessage::INPUT_MSG_SPLITUTIL_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << px << ".";
	}
	ImbalanceMatrix result = calculateProcessToProcessImbalanceMatrixLocal(g, myCluster, vertexImbalance, numberOfProcesses);
	for(int i = 1; i < numberOfProcesses; i++) {
		OutputMessageSplitUtil omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfProcesses) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ".";
		// MERGE IMBALANCE MATRIX
		result += omsg.imbalanceMatrix;
	}
	return result;
}

ImbalanceMatrix SplitgraphUtil::calculateProcessToProcessImbalanceMatrixLocal(SignedGraph& g, ClusterArray& myCluster,
		std::vector< pair<long, double> >& vertexImbalance, const int& numberOfProcesses) {

	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	BOOST_LOG_TRIVIAL(info) << "[calculateProcessToProcessImbalanceMatrixLocal] Obtaining name_map...";
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g.graph));

	long nc = numberOfProcesses;
	boost::mpi::communicator world;
	int myRank = world.rank();
	ParallelGraph::edge_descriptor e;
	// local subgraph creation
	// LocalSubgraph lsg = make_local_subgraph(*(g.graph));
	// Process to process matrix containing positive / negative contribution to imbalance
	ImbalanceMatrix clusterImbMatrix(nc);
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g.graph));
	// Vector containing each vertex contribution to imbalance: vertexImbalance
	// TODO VERIFICAR A MELHOR FORMA DE REPROGRAMAR A LEITURA DE ARESTAS COM A PARALLEL BGL

	BGL_FORALL_VERTICES(v, *(g.graph), ParallelGraph) {  // For each vertex v
		int i = get(name_map, v); // v.local;   // TODO TROCADO PELO GLOBAL
		BOOST_LOG_TRIVIAL(info) << "[calculateProcessToProcessImbalanceMatrixLocal] Processing local vertex " << v.local << " from process " << v.owner;
		 BOOST_LOG_TRIVIAL(info) << "[calculateProcessToProcessImbalanceMatrixLocal] Processing global vertex " << i << " from process " << v.owner;
		if(owner(v) != myRank) {
			BOOST_LOG_TRIVIAL(error) << "**** VERTICE DE OUTRO PROCESSO! ****";
		} else {
			long ki = myCluster[i];
			ParallelGraph::out_edge_iterator f, l;
			double positiveSum = double(0.0), negativeSum = double(0.0);
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(v, *(g.graph)); f != l; ++f) {
				e = *f;
				// int v_targ_local = target(e, *(g.graph)).local;
				// BOOST_LOG_TRIVIAL(info) << "[calculateProcessToProcessImbalanceMatrixLocal] Processing neighbor " << v_targ_local
				//		<< " from processor " << target(e, *(g.graph)).owner;
				int targ = get(name_map, target(e, *(g.graph)));  // TODO TROCADO PELO GLOBAL
				double weight = ew[e].weight;
				bool sameCluster = (myCluster[targ] == ki);
				if(weight < 0 and sameCluster) {  // negative edge
					// i and j are in the same cluster
					negativeSum += fabs(weight);
					clusterImbMatrix.neg(ki, myCluster[targ]) += fabs(weight);
				} else if(weight > 0 and (not sameCluster)) {  // positive edge
					// i and j are NOT in the same cluster
					positiveSum += weight;
					clusterImbMatrix.pos(ki, myCluster[targ]) += fabs(weight);
				}
			}
			vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
		}
	}
	BOOST_LOG_TRIVIAL(info) << "[calculateProcessToProcessImbalanceMatrixLocal] Done.";
	return clusterImbMatrix;
}

void SplitgraphUtil::updateProcessToProcessImbalanceMatrix(SignedGraph& g,
			const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix, const int& numberOfProcesses) {
	boost::mpi::communicator world;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] updateProcessToProcessImbalanceMatrix invoked.";
	for(int px = 1; px < numberOfProcesses; px++) {
		// worker processes, obtain answer via MPI
		InputMessageSplitUtil imsg(numberOfProcesses);
		imsg.previousSplitgraphClusterArray = previousSplitgraphClusterArray;
		imsg.newSplitgraphClusterArray = newSplitgraphClusterArray;
		imsg.listOfModifiedVertices = listOfModifiedVertices;
		// imsg.processClusterImbMatrix = ImbalanceMatrix(numberOfProcesses);  // zero matrix
		//imsg.numberOfProcesses = numberOfProcesses;
		imsg.functionRequested = InputMessageSplitUtil::updateProcessToProcessImbalanceMatrix;
		world.send(px, MPIMessage::INPUT_MSG_SPLITUTIL_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << px << ".";
	}
	updateProcessToProcessImbalanceMatrixLocal(g, previousSplitgraphClusterArray, newSplitgraphClusterArray,
					listOfModifiedVertices, processClusterImbMatrix, numberOfProcesses);
	for(int i = 1; i < numberOfProcesses; i++) {
		OutputMessageSplitUtil omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfProcesses) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ".";
		// MERGE IMBALANCE MATRIX
		processClusterImbMatrix += omsg.imbalanceMatrix;
	}
}

void SplitgraphUtil::updateProcessToProcessImbalanceMatrixLocal(SignedGraph& g,
			const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			ImbalanceMatrix& processClusterImbMatrix, const int& numberOfProcesses) {

	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	BOOST_LOG_TRIVIAL(info) << "[updateProcessToProcessImbalanceMatrixLocal] Obtaining name_map...";
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g.graph));

	BOOST_LOG_TRIVIAL(info) << "[updateProcessToProcessImbalanceMatrixLocal] Obtaining global_index...";
	mpi_process_group pg = g.graph->process_group();
	boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
	  global_index(pg, num_vertices(*(g.graph)),
				   get(vertex_index, *(g.graph)), get(vertex_global, *(g.graph)));

	long nc = numberOfProcesses;
	ParallelGraph::edge_descriptor e;
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g.graph));
	boost::mpi::communicator world;
	int myRank = world.rank();
	// local subgraph creation
	// LocalSubgraph lsg = make_local_subgraph(*(g.graph));

	// TODO VERIFICAR A MELHOR FORMA DE REPROGRAMAR A LEITURA DE ARESTAS COM A PARALLEL BGL
	// TODO: ALTERAR FORMA DE BUSCAR AS ARESTAS A PARTIR DO NUMERO DO VERTICE, POIS O NUMERO I EH ID GLOBAL DO VERTICE
	// FIXME CADA PROCESSO DEVE ATUALIZAR A PARTE QUE LHE CABE NA MATRIZ, DE MANEIRA DISTRIBUIDA (DISTRIBUTED PROPERTY MAP) ?
	// For each vertex i in listOfModifiedVertices
	// for(long item = 0; item < listOfModifiedVertices.size(); item++) {
	BGL_FORALL_VERTICES(v, *(g.graph), ParallelGraph) {
		int i = get(name_map, v);
		// long i = listOfModifiedVertices[item];
		long old_ki = previousSplitgraphClusterArray[i];
		long new_ki = newSplitgraphClusterArray[i];
		if(owner(v) != myRank) {
			BOOST_LOG_TRIVIAL(error) << "**** VERTICE DE OUTRO PROCESSO! ****";
		}
		else
		if(old_ki != new_ki) {
			// BOOST_LOG_TRIVIAL(debug) << "[updateProcessToProcessImbalanceMatrixLocal] Processing vertex " << i << " from process " << v.owner;
			ParallelGraph::out_edge_iterator f, l;
			// For each out edge of v
			for (boost::tie(f, l) = out_edges(v, *(g.graph)); f != l; ++f) {
				e = *f;
				int targ = get(name_map, target(e, *(g.graph)));   // TODO TROCADO PELO GLOBAL
				// BOOST_LOG_TRIVIAL(info) << "[updateProcessToProcessImbalanceMatrixLocal] Processing neighbor " << targ
				// 					<< " from processor " << target(e, *(g.graph)).owner;
				double weight = ew[e].weight;
				// Etapa 1: subtracao dos imbalances antigos
				long old_kj = previousSplitgraphClusterArray[targ];
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
				long new_kj = newSplitgraphClusterArray[targ];
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
		}
		// NAO EH NECESSARIO ATUALIZAR O VERTEX IMBALANCE ABAIXO, POIS A BUSCA LOCAL (VND) EH FIRST IMPROVEMENT
		// vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
	BOOST_LOG_TRIVIAL(debug) << "updateProcessToProcessImbalanceMatrixLocal done.";
}

std::vector<long> SplitgraphUtil::getListOfVeticesInCluster(ParallelBGLSignedGraph& g, const Clustering& globalClustering,
		long clusterNumber) {

	long n = g.getGlobalN();
	ClusterArray globalClusterArray = globalClustering.getClusterArray();
	std::vector<long> vertexList;
	for(long vx = 0; vx < n; vx++) {
		if(globalClusterArray[vx] == clusterNumber) {
			// BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Choosing vertex " << vx << " to be moved, it belongs to cluster" << clusterNumber << ".";
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

std::vector<Imbalance> SplitgraphUtil::calculateProcessInternalImbalance(SignedGraph& g,
		ClusterArray& splitGraphCluster, ClusterArray& globalCluster, int numberOfProcesses) {
	boost::mpi::communicator world;
	std::vector<Imbalance> processInternalImbalance(numberOfProcesses, Imbalance());
	for(int px = 1; px < numberOfProcesses; px++) {
		// worker processes, obtain answer via MPI
		// TODO merge the other matrices
		InputMessageSplitUtil imsg;
		imsg.splitGraphCluster = splitGraphCluster;
		imsg.globalCluster = globalCluster;
		imsg.numberOfProcesses = numberOfProcesses;
		imsg.functionRequested = InputMessageSplitUtil::calculateProcessInternalImbalance;
		world.send(px, MPIMessage::INPUT_MSG_SPLITUTIL_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << px << ".";
	}
	processInternalImbalance[0] = calculateProcessInternalImbalanceLocal(&g, globalCluster);
	for(int i = 1; i < numberOfProcesses; i++) {
		OutputMessageSplitUtil omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfProcesses) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ".";
		// MERGE IMBALANCE VALUES
		processInternalImbalance[procNum] = omsg.imbalance;
	}
	stringstream ss;
	ss << "Calulated InternalImbalanceVector: ";
	for(int x = 0; x < numberOfProcesses; x++) {
		ss << processInternalImbalance[x].getValue() << " ";
	}
	BOOST_LOG_TRIVIAL(info) << ss.str();
	return processInternalImbalance;
}

Imbalance SplitgraphUtil::calculateProcessInternalImbalance(ParallelBGLSignedGraph *g, Clustering& c, unsigned int processNumber) {
	ClusterArray clArray = c.getClusterArray();
	return calculateProcessInternalImbalance(g, clArray, processNumber);
}

// FUNCAO MODIFICADA PARA A BOOST PARALLEL BGL, OBTENDO O IMBALANCE DE CADA PROCESSO VIA CONSULTA POR MENSAGEM MPI
Imbalance SplitgraphUtil::calculateProcessInternalImbalance(ParallelBGLSignedGraph *g, ClusterArray& c, unsigned int processNumber) {
	boost::mpi::communicator world;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] calculateProcessInternalImbalance invoked for process " << processNumber;
	if(processNumber == 0) {  // leader process
		return calculateProcessInternalImbalanceLocal(g, c);
	} else {  // worker processes, obtain answer via MPI
		InputMessageSplitUtil imsg;
		imsg.globalCluster = c;
		imsg.functionRequested = InputMessageSplitUtil::calculateProcessInternalImbalance;
		world.send(processNumber, MPIMessage::INPUT_MSG_SPLITUTIL_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << processNumber << ".";

		OutputMessageSplitUtil omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ".";
		return omsg.imbalance;
	}
}

Imbalance SplitgraphUtil::calculateProcessInternalImbalanceLocal(ParallelBGLSignedGraph *g, ClusterArray& c) {

	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	BOOST_LOG_TRIVIAL(info) << "Obtaining name_map...";
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g->graph));

	BOOST_LOG_TRIVIAL(info) << "Obtaining global_index...";
	mpi_process_group pg = g->graph->process_group();
	boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
	  global_index(pg, num_vertices(*(g->graph)),
				   get(vertex_index, *(g->graph)), get(vertex_global, *(g->graph)));

	double positiveSum = 0.0, negativeSum = 0.0;
	int n = g->getGlobalN();
	ClusterArray myCluster = c;
	// local subgraph creation
	LocalSubgraph lsg = make_local_subgraph(*(g->graph));
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, lsg);
	LocalSubgraph::edge_descriptor e;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] calculateProcessInternalImbalanceLocal...";

	// For each vertex i that belongs to process 'processNumber'
	BGL_FORALL_VERTICES(v, lsg, LocalSubgraph) {  // For each vertex v
		int i = get(name_map, v);  // TODO trocado o local pelo global
		long ki = myCluster[i];
		LocalSubgraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(v, lsg); f != l; ++f) {
			e = *f;
			double weight = ew[e].weight;
			long j = get(name_map, target(*f, lsg));  // TODO trocado o local pelo global
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
	return Imbalance(positiveSum, negativeSum);
}

// FUNCAO MODIFICADA PARA A BOOST PARALLEL BGL, OBTENDO OS DADOS DE CADA PROCESSO VIA CONSULTA POR MENSAGEM MPI
std::vector<Coordinate> SplitgraphUtil::obtainListOfImbalancedClusters(SignedGraph& g,
		Clustering& globalClustering) {
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] obtainListOfImbalancedClusters invoked.";
	mpi::communicator world;
	int numberOfProcesses = world.size();
	std::vector<Coordinate> result;
	for(int px = 1; px < numberOfProcesses; px++) {
		// worker processes, obtain answer via MPI
		InputMessageSplitUtil imsg;
		imsg.globalClustering = globalClustering;
		imsg.functionRequested = InputMessageSplitUtil::obtainListOfImbalancedClusters;
		world.send(px, MPIMessage::INPUT_MSG_SPLITUTIL_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message sent to process " << px << ".";
	}
	// leader process
	std::vector<Coordinate> clusterList = obtainListOfImbalancedClustersLocal(g, globalClustering);
	// merge to result vector
	result.insert(result.end(), clusterList.begin(), clusterList.end());
	for(int i = 1; i < numberOfProcesses; i++) {
		OutputMessageSplitUtil omsg;
		mpi::status stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
		int tag = stat.tag();
		int procNum = stat.source();
		if(tag == MPIMessage::TERMINATE_MSG_TAG) {  // processing error
			BOOST_LOG_TRIVIAL(error) << "Error message received from process " << procNum;
			i++;
			while(i < numberOfProcesses) {  // receive the remaining messages to empty the buffer
				stat = world.recv(mpi::any_source, mpi::any_tag, omsg);
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << stat.source();
				i++;
			}
			throw std::invalid_argument( "Error message received from slave. See error logs." );
		}
		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Message received from process " << procNum << ".";
		// MERGE IMBALANCE MATRIX
		result.insert(result.end(), omsg.listOfImbalancedClusters.begin(), omsg.listOfImbalancedClusters.end());
	}
	return result;
}

std::vector<Coordinate> SplitgraphUtil::obtainListOfImbalancedClustersLocal(SignedGraph& g,
		Clustering& globalClustering) {

	typedef property_map<ParallelGraph, vertex_index_t>::type VertexIndexMap;
	typedef property_map<ParallelGraph, vertex_global_t>::type VertexGlobalMap;
	BOOST_LOG_TRIVIAL(info) << "Obtaining name_map...";
	property_map<ParallelGraph, vertex_name_t>::type name_map = get(vertex_name, *(g.graph));

	BOOST_LOG_TRIVIAL(info) << "Obtaining global_index...";
	mpi_process_group pg = g.graph->process_group();
	boost::parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
	  global_index(pg, num_vertices(*(g.graph)),
				   get(vertex_index, *(g.graph)), get(vertex_global, *(g.graph)));

	long n = g.getN();
	// local subgraph creation
	// LocalSubgraph lsg = make_local_subgraph(*(g.graph));

	ParallelGraph::edge_descriptor e;
	// Cluster to cluster matrix containing positive / negative contribution to imbalance
	long nc = globalClustering.getNumberOfClusters();
	ClusterArray globalCluster = globalClustering.getClusterArray();
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);
	matrix<double> clusterEdgeSumMatrix = zero_matrix<double>(nc, nc);
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g.graph));

	BGL_FORALL_VERTICES(v, *(g.graph), ParallelGraph) {  // For each vertex v
		int i = get(name_map, v); // v.local;  // TROCADO POR GLOBAL
		long ki = globalCluster[i];
		ParallelGraph::out_edge_iterator f, l;
		// For each out edge of i
		for (boost::tie(f, l) = out_edges(v, *(g.graph)); f != l; ++f) {
			e = *f;
			int targ = get(name_map, target(e, *(g.graph)));   // TROCADO POR GLOBAL
			double weight = ew[e].weight;
			bool sameCluster = (globalCluster[targ] == ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are in the same cluster
				clusterImbMatrix(ki, globalCluster[targ]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				clusterImbMatrix(ki, globalCluster[targ]) += fabs(weight);
			}
			clusterEdgeSumMatrix(ki, globalCluster[targ]) += fabs(weight);
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


std::vector<Coordinate> SplitgraphUtil::obtainListOfOverloadedProcesses(ParallelBGLSignedGraph& g,
		const ProcessClustering& processClustering, long maximumNumberOfVertices) {
	long n = g.getGlobalN();
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

std::vector<Coordinate> SplitgraphUtil::obtainListOfOverloadedProcesses(ParallelBGLSignedGraph& g,
		const ProcessClustering& processClustering) {

	long n = g.getGlobalN();
	const Clustering& splitgraphClustering = processClustering.getSplitgraphClustering();
	int numberOfProcesses = splitgraphClustering.getNumberOfClusters();
	// A vertex-overloaded process is a process with more than (n / numberOfProcesses) vertices.
	long numberOfEquallyDividedVertices = (long)ceil(n / (double)numberOfProcesses);
	return obtainListOfOverloadedProcesses(g, processClustering, numberOfEquallyDividedVertices);
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
	LocalSubgraph::edge_descriptor e;
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, *(g.graph));
	double maxImbalance = -DBL_MAX;
	long maxVertex = 0;

	// TODO VERIFICAR A MELHOR FORMA DE REPROGRAMAR A LEITURA DE ARESTAS COM A PARALLEL BGL
	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = globalCluster[i];
		LocalSubgraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// only processes vertexes belonging to the process pair
		if((splitGraphCluster[i] == processPair.x) or (splitGraphCluster[i] == processPair.y)) {
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(vertex(i, *(g.graph)), *(g.graph)); f != l; ++f) {
				e = *f;
				Vertex src = source(e, *(g.graph)).local, targ = target(e, *(g.graph)).local;
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

	// local subgraph creation
	LocalSubgraph lsg = make_local_subgraph(*(g->graph));
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, lsg);
	LocalSubgraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCount = 0, positiveEdgeCount = 0, negativeEdgeCount = 0;

	// validates the edges from vertex u to cliqueC
	LocalSubgraph::vertex_descriptor vx_u;
	BGL_FORALL_VERTICES(vx, lsg, LocalSubgraph) {  // For each vertex u
		if(vx.local == u) {
			vx_u = vx;
			break;
		}
	}
	LocalSubgraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(vx_u, lsg); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, lsg).local;
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

	// local subgraph creation
	LocalSubgraph lsg = make_local_subgraph(*(g->graph));
	boost::property_map<ParallelGraph, edge_properties_t>::type ew = boost::get(edge_properties, lsg);
	LocalSubgraph::edge_descriptor e;
	// every vertex in cliqueC must be connected to vertex u
	// counts the number of clique elements connected to u
	long totalEdgeCountU = 0, positiveEdgeCountU = 0, negativeEdgeCountU = 0;
	long totalEdgeCountV = 0, positiveEdgeCountV = 0, negativeEdgeCountV = 0;

	// validates the edges from vertex u to cliqueC
	LocalSubgraph::vertex_descriptor vx_u;
	BGL_FORALL_VERTICES(vx, lsg, LocalSubgraph) {  // For each vertex u
		if(vx.local == u) {
			vx_u = vx;
			break;
		}
	}
	LocalSubgraph::out_edge_iterator f2, l2;
	for (boost::tie(f2, l2) = out_edges(vx_u, lsg); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, lsg).local;
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
	LocalSubgraph::vertex_descriptor vx_v;
	BGL_FORALL_VERTICES(vx, lsg, LocalSubgraph) {  // For each vertex v
		if(vx.local == v) {
			vx_v = vx;
			break;
		}
	}
	for (boost::tie(f2, l2) = out_edges(vx_v, lsg); f2 != l2; ++f2) {
		e = *f2;
		long j = target(*f2, lsg).local;
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


void SplitgraphUtil::validaSplitgraphArray(SignedGraph &g, ProcessClustering& processClustering, Clustering& globalClustering) {
	/*  DISABLED - USE ONLY FOR DEBUGGING PURPOSES
	const ClusterArray splitgraphClusterArray = processClustering.getClusterArray();
	const std::vector<unsigned int> clusterProcessOrigin = globalClustering.getClusterProcessOrigin();
	long nc = globalClustering.getNumberOfClusters();
	long n = g.getN();
	if(nc != clusterProcessOrigin.size()) {
		BOOST_LOG_TRIVIAL(error) << "nc != clusterProcessOrigin.size()";
	}

	for(long k = 0; k < nc; k++) {
		std::vector<long> vertexList = this->getListOfVeticesInCluster(g, globalClustering, k);
		unsigned int processOrigin = clusterProcessOrigin[k];
		for(std::vector<long>::iterator it = vertexList.begin(); it != vertexList.end(); it++) {
			long v = *it;
			if(processOrigin != splitgraphClusterArray[v]) {
				BOOST_LOG_TRIVIAL(error) << "Vertex " << v << " (from global cluster " << k << ") is in process " << splitgraphClusterArray[v]
					<< " but its cluster is in process " << processOrigin;
			}
		}
	} */
}

} /* namespace ils */
} /* namespace resolution */
