/*
 * ImbalanceSubgraphParallelILS.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: mlevorato
 */

#include "./include/ImbalanceSubgraphParallelILS.h"
#include "util/include/ProcessUtil.h"
#include "util/include/MPIMessage.h"
#include "util/parallel/include/MPIUtil.h"
#include "../include/CUDAILS.h"
#include "problem/include/RCCProblem.h"

#include "../../construction/include/GainFunctionFactory.h"
#include "../../construction/include/GainFunction.h"
#include "graph/include/NeighborhoodSearchFactory.h"
#include "graph/include/ParallelNeighborhoodSearch.h"
#include "graph/include/SequentialNeighborhoodSearch.h"

#include <iomanip>
#include <cstring>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <cfloat>
#include <algorithm>

// TODO CHANGE TO APPLICATION CLI PARAMETER
#define PERCENTAGE_OF_PROCESS_PAIRS 0.2
#define PERCENTAGE_OF_MOST_IMBALANCED_VERTICES_TO_BE_MOVED 0.25

namespace ublas = boost::numeric::ublas::detail;

namespace resolution {
namespace ils {

using namespace resolution::construction;
using namespace boost::algorithm;
using namespace std;
using namespace util;
using namespace util::parallel;

ImbalanceSubgraphParallelILS::ImbalanceSubgraphParallelILS(const int& allocationStrategy, const int& slaves, const int& searchSlaves,
		const bool& split, const bool& cuda) : ILS(),
		machineProcessAllocationStrategy(allocationStrategy), numberOfSlaves(slaves), numberOfSearchSlaves(searchSlaves),
		splitGraph(split), cudaEnabled(cuda), vertexImbalance(), timeSpentAtIteration() {

}

ImbalanceSubgraphParallelILS::~ImbalanceSubgraphParallelILS() {
	// TODO Auto-generated destructor stub
}

Clustering ImbalanceSubgraphParallelILS::executeILS(ConstructClustering *construct, VariableNeighborhoodDescent *vnd,
		SignedGraph *g, const int& iter, const int& iterMaxILS, const int& perturbationLevelMax,
		ClusteringProblem& problem, ExecutionInfo& info) {

	stringstream constructivePhaseResults;
	stringstream iterationResults;
	stringstream iterationTimeSpent;
	unsigned int numberOfProcesses = numberOfSlaves + 1;

	// *** STEP A => PREPROCESSING PHASE OF DISTRIBUTED METAHEURISTIC: New split graph partitioning here
	// Triggers local processing time calculation
	boost::timer::cpu_timer timer;
	timer.start();
	boost::timer::cpu_times start_time = timer.elapsed();

	// TODO PARAMETRIZAR DIVISAR POR VERTICE OU POR SOMA DE ARESTAS
	Clustering splitgraphClustering = preProcessSplitgraphPartitioning(g, problem, true);
	BOOST_LOG_TRIVIAL(info) << "Initial split graph processor partitioning I(P) = " << splitgraphClustering.getImbalance().getValue();

	// *** STEP B => Calls the individual ILS processing for each subgraph, invoking Parallel ILS with MPI
	// WARNING: THE SPLITGRAPHCLUSTERARRAY IS A SEPARATE DATA STRUCTURE, DIFFERENT THAN THE CURRENT CLUSTERING
	// IT CONTROLS THE PARTITIONING BETWEEN PROCESSORS (SPLIT GRAPH)
	ClusterArray initialSplitgraphClusterArray = splitgraphClustering.getClusterArray();
	matrix<double> processClusterImbMatrix = calculateProcessToProcessImbalanceMatrix(*g, initialSplitgraphClusterArray);
	Clustering Cc = distributeSubgraphsBetweenProcessesAndRunILS(construct,
			vnd, g, iter, iterMaxILS, perturbationLevelMax, problem, info, initialSplitgraphClusterArray, processClusterImbMatrix);

	timer.stop();
	boost::timer::cpu_times end_time = timer.elapsed();
	double timeSpentInConstruction = (end_time.wall - start_time.wall) / double(1000000000);
	timeSpentInILS = timeSpentInConstruction;
	timer.resume();
	start_time = timer.elapsed();

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Initial solution found: I(P) = " <<
			Cc.getImbalance().getValue();
	constructivePhaseResults << 1 << "," << Cc.getImbalance().getValue() << ","
					<< Cc.getImbalance().getPositiveValue()
					<< "," << Cc.getImbalance().getNegativeValue() << "," << Cc.getNumberOfClusters()
					<< "," << fixed << setprecision(4) << timeSpentInConstruction << "\n";
	constructivePhaseResults << "Average initial I(P)," << fixed << setprecision(4) << Cc.getImbalance().getValue()
						<< "\n";
	Cc.printClustering(g->getN());
	Clustering bestClustering(Cc);
	Clustering bestSplitgraphClustering(splitgraphClustering);

	// STEP 2: TRY TO IMPROVE THE GLOBAL SOLUTION THROUGH A DISTRIBUTED METAHEURISTIC WHICH
	// EXCHANGES VERTICES BETWEEN PROCESSES
	int iterCount = 0;
	bool foundBetterSolution = false;
	do {
		// 1. O processo mestre será o responsável por manter duas estruturas de dados de controle do imbalance:
		//		Um vetor com os vértices que mais contribuem para I(P) (soma de imbalance de cada vertice)
		//		Uma matriz com a soma do imbalance entre processos
		// => Calculates the imbalance info between process partitions and of each vertex individually
		ClusterArray bestSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
		foundBetterSolution = false;

		// 2. O processo mestre elege dois processos condidatos a participar de uma movimentação de vértices
		// a) O processo mestre escolhe os dois processos com maior valor de imbalance entre si (usando a matriz 1b)
		// Chooses the top 20% biggest elements in the matrix
		// Coordinate processPair = findMaximumElementInMatrix(processClusterImbMatrix);
		std::vector< Coordinate > imbMatrixElementList = getMatrixElementsAsList(processClusterImbMatrix);
		int numberOfProcessCombinations = (int)floor((numberOfProcesses * (numberOfProcesses - 1)) / double(2.0));  // symmetric matrix
		int desiredNumberOfProcessPairs = (int)ceil(numberOfProcessCombinations * PERCENTAGE_OF_PROCESS_PAIRS);
		// obtain and remove the next biggest element from the list
		for(int processPairCount = 0; processPairCount < desiredNumberOfProcessPairs; processPairCount++) {
			std::vector< Coordinate >::iterator pos = std::max_element(imbMatrixElementList.begin(),
					imbMatrixElementList.end(), coordinate_ordering_asc());
			Coordinate processPair = *pos;
			imbMatrixElementList.erase(pos);

			BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] The next process pair with biggest imbalance is (" << processPair.x
					<< ", " << processPair.y << ") with imbalance = " << processClusterImbMatrix(processPair.x, processPair.y);
			// TODO Implementar aqui o VND com 2 ou mais estruturas de vizinhanca
			// TODO Alem do move 1-opt, implementar o move de um cluster inteiro

			// ************************   Move   1 - opt  ***********************************
			// b) O processo mestre escolhe um vértice a ser movimentado (move) de um processo para o outro (usando o vetor 1 A)
			// select a vertex v to move - v is the vertex that causes the biggest imbalance
			// TODO the biggest imbalance in the whole graph or only between the pair of processes x and y?
			// for now it is the vertex v that belongs to the processPair and causes the biggest imbalance in the whole graph
			// Tries to move the top 25% of most imbalanced vertices
			long numberOfVerticesInProcessPair = bestSplitgraphClustering.getClusterSize(processPair.x)
					+ bestSplitgraphClustering.getClusterSize(processPair.y);
			long numberOfImbalancedVerticesToBeMoved = (int)ceil(numberOfVerticesInProcessPair * PERCENTAGE_OF_MOST_IMBALANCED_VERTICES_TO_BE_MOVED);
			ClusterArray splitgraphClusterList = bestSplitgraphClustering.getClusterArray();
			for(long niv = 0; niv < numberOfImbalancedVerticesToBeMoved; niv++) {
				long v = findMostImbalancedVertexInProcessPair(*g, splitgraphClusterList, processPair);
				int currentProcess = splitgraphClusterList[v];
				// invalidates vertex v from the cluster array (== removes from the list), so that it is not chosen again
				splitgraphClusterList[v] = Clustering::NO_CLUSTER;
				BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Choosing vertex " << v << " as the next most imbalanced one.";

				// TODO PENSAR EM UMA HEURISTICA PARA ESCOLHER O PROCESSO DE DESTINO DO VERTICE V COM MAIOR IMBALANCE
				// OU ENTAO TENTAR MOVER O VERTICE V PARA QUALQUER OUTRO PROCESSO (ITERAR)
				ClusterArray currentSplitgraphClusterArray = bestSplitgraphClustering.getClusterArray();
				ClusterArray tempSplitgraphClusterArray = currentSplitgraphClusterArray;
				for(int px = 0; px < numberOfProcesses; px++) {
					if(px != currentProcess) {
						// Realiza a movimentacao de um vertice (1-opt move)
						tempSplitgraphClusterArray[v] = px;
						std::vector<long> listOfModifiedVertices;
						listOfModifiedVertices.push_back(v);
						// TODO RECALCULAR A MATRIZ ABAIXO DE FORMA INCREMENTAL, reutilizando a matrix processClusterImbMatrix
						ublas::matrix<double> tempProcessClusterImbMatrix2 = calculateProcessToProcessImbalanceMatrix(*g, tempSplitgraphClusterArray);
						ublas::matrix<double> tempProcessClusterImbMatrix = processClusterImbMatrix;
						// Realiza o calculo do delta da matriz de imbalance entre processos
						updateProcessToProcessImbalanceMatrix(*g, currentSplitgraphClusterArray, tempSplitgraphClusterArray,
								listOfModifiedVertices, tempProcessClusterImbMatrix);
						// Valida se a matriz incremental e a full sao iguais
						bool igual = ublas::equals(tempProcessClusterImbMatrix, tempProcessClusterImbMatrix2, 0.001, 0.1);
						BOOST_LOG_TRIVIAL(info) << "*** As matrizes de delta sao iguais: " << igual;

						// c) O processo mestre solicita que cada processo participante execute um novo ILS tendo como base a configuracao de cluster modificada
						Clustering newGlobalClustering = distributeSubgraphsBetweenProcessesAndRunILS(construct, vnd, g,
									iter, iterMaxILS, perturbationLevelMax, problem, info, tempSplitgraphClusterArray, tempProcessClusterImbMatrix);
						BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] New solution found: I(P) = " << newGlobalClustering.getImbalance().getValue();
						if(newGlobalClustering.getImbalance().getValue() < bestClustering.getImbalance().getValue()) {
							BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Improved solution found! I(P) = " << newGlobalClustering.getImbalance().getValue();
							bestClustering = newGlobalClustering;

							// Atualiza a matriz de imbalance entre processos apos aceitar uma nova solucao
							processClusterImbMatrix = tempProcessClusterImbMatrix;

							bestSplitgraphClustering = Clustering(tempSplitgraphClusterArray, *g, problem);
							foundBetterSolution = true;
							break;
						}
					}
				}
				if(foundBetterSolution) {
					break;
				}
				// d) Dúvida: rodar um novo ILS inteiro para cada movimentação ou apenas uma busca local?
				// e) Cada processo envia o resultado (novo clustering) para o mestre e este calcula se vale a pena manter a alteração
				// f) O cálculo feito pelo mestre baseia-se na variação do imbalance dentro de cada processo (delta I(P) int) mais o novo imbalance externo, entre os processos (I(P) ext)
				// 		(a variação dentro de cada processo já é fornecida por cada processo, via mensagem – diferença do I(P) anterior)
				// 		(a variação do imbalance entre os processos é calculada varrendo a vizinhança do vértice movido)
				// g) Se gerar um imbalance pior, o movimento pode ser rejeitado, gerando uma volta da configuração anterior, por isso cada processo deve guardar o clustering anterior
			}
			if(foundBetterSolution) {
				break;
			}
		}
		// Last step. Stops the timer and stores the elapsed time
		timer.stop();
		end_time = timer.elapsed();

		// Write the results into ostream os, using csv format
		// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
		timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
		// iteration results are only being shown when the solution improved
		iterationResults << (iterCount+1) << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
				<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
				<< "," << fixed << setprecision(4) << timeSpentInILS << "\n";
		timer.resume();
		start_time = timer.elapsed();
		// if elapsed time is bigger than timeLimit, break
		if(timeSpentInILS >= vnd->getTimeLimit()) {
			BOOST_LOG_TRIVIAL(info) << "Time limit exceeded." << endl;
			break;
		}
		iterCount++;
	} while(foundBetterSolution);

	// prints all the info about time spent in each distributed MH invocation
	for(int it = 0; it < timeSpentAtIteration.size(); it++) {
		for(int px = 0; px < numberOfProcesses; px++) {
			iterationTimeSpent << timeSpentAtIteration[it][px] << ", ";
		}
		iterationTimeSpent << "\n";
	}

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Number of improvements to initial solution: " << iterCount;
	// *** ONLY IN THE END OF THE DISTRIBUTED METAHEURISTIC: runs a ILS iteration (or a local search?) over the merged global solution

	// 6. Runs a local search over the merged global solution, trying to improve it
	// TODO Run a full ILS iteration over the split graph initial solution
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Applying local search over the initial solution.";
	NeighborhoodSearch* neigborhoodSearch;
	GainFunctionFactory functionFactory(g);
	NeighborhoodSearchFactory nsFactory(machineProcessAllocationStrategy, this->numberOfSlaves, this->numberOfSearchSlaves);
	neigborhoodSearch = &nsFactory.build(NeighborhoodSearchFactory::SEQUENTIAL);
	VariableNeighborhoodDescent vnd2(*neigborhoodSearch, vnd->getRandomSeed(), vnd->getNeighborhoodSize(),
			false, vnd->getTimeLimit());
	bestClustering = vnd2.localSearch(g, bestClustering, 0, problem, timeSpentInILS, info.processRank);
	CBest = bestClustering;
	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after applying local search: I(P) = " << bestClustering.getImbalance().getValue();

	ConstructClustering NoConstruct(functionFactory.build(construct->getGainFunctionType()), vnd->getRandomSeed(),
			construct->getAlpha(), &bestClustering);
	CUDAILS cudails;
	bestClustering = cudails.executeILS(&NoConstruct, vnd, g, iter, 1, perturbationLevelMax, problem, info);

	// 7. Stops the timer and stores the elapsed time
	timer.stop();
	end_time = timer.elapsed();

	// 8. Write the results into ostream os, using csv format
	// Format: iterationNumber,imbalance,imbalance+,imbalance-,time(s),boolean
	timeSpentInILS += (end_time.wall - start_time.wall) / double(1000000000);
	iterationResults << (iterCount+2) << "," << bestClustering.getImbalance().getValue() << "," << bestClustering.getImbalance().getPositiveValue()
			<< "," << bestClustering.getImbalance().getNegativeValue() << "," << bestClustering.getNumberOfClusters()
			<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) << "\n";
	Imbalance bestValue = bestClustering.getImbalance();
	int iterationValue = 1;
	iterationResults << "Best value," << fixed << setprecision(4) << bestValue.getValue()
				<< "," << bestValue.getPositiveValue()
				<< "," << bestValue.getNegativeValue()
				<< setprecision(0)
				<< "," << bestClustering.getNumberOfClusters()
				<< "," << (iterationValue+1)
				<< "," << fixed << setprecision(4) << (timeSpentInILS + timeSpentInConstruction) // timeSpentOnBestSolution
				<< "," << iter
				<< "," << numberOfTestedCombinations << "\n";

	BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Solution found after 1 ILS parallel graph iteration: I(P) = " << bestClustering.getImbalance().getValue();

	// Saves the iteration results to csv file
	generateOutputFile(problem, iterationResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("iterations-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Saves the initial solutions data to file
	generateOutputFile(problem, constructivePhaseResults, info.outputFolder, info.fileId, info.executionId, info.processRank, string("initialSolutions-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);
	// Save the time spent at each iteration by each process to csv file
	generateOutputFile(problem, iterationTimeSpent, info.outputFolder, info.fileId, info.executionId, info.processRank, string("timespent-splitgraph"), construct->getAlpha(), vnd->getNeighborhoodSize(), iter);

	return bestClustering;
}

Clustering ImbalanceSubgraphParallelILS::preProcessSplitgraphPartitioning(SignedGraph *g, ClusteringProblem& problem, bool partitionByVertex) {
	// 1. Divide the graph into (numberOfSlaves + 1) non-overlapping parts

	unsigned int numberOfProcesses = numberOfSlaves + 1;
	string strNP = boost::lexical_cast<std::string>(numberOfProcesses);
	BOOST_LOG_TRIVIAL(info) << "Invoking partitioning for spit graph (numberOfProcesses = " << strNP << ")...";

	/**  ETAPA DE PRE-PROCESSAMENTO - GERACAO DA PARTICAO INICIAL DO GRAFO:
	 n: numero de vertices
	 p: numero de processos
	P0: Gr = G (para o processo 0);
	Loop sequencial, cada processador executa o loop, um processador pi por vez:
	{
		Si = vazio;
		Escolha um nó com menor cadinalidade de arestas negativas entre do grafo residual Gr (o que vai aumentar menos o imbalance)
		(opção 2: maior cardinalidade de arestas positivas)

		Enquanto a cardinalidade do conjunto Si for menor que o n / p:
		{
			Escolha um nó ni pert (Gr \ Si) com maior cardinalidade de arestas positivas entre ni e Si;
			Si = Si + ni;
		}
		Gr = Gr – Si
		Envia Gr para o processo pi+1
	}
	Rodar o ILS dentro de cada subgrafo Gi;
	Merge das soluçoes;
	Busca 1 iteração do ILS sobre a solucao global;  */

	// *** INITIALIZATION: Gr = G (in the begginning, the residual graph equals the whole graph);
	long n = g->getN();
	long m = g->getM();
	// initializes a vector containing vertex degree info (vertex number and its pos and neg degrees)
	std::vector<VertexDegree> Gr;
	for(long i = 0; i < n; i++) {
		Gr.push_back(VertexDegree(i, 0.0, 0.0));
	}
	long desiredCardinality = 0;
	// a list containing the partition number each vertex belongs to
	ClusterArray splitgraphClusterArray(n);
	if(partitionByVertex) {
		desiredCardinality = long(floor(n / (double)numberOfProcesses));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " vertices.";
	} else {
		desiredCardinality = long(floor(m / (double)numberOfProcesses));
		BOOST_LOG_TRIVIAL(info) << "Desired cardinality of each partition is " << desiredCardinality << " edges.";
	}

	// For each process pi
	for(int pi = 0; pi < numberOfProcesses; pi++) {
		BOOST_LOG_TRIVIAL(debug) << "Processing partition for processor " << pi;
		// *** STEP 1: Chooses a node with the smallest negative cardinality INSIDE the residual graph Gr
		//  (will probably result in the smallest increase in the imbalance)
		Clustering GrCluster(*g);  // uses Clustering data structure to control the elements in the set
		assert(Gr.size() > 0);
		GrCluster.addCluster(*g, problem, Gr[0].id, false);
		for(int i = 1; i < Gr.size(); i++) {
			GrCluster.addNodeToCluster(*g, problem, Gr[i].id, 0, false);
		}
		BOOST_LOG_TRIVIAL(debug) << "Created initial residual graph list of size " << GrCluster.getClusterSize(0);

		// updates the negative edge sum of each vertex ni in Gr (neg edge sum between ni and Si)
		for(int i = 0; i < Gr.size(); i++) {
			long ni = Gr[i].id;
			Gr[i].negativeDegree = fabs(g->getNegativeEdgeSumBetweenVertexAndClustering(ni, GrCluster.getClusterArray()));
		}
		// Obtains the vertex in Gr which has the SMALLEST absolute negative edge sum
		std::vector<VertexDegree>::iterator it_min = std::min_element(Gr.begin(), Gr.end(), neg_degree_ordering_asc());  // O(n)
		long vertex_y = it_min->id;
		BOOST_LOG_TRIVIAL(debug) << "The vertex in Gr with the smallest neg edge sum is " << vertex_y << " with sum = " << it_min->negativeDegree;
		// the Si clustering will have only 1 cluster, cluster zero
		Clustering Si(*g);  // Si = empty, uses Clustering data structure to control the elements in the set
		Si.addCluster(*g, problem, vertex_y, false);
		Gr.erase(it_min);

		ClusterArray cTemp = Si.getClusterArray();
		long numberOfEdges = 0;
		while(Gr.size() > 0) {
			// updates the positive edge sum of each vertex ni in Gr (pos edge sum between ni and Si)
			for(int i = 0; i < Gr.size(); i++) {
				long ni = Gr[i].id;
				Gr[i].positiveDegree = g->getPositiveEdgeSumBetweenVertexAndClustering(ni, Si.getClusterArray());
			}
			std::vector<VertexDegree>::iterator it_max = std::max_element(Gr.begin(), Gr.end(), pos_degree_ordering_asc());  // O(n)
			BOOST_LOG_TRIVIAL(debug) << "ClusterSize = " << Si.getClusterSize(0) << ": The vertex in (Gr - Si) with the biggest pos edge sum is " << it_max->id << " with sum = " << it_max->positiveDegree;

			// Chooses a node ni belonging to (Gr \ Si) with the biggest cardinality of positive edges between ni and Si;
			long ni = it_max->id;
			// Si = Si + ni
			Si.addNodeToCluster(*g, problem, ni, 0, false);
			// Gr = Gr - ni
			Gr.erase(it_max);
			numberOfEdges += g->getOutDegree(ni);

			if(partitionByVertex) {
				if(Si.getClusterSize(0) >= desiredCardinality)  break;
			} else {
				if(numberOfEdges >= desiredCardinality)  break;
			}
		}
		BOOST_LOG_TRIVIAL(debug) << "Partition ready.";
		// Gr = Gr – Si  TODO: subtrair a lista partialGraph da lista residualGraph -> ja feito acima
		// Envia Si para o processo pi+1
		ClusterArray SiArray = Si.getClusterArray();
		for(int i = 0; i < n; i++) {
			if(SiArray[i] >= 0) {
				splitgraphClusterArray[i] = pi;
			}
		}
	}

	// processes the initial solution generated by the split graph partitioning
	BOOST_LOG_TRIVIAL(info) << "Generated a cluster array of " << splitgraphClusterArray.size() << " vertices.";
	Clustering Cc(splitgraphClusterArray, *g, problem);

	return Cc;
}

Clustering ImbalanceSubgraphParallelILS::distributeSubgraphsBetweenProcessesAndRunILS(ConstructClustering *construct,
		VariableNeighborhoodDescent *vnd, SignedGraph *g, const int& iter, const int& iterMaxILS,
		const int& perturbationLevelMax, ClusteringProblem& problem, ExecutionInfo& info,
		ClusterArray& splitgraphClusterArray, matrix<double>& processClusterImbMatrix) {

	// Creates the subgraphs for each processor, based on the splitgraphClusterArray structure
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (i = 0) does part of the work too
	std::vector<int> slaveList;
	MPIUtil::populateListOfMasters(machineProcessAllocationStrategy, slaveList, info.processRank, numberOfSlaves, numberOfSearchSlaves);

	long n = g->getN();
	unsigned int numberOfProcesses = numberOfSlaves + 1;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
		if(k < 0) {  // reuses CC problem's best solution in the constructive phase
			CCclustering = *(construct->getCCclustering());
		}
	}

	std::vector< std::vector< long > > verticesInCluster(numberOfProcesses, std::vector< long >());
	for(long i = 0; i < n; i++) {
		long k = splitgraphClusterArray[i];
		verticesInCluster[k].push_back(i);
	}

	// Creates numberOfProcesses subgraphs
	std::vector<SubGraph> subgraphList;
	// each subgraph will have a subset of the main graph's nodes and edges, based on the previous clustering
	for(int p = 0; p < numberOfProcesses; p++) {
		SubGraph sg = (g->graph).create_subgraph(); //verticesInCluster[k].begin(), verticesInCluster[k].end());
		subgraphList.push_back(sg);  // --> SUBGRAPH COPY CTOR NOT WORKING!!!
		for(std::vector<long>::iterator it = verticesInCluster[p].begin(); it != verticesInCluster[p].end(); it++) {
			add_vertex(*it, subgraphList.back());
			// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << p;
		}

		BOOST_LOG_TRIVIAL(info) << "[Parallel ILS SplitGraph] Subgraph " << p << ": num_edges = " <<
				num_edges(subgraphList.back()) << " , num_vertices = " << num_vertices(subgraphList.back());
	}

	for(int i = 0; i < numberOfSlaves; i++) {
		// 2. Distribute numberOfSlaves graph parts between the ILS Slave processes
		InputMessageParallelILS imsg(g->getId(), g->getGraphFileLocation(), iter, construct->getAlpha(), vnd->getNeighborhoodSize(),
							problem.getType(), construct->getGainFunctionType(), info.executionId, info.fileId, info.outputFolder, vnd->getTimeLimit(),
							numberOfSlaves, numberOfSearchSlaves, vnd->isFirstImprovementOnOneNeig(), iterMaxILS, perturbationLevelMax, k);
		if(k < 0) {
			imsg.setClustering(&CCclustering);
		}
		imsg.setVertexList(verticesInCluster[i+1]);
		world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_ILS_TAG, imsg);
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Message sent to process " << slaveList[i];
	}
	// 2.1. the leader does its part of the work: runs ILS using the first part of the divided graph
	CUDAILS cudails;
	SignedGraph sg(g->graph, verticesInCluster[0]);
	// sg.graph = (g->graph).create_subgraph(verticesInCluster[0].begin(), verticesInCluster[0].end());

	/*
	for(std::vector<long>::iterator it = verticesInCluster[0].begin(); it != verticesInCluster[0].end(); it++) {
		add_vertex(*it, sg.graph);
		// BOOST_LOG_TRIVIAL(info) << "Inserting vertex " << *it << " in cluster " << k;
	} */
	BOOST_LOG_TRIVIAL(debug) << "Processing subgraph with n =  " << num_vertices(sg.graph) << ", " << "e =  " << num_edges(sg.graph);

	// rebuilds construct clustering objects based on partial graph 'sg'
	GainFunctionFactory functionFactory(&sg);
	ConstructClustering defaultConstruct(functionFactory.build(construct->getGainFunctionType()),
			construct->getRandomSeed(), construct->getAlpha());
	ConstructClustering noConstruct(functionFactory.build(construct->getGainFunctionType()),
			construct->getRandomSeed(), construct->getAlpha(), construct->getCCclustering());
	ConstructClustering* construct2 = &defaultConstruct;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		int RCCk = rp.getK();
		if(RCCk < 0) {
			construct2 = &noConstruct;
		}
	}
	// Leader processing his own partition
	Clustering leaderClustering = cudails.executeILS(construct2, vnd, &sg, iter, iterMaxILS, perturbationLevelMax, problem, info);
	// Structure containing processing time of each MH slave
	std::vector<double> timeSpent(numberOfProcesses, double(0.0));
	timeSpent[0] = cudails.getTotalTimeSpent();

	// 3. the leader receives the processing results
	OutputMessage omsg;
	BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Waiting for slaves return messages.";
	// Global cluster array
	ClusterArray globalClusterArray(g->getN(), 0);
	long clusterOffset = 0;
	double internalImbalanceSum = 0.0;
	for(int i = 0; i < numberOfSlaves; i++) {
		mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_ILS_TAG, omsg);
		BOOST_LOG_TRIVIAL(debug) << "[Parallel ILS SplitGraph] Message received from process " << stat.source() << ". Obj = " <<
				omsg.clustering.getImbalance().getValue();
		// process the result of the execution of process i
		// sums the total number of tested combinations
		numberOfTestedCombinations += omsg.numberOfTestedCombinations;
		// stores the time spent by each process
		timeSpent[i+1] = omsg.timeSpent;
		// 4. Merge the partial solutions into a global solution for the whole graph
		ClusterArray localClusterArray = omsg.clustering.getClusterArray();
		internalImbalanceSum += omsg.clustering.getImbalance().getValue();
		assert(localClusterArray.size() == omsg.globalVertexId.size());
		for(long v = 0; v < localClusterArray.size(); v++) {
			// Obtains vertex v's number in the global graph
			long vglobal = omsg.globalVertexId[v];
			globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
			/*
			BOOST_LOG_TRIVIAL(info) << "Vertex " << v << " in local solution becomes vertex " << vglobal <<
					" in global solution and belongs to cluster " << localClusterArray[v] <<
					", that is, global cluster " << globalClusterArray[vglobal]; */
		}
		long nc = omsg.clustering.getNumberOfClusters();
		clusterOffset += nc;
	}
	// includes the leader's processing result as well
	// builds a global cluster array, containing each vertex'es true id in the global / full parent graph
	std::pair< graph_traits<SubGraph>::vertex_iterator, graph_traits<SubGraph>::vertex_iterator > v_it = vertices(sg.graph);
	std::vector<long> globalVertexId;
	for(graph_traits<SubGraph>::vertex_iterator it = v_it.first; it != v_it.second; it++) {
		globalVertexId.push_back(sg.graph.local_to_global(*it));
	}
	// 4. Merges the leader's partial solution into a global solution for the whole graph
	ClusterArray localClusterArray = leaderClustering.getClusterArray();
	for(long v = 0; v < localClusterArray.size(); v++) {
		// Obtains vertex v's number in the global graph
		long vglobal = globalVertexId[v];
		globalClusterArray[vglobal] = clusterOffset + localClusterArray[v];
		/*
		BOOST_LOG_TRIVIAL(info) << "Vertex " << v << " in local solution becomes vertex " << vglobal <<
								" in global solution and belongs to cluster " << localClusterArray[v] <<
								", that is, global cluster " << globalClusterArray[vglobal]; */
	}
	long nc = leaderClustering.getNumberOfClusters();
	clusterOffset += nc;
	internalImbalanceSum += leaderClustering.getImbalance().getValue();

	// Calculates the external imbalance sum (between processes)
	double externalImbalanceSum = 0.0;
	for (unsigned i = 0; i < processClusterImbMatrix.size1(); ++i) {
		for (unsigned j = 0; j < processClusterImbMatrix.size2(); ++j) {
			if(i != j) {
				externalImbalanceSum += processClusterImbMatrix(i, j);
			}
		}
	}

	// 5. Builds the clustering with the merge of each process result
	// // TODO EVITAR O RECALCULO DA FUNCAO OBJETIVO NA LINHA ABAIXO
	Clustering globalClustering(globalClusterArray, *g, problem);

	BOOST_LOG_TRIVIAL(info) << "Separate internal imbalance I(P) = " << (internalImbalanceSum) << " and external I(P) = " << (externalImbalanceSum);
	BOOST_LOG_TRIVIAL(info) << "Separate clustering imbalance I(P) = " << (internalImbalanceSum + externalImbalanceSum)
			<< " versus full calculated I(P) = " << globalClustering.getImbalance().getValue();

	// 6. Stores the time spent in this iteration of distributed MH
	timeSpentAtIteration.push_back(timeSpent);

	return globalClustering;
}

matrix<double> ImbalanceSubgraphParallelILS::calculateProcessToProcessImbalanceMatrix(SignedGraph& g, ClusterArray& myCluster) {

	long nc = numberOfSlaves + 1;
	long n = g.getN();
	DirectedGraph::edge_descriptor e;
	// Process to process matrix containing positive / negative contribution to imbalance
	matrix<double> clusterImbMatrix = zero_matrix<double>(nc, nc);
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	// Vector containing each vertex contribution to imbalance: vertexImbalance

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		DirectedGraph::out_edge_iterator f, l;
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
				clusterImbMatrix(ki, myCluster[targ.id]) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				positiveSum += weight;
				clusterImbMatrix(ki, myCluster[targ.id]) += fabs(weight);
			}
		}
		vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
	return clusterImbMatrix;
}

void ImbalanceSubgraphParallelILS::updateProcessToProcessImbalanceMatrix(SignedGraph& g,
			const ClusterArray& previousSplitgraphClusterArray,
			const ClusterArray& newSplitgraphClusterArray, const std::vector<long>& listOfModifiedVertices,
			matrix<double>& processClusterImbMatrix) {

	long nc = numberOfSlaves + 1;
	long n = g.getN();
	DirectedGraph::edge_descriptor e;
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);

	// For each vertex i in listOfModifiedVertices
	for(long item = 0; item < listOfModifiedVertices.size(); item++) {
		long i = listOfModifiedVertices[item];
		long old_ki = previousSplitgraphClusterArray[i];
		long new_ki = newSplitgraphClusterArray[i];
		DirectedGraph::out_edge_iterator f, l;
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
				processClusterImbMatrix(old_ki, old_kj) -= fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j were NOT in the same cluster
				processClusterImbMatrix(old_ki, old_kj) -= fabs(weight);
			}
			// Etapa 2: acrescimo dos novos imbalances
			long new_kj = newSplitgraphClusterArray[targ.id];
			sameCluster = (new_kj == new_ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are now in the same cluster
				processClusterImbMatrix(new_ki, new_kj) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				processClusterImbMatrix(new_ki, new_kj) += fabs(weight);
			}
		}
		// For each in edge of i
		DirectedGraph::in_edge_iterator f2, l2;
		for (boost::tie(f2, l2) = in_edges(i, g.graph); f2 != l2; ++f2) {  // in edges of i
			e= *f2;
			double weight = ew[e].weight;
			int j = source(*f2, g.graph);
			// Etapa 1: subtracao dos imbalances antigos
			long old_kj = previousSplitgraphClusterArray[j];
			bool sameCluster = (old_kj == old_ki);
			// TODO VERIFICAR SE DEVE SER RECALCULADO TAMBEM O PAR OPOSTO DA MATRIZ: (KJ, KI)
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j were in the same cluster
				processClusterImbMatrix(old_kj, old_ki) -= fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j were NOT in the same cluster
				processClusterImbMatrix(old_kj, old_ki) -= fabs(weight);
			}
			// Etapa 2: acrescimo dos novos imbalances
			long new_kj = newSplitgraphClusterArray[j];
			sameCluster = (new_kj == new_ki);
			if(weight < 0 and sameCluster) {  // negative edge
				// i and j are now in the same cluster
				processClusterImbMatrix(new_kj, new_ki) += fabs(weight);
			} else if(weight > 0 and (not sameCluster)) {  // positive edge
				// i and j are NOT in the same cluster
				processClusterImbMatrix(new_kj, new_ki) += fabs(weight);
			}
		}
		// NAO EH NECESSARIO ATUALIZAR O VERTEX IMBALANCE ABAIXO, POIS A BUSCA LOCAL (VND) EH FIRST IMPROVEMENT
		// vertexImbalance.push_back(std::make_pair(i, positiveSum + negativeSum));
	}
}

Coordinate ImbalanceSubgraphParallelILS::findMaximumElementInMatrix(matrix<double> &mat) {
	int x = 0, y = 0;
	double maxValue = -DBL_MAX;

	// TODO INCLUDE HANDLING FOR MATH PRECISION DOUBLE (EPSILON)
	for(int i = 0; i < mat.size1(); i++) {
		for(int j = 0; j < mat.size2(); j++) {
			if(mat(i, j) > maxValue) {
				maxValue = mat(i, j);
				x = i;
				y = j;
			}
		}
	}
	return Coordinate(x, y);
}

std::vector< Coordinate > ImbalanceSubgraphParallelILS::getMatrixElementsAsList(matrix<double> &mat) {
	std::vector< Coordinate > returnList;

	for(int i = 0; i < mat.size1(); i++) {
		for(int j = 0; j < mat.size2(); j++) {
			returnList.push_back(Coordinate(i, j, mat(i, j)));
		}
	}
	return returnList;
}

long ImbalanceSubgraphParallelILS::findMostImbalancedVertexInProcessPair(SignedGraph& g,
		ClusterArray& myCluster, Coordinate processPair) {

	long nc = myCluster.size();
	long n = g.getN();
	DirectedGraph::edge_descriptor e;
	boost::property_map<DirectedGraph, edge_properties_t>::type ew = boost::get(edge_properties, g.graph);
	double maxImbalance = -DBL_MAX;
	long maxVertex = 0;

	// For each vertex i
	for(long i = 0; i < n; i++) {
		long ki = myCluster[i];
		DirectedGraph::out_edge_iterator f, l;
		double positiveSum = double(0.0), negativeSum = double(0.0);
		// only processes vertexes belonging to the process pair
		if((ki == processPair.x) or (ki == processPair.y)) {
			// For each out edge of i
			for (boost::tie(f, l) = out_edges(i, g.graph); f != l; ++f) {
				e = *f;
				Vertex src = source(e, g.graph), targ = target(e, g.graph);
				double weight = ew[e].weight;
				long kj = myCluster[targ.id];
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

} /* namespace ils */
} /* namespace resolution */
