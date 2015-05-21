/*
 * ParallelNeighborhoodSearch.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#include "include/ParallelNeighborhoodSearch.h"
#include "include/SequentialNeighborhoodSearch.h"
#include "../util/include/MPIMessage.h"
#include "../util/parallel/include/MPIUtil.h"
#include "../problem/include/CCProblem.h"
#include "../problem/include/RCCProblem.h"

#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>

#include <limits>

namespace clusteringgraph {

using namespace boost::mpi;
using namespace util;
using namespace util::parallel;
using namespace std;

ParallelNeighborhoodSearch::ParallelNeighborhoodSearch(int allocStrategy,
		unsigned int offset, unsigned int numproc) :
		machineProcessAllocationStrategy(allocStrategy), numberOfSlaves(offset), numberOfSearchSlaves(numproc) {

}

ParallelNeighborhoodSearch::~ParallelNeighborhoodSearch() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
        Clustering* clustering, ClusteringProblem& problem,
        double timeSpentSoFar, double timeLimit, unsigned long randomSeed,
        int myRank, bool firstImprovementOnOneNeig) {

	// Resets the number of combinations tested on neighborhood search
	numberOfTestedCombinations = 0;
	// max number of clusters (RCC Problem Only)
	long k = 0;
	if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
		RCCProblem& rp = static_cast<RCCProblem&>(problem);
		k = rp.getK();
	}

	// Splits the processing in (n / numberOfSearchSlaves) chunks,
	// to be consumed by numberOfSearchSlaves processes
	unsigned long sizeOfChunk = g->getN() / (numberOfSearchSlaves + 1);
	unsigned long remainingVertices = g->getN() % (numberOfSearchSlaves + 1);
	mpi::communicator world;
	// the leader distributes the work across the processors
	// the leader itself (myRank) does part of the work too
	int i = 0;
	std::vector<int> slaveList;
	if(sizeOfChunk > 0) {
		MPIUtil::populateListOfVNDSlaves(machineProcessAllocationStrategy, slaveList, myRank, numberOfSlaves, numberOfSearchSlaves);
		// Sends the parallel search (VND) message to the search slaves via MPI
		for(i = 0; i < numberOfSearchSlaves; i++) {
			InputMessageParallelVND imsgpvns(g->getId(), l, g->getGraphFileLocation(), *clustering, problem.getType(),
					timeSpentSoFar, timeLimit, i * sizeOfChunk, (i + 1) * sizeOfChunk - 1, numberOfSlaves,
					numberOfSearchSlaves, k, false);
			world.send(slaveList[i], MPIMessage::INPUT_MSG_PARALLEL_VND_TAG, imsgpvns);
			BOOST_LOG_TRIVIAL(trace) << "VND Message sent to process " << slaveList[i] << "; [" << i * sizeOfChunk
					<< ", " << (i + 1) * sizeOfChunk - 1 << "]" << endl;
		}
	}
	// the leader (me) does its part of the work too
	BOOST_LOG_TRIVIAL(trace) << "Total number of vertices is " << g->getN() << endl;
	BOOST_LOG_TRIVIAL(trace) << "VND Parallelization status: The master plus " << numberOfSearchSlaves <<
			" search slaves will process " << sizeOfChunk << " vertices each one." << endl;
	BOOST_LOG_TRIVIAL(trace) << "RemainingVertices is " << remainingVertices << endl;

	double bestValue = numeric_limits<double>::infinity();
	Clustering bestClustering;
	// a busca abaixo, que eh feita pelo mestre, sera interrompida se ele receber uma solucao melhor de um dos escravos
	bestClustering = this->searchNeighborhood(l, g, clustering,
			problem, timeSpentSoFar, timeLimit, randomSeed, myRank,
			i * sizeOfChunk, g->getN() - 1, firstImprovementOnOneNeig, k);
	bestValue = bestClustering.getImbalance().getValue();
	if((firstImprovementOnOneNeig and l == 1) or (l == 2)) {
        	InputMessageParallelVND imsg;
                BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VND] First improvement on PVND: interrupting other VND slaves.";
                for(int j = 0; j < numberOfSearchSlaves; j++) {
                        world.send(slaveList[j], MPIMessage::INTERRUPT_MSG_PARALLEL_VND_TAG, imsg);
                }
        }
	// TODO enviar msg de interrupcao se a busca local do mestre acabar antes das dos demais escravos
	BOOST_LOG_TRIVIAL(debug) << "Waiting for slaves return messages...\n";
	// the leader receives the processing results, if that is the case
	// TODO implement global first improvement (l = 2-opt) in parallel VND here!
	if(sizeOfChunk > 0) {
		OutputMessage omsg;
		for(i = 0; i < numberOfSearchSlaves; i++) {
			mpi::status stat = world.recv(mpi::any_source, MPIMessage::OUTPUT_MSG_PARALLEL_VND_TAG, omsg);
			BOOST_LOG_TRIVIAL(debug) << "Message received from process " << stat.source() << ". Obj = " <<
					omsg.clustering.getImbalance().getValue();
			BOOST_LOG_TRIVIAL(trace) << omsg.clustering.toString(g->getN());
			// processes the result of the execution of process p(i)
			// sums the number of tested combinations
			numberOfTestedCombinations += omsg.numberOfTestedCombinations;
			// checks if the obj value improved
			if(omsg.clustering.getImbalance().getValue() < bestValue) {
				Clustering clustering = omsg.clustering;
				bestClustering = clustering;
				bestValue = omsg.clustering.getImbalance().getValue();
				BOOST_LOG_TRIVIAL(debug) << "*** [Parallel VND] Better value found for objective function in node "
						<< stat.source() << ": " <<
						omsg.clustering.getImbalance().getValue() << endl;
			}
		}
	}
	BOOST_LOG_TRIVIAL(debug) << "[Parallel VND] Best solution found: Obj = " << bestClustering.getImbalance().getValue() << endl;
	bestClustering.printClustering(g->getN());
	return bestClustering;
}

Clustering ParallelNeighborhoodSearch::searchNeighborhood(int l, SignedGraph* g,
		Clustering* clustering, ClusteringProblem& problem, double timeSpentSoFar,
		double timeLimit, unsigned long randomSeed, int myRank, unsigned long initialSearchIndex,
		unsigned long finalSearchIndex, bool firstImprovementOnOneNeig, unsigned long k) {

	assert(initialSearchIndex < g->getN());
	assert(finalSearchIndex < g->getN());

	if (l == 1) {  // 1-opt
		// Parallel search always does best improvement in 1-opt
		// Therefore, parameter firstImprovementOnOneNeig is ignored => TODO: DISABLED!!! Now the parameter is in use!
		return this->search1opt(g, clustering, problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, firstImprovementOnOneNeig, k);
	} else {  // 2-opt
		// TODO implement global first improvement in parallel 2-opt
		// first improvement found in one process must break all other processes' loop
		// IMPORTANT: Parallel VND does first improvement on 2-opt if CC Problem
		// Or best improvement on 2-opt if RCC Problem
		bool firstImprovementOn2Opt = true;
		if(problem.getType() == ClusteringProblem::CC_PROBLEM) {
			firstImprovementOn2Opt = true;
		} else if(problem.getType() == ClusteringProblem::RCC_PROBLEM) {
			// TODO check if this is still necessary with GRASP RCC
			// firstImprovementOn2Opt = false;
		}
		return this->search2opt(g, clustering, &problem, timeSpentSoFar, timeLimit, randomSeed,
				myRank, initialSearchIndex, finalSearchIndex, firstImprovementOn2Opt, k);
	}
}


} /* namespace graph */
