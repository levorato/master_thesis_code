/*
 * MPIUtil.h
 *
 *  Created on: 25/02/2014
 *      Author: Mario Levorato
 */

#ifndef MPIUTIL_H_
#define MPIUTIL_H_

#include <vector>
#include <boost/log/trivial.hpp>
#include <boost/serialization/access.hpp>

namespace util {
namespace parallel {

using namespace boost;
using namespace std;

class MPIUtil {
public:
	MPIUtil();
	virtual ~MPIUtil();

	/**
	 * Returns true if a process rank is considered a master process.
	 */
	static bool isMaster(int machineProcessAllocationStrategy, const unsigned int& myRank, const unsigned int& numberOfMasters,
			const unsigned int& numberOfSearchSlaves) {
		if(machineProcessAllocationStrategy == ALL_MASTERS_FIRST) {
			return myRank <= numberOfMasters;
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == MASTER_AND_VND_SLAVES_TOGETHER)
			return myRank % (numberOfSearchSlaves + 1) == 0;
		}
	}

	/**
	 * Returns true if a process rank is considered a VND Slave process.
	 */
	static bool isVNDSlave(int machineProcessAllocationStrategy, const unsigned int& myRank, const unsigned int& numberOfMasters,
			const unsigned int& numberOfSearchSlaves) {
		return not isMaster(machineProcessAllocationStrategy, myRank, numberOfMasters, numberOfSearchSlaves);
	}

	static unsigned int calculateNumberOfSearchSlaves(const unsigned int& np, const unsigned int& numberOfMasters) {
		// Number of remaining slaves after using 1 (leader) + 'numberOfMasters' processes for parallel GRASP execution
		unsigned int remainingSlaves = np - numberOfMasters - 1;
		// Divides the remaining slaves that will be used in parallel VND processing
		unsigned int numberOfSearchSlaves = remainingSlaves / (numberOfMasters + 1);
		BOOST_LOG_TRIVIAL(info) << "The number of VND search slaves per master is " << numberOfSearchSlaves;

		return numberOfSearchSlaves;
	}

	static void populateListOfMasters(int machineProcessAllocationStrategy,
			std::vector<int>& masterList, const unsigned int& myRank,
			const unsigned int& numberOfMasters, const unsigned int& numberOfSearchSlaves) {
		if(machineProcessAllocationStrategy == ALL_MASTERS_FIRST) {
			for(int i = 1; i <= numberOfMasters; i++) {
				masterList.push_back(i);
			}
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == MASTER_AND_VND_SLAVES_TOGETHER)
			for(int i = 1; i <= numberOfMasters; i++) {
				masterList.push_back(i * (numberOfSearchSlaves + 1));
			}
		}
	}

	static void populateListOfVNDSlaves(int machineProcessAllocationStrategy, std::vector<int>& slaveList, const unsigned int& myRank,
			const unsigned int& numberOfMasters, const unsigned int& numberOfSearchSlaves) {
		if(machineProcessAllocationStrategy == ALL_MASTERS_FIRST) {
			// The formulas below determine the first and last search slaves of this master process
			// IMPORTANT! Process mapping:
			// * 1..numberOfSlaves => p(i) master processes (execute parellel GRASP iterations with MPI)
			// * numberOfSlaves+1..numberOfSearchSlaves => VND slave processes for p(0) (execute parallel VND for p(0))
			// * numberOfSlaves+1+i*numberOfSearchSlaves..numberOfSlaves+(i+1)*numberOfSearchSlaves => VND slave processes for p(i)
			// Here, p(i) is represented by myRank.
			int firstSlave = numberOfMasters + 1 + myRank * numberOfSearchSlaves;
			int lastSlave = numberOfMasters + (myRank + 1) * numberOfSearchSlaves;
			for(int i = firstSlave; i <= lastSlave; i++) {
				slaveList.push_back(i);
			}
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == MASTER_AND_VND_SLAVES_TOGETHER)
			// The formulas below determine the first and last search slaves of this master process
			// IMPORTANT! Process mapping:
			// * myRank * (numberOfSearchSlaves+1) == 0 => p(i) master processes (execute parellel GRASP iterations with MPI)
			// * 1..numberOfSearchSlaves => VND slave processes for p(0) (execute parallel VND for p(0))
			// * [myRank+1]..[myRank + numberOfSearchSlaves] => VND slave processes for p(i), where i = myRank
			int firstSlave = myRank + 1;
			int lastSlave = myRank + numberOfSearchSlaves;
			for(int i = firstSlave; i <= lastSlave; i++) {
				slaveList.push_back(i);
			}
		}
	}

	// Machine x process allocation stategies
	static const int ALL_MASTERS_FIRST = 0;
	static const int MASTER_AND_VND_SLAVES_TOGETHER = 1;

	enum SearchName {SEQUENTIAL_SEARCH, PARALLEL_SEARCH};
};


class MPIInitParams {
public:
	MPIInitParams() : numberOfMasters(0), numberOfSearchSlavesPerMaster(0),
			machineProcessAllocationStrategy(MPIUtil::MASTER_AND_VND_SLAVES_TOGETHER),
			searchType(0), isParallelGraph(false), N(0)  { }

	MPIInitParams(int masters, int slavesPerMaster, int machineAllocationStrategy, int searchTp, bool pgraph, long numVertices = 0) :
			numberOfMasters(masters),
			numberOfSearchSlavesPerMaster(slavesPerMaster),
			machineProcessAllocationStrategy(machineAllocationStrategy),
			searchType(searchTp), isParallelGraph(pgraph), N(numVertices)  { }

	int numberOfMasters;
	int numberOfSearchSlavesPerMaster;
	// The following variable defines the strategy for machine x process allocation
	int machineProcessAllocationStrategy;
	int searchType;
	bool isParallelGraph;
	long N;  // the number of vertices of the global graph
private:

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & numberOfMasters;
		ar & numberOfSearchSlavesPerMaster;
		ar & machineProcessAllocationStrategy;
		ar & searchType;
		ar & isParallelGraph;
		ar & N;
	}
};

} /* namespace parallel */
} /* namespace util */
#endif /* MPIUTIL_H_ */
