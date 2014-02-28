/*
 * MPIUtil.h
 *
 *  Created on: 25/02/2014
 *      Author: czt0
 */

#ifndef MPIUTIL_H_
#define MPIUTIL_H_

#include <vector>
#include <boost/log/trivial.hpp>

namespace util {
namespace parallel {

using namespace boost;
using namespace std;

class MPIUtil {
public:
	MPIUtil();
	virtual ~MPIUtil();

	/**
	 * Returns true if a process rank is considered a GRASP Slave process.
	 */
	static bool isGRASPSlave(const unsigned int& myRank, const unsigned int& numberOfSlaves,
			const unsigned int& numberOfSearchSlaves) {
		if(MACHINE_PROCESS_ALLOCATION_STRATEGY == ALL_GRASP_SLAVES_FIRST) {
			return myRank <= numberOfSlaves;
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER)
			return myRank % (numberOfSearchSlaves + 1) == 0;
		}
	}

	/**
	 * Returns true if a process rank is considered a VNS Slave process.
	 */
	static bool isVNSSlave(const unsigned int& myRank, const unsigned int& numberOfSlaves,
			const unsigned int& numberOfSearchSlaves) {
		return not isGRASPSlave(myRank, numberOfSlaves, numberOfSearchSlaves);
	}

	static unsigned int calculateNumberOfSearchSlaves(const unsigned int& np, const unsigned int& numberOfSlaves) {
		// Number of remaining slaves after using 1 (leader) + 'numberOfSlaves' processes for parallel GRASP execution
		unsigned int remainingSlaves = np - numberOfSlaves - 1;
		// Divides the remaining slaves that will be used in parallel VNS processing
		unsigned int numberOfSearchSlaves = remainingSlaves / (numberOfSlaves + 1);
		BOOST_LOG_TRIVIAL(info) << "The number of VNS search slaves per master is " << numberOfSearchSlaves;

		return numberOfSearchSlaves;
	}

	static void populateListOfGRASPSlaves(std::vector<int>& slaveList, const unsigned int& myRank,
			const unsigned int& numberOfSlaves, const unsigned int& numberOfSearchSlaves) {
		if(MACHINE_PROCESS_ALLOCATION_STRATEGY == ALL_GRASP_SLAVES_FIRST) {
			for(int i = 1; i <= numberOfSlaves; i++) {
				slaveList.push_back(i);
			}
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER)
			for(int i = 1; i <= numberOfSlaves; i++) {
				slaveList.push_back(i * (numberOfSearchSlaves + 1));
			}
		}
	}

	static void populateListOfVNSSlaves(std::vector<int>& slaveList, const unsigned int& myRank,
			const unsigned int& numberOfSlaves, const unsigned int& numberOfSearchSlaves) {
		if(MACHINE_PROCESS_ALLOCATION_STRATEGY == ALL_GRASP_SLAVES_FIRST) {
			// The formulas below determine the first and last search slaves of this grasp slave process
			// IMPORTANT! Process mapping:
			// * 1..numberOfSlaves => p(i) GRASP slave processes (execute parellel GRASP iterations with MPI)
			// * numberOfSlaves+1..numberOfSearchSlaves => VNS slave processes for p(0) (execute parallel VNS for p(0))
			// * numberOfSlaves+1+i*numberOfSearchSlaves..numberOfSlaves+(i+1)*numberOfSearchSlaves => VNS slave processes for p(i)
			// Here, p(i) is represented by myRank.
			int firstSlave = numberOfSlaves + 1 + myRank * numberOfSearchSlaves;
			int lastSlave = numberOfSlaves + (myRank + 1) * numberOfSearchSlaves;
			for(int i = firstSlave; i <= lastSlave; i++) {
				slaveList.push_back(i);
			}
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER)
			// The formulas below determine the first and last search slaves of this grasp slave process
			// IMPORTANT! Process mapping:
			// * myRank * (numberOfSearchSlaves+1) == 0 => p(i) GRASP slave processes (execute parellel GRASP iterations with MPI)
			// * 1..numberOfSearchSlaves => VNS slave processes for p(0) (execute parallel VNS for p(0))
			// * [myRank+1]..[myRank + numberOfSearchSlaves] => VNS slave processes for p(i), where i = myRank
			int firstSlave = myRank + 1;
			int lastSlave = myRank + numberOfSearchSlaves;
			for(int i = firstSlave; i <= lastSlave; i++) {
				slaveList.push_back(i);
			}
		}
	}

	// Machine x process allocation stategies
	static const int ALL_GRASP_SLAVES_FIRST = 0;
	static const int GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER = 1;

	// The following variable defines the strategy for machine x process allocation
	static const int MACHINE_PROCESS_ALLOCATION_STRATEGY = ALL_GRASP_SLAVES_FIRST;
};

} /* namespace parallel */
} /* namespace util */
#endif /* MPIUTIL_H_ */
