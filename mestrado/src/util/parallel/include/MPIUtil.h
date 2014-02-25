/*
 * MPIUtil.h
 *
 *  Created on: 25/02/2014
 *      Author: czt0
 */

#ifndef MPIUTIL_H_
#define MPIUTIL_H_

#include <boost/log/trivial.hpp>

namespace util {
namespace parallel {

using namespace boost;

class MPIUtil {
public:
	MPIUtil();
	virtual ~MPIUtil();

	/**
	 * Returns true if a process rank is considered a GRASP Slave process.
	 */
	static bool isGRASPSlave(const& unsigned int myRank, const& unsigned int numberOfSearchSlaves) {
		if(MACHINE_PROCESS_ALLOCATION_STRATEGY == ALL_GRASP_SLAVES_FIRST) {
			return myRank <= numberOfSlaves;
		} else { // if(MACHINE_PROCESS_ALLOCATION_STRATEGY == GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER)
			return myRank % (numberOfSearchSlaves + 1) == 0;
		}
	}

	/**
	 * Returns true if a process rank is considered a VNS Slave process.
	 */
	static bool isVNSSlave(const& unsigned int myRank, const& unsigned int numberOfSearchSlaves) {
		return not isGRASPSlave(myRank, numberOfSearchSlaves);
	}

	static unsigned int calculateNumberOfSearchSlaves(const unsigned int& np, const unsigned int& numberOfSlaves) {
		// Number of remaining slaves after using 1 (leader) + 'numberOfSlaves' processes for parallel GRASP execution
		unsigned int remainingSlaves = np - numberOfSlaves - 1;
		// Divides the remaining slaves that will be used in parallel VNS processing
		unsigned int numberOfSearchSlaves = remainingSlaves / (numberOfSlaves + 1);
		BOOST_LOG_TRIVIAL(info) << "The number of VNS search slaves per master is " << numberOfSearchSlaves;

		return numberOfSearchSlaves;
	}

	void populateListOfGRASPSlaves(std::vector<int>& slaveList, const& unsigned int myRank,
				const& unsigned int numberOfSearchSlaves) {
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



	// Machine x process allocation stategies
	static const int ALL_GRASP_SLAVES_FIRST = 0;
	static const int GRASP_SLAVE_AND_VNS_SLAVES_TOGETHER = 0;

	// The following variable defines the strategy for machine x process allocation
	static const int MACHINE_PROCESS_ALLOCATION_STRATEGY = ALL_GRASP_SLAVES_FIRST;
};

} /* namespace parallel */
} /* namespace util */
#endif /* MPIUTIL_H_ */
