/*
 * NeighborhoodSearchFactory.cpp
 *
 *  Created on: 09/07/2013
 *      Author: mlevorato
 */

#include "include/NeighborhoodSearchFactory.h"

namespace resolution {
namespace grasp {

NeighborhoodSearchFactory::NeighborhoodSearchFactory(unsigned long numberOfSlaves,
		unsigned long numberOfSearchSlaves) : sequentialNeighborhoodSearch(),
				parallelNeighborhoodSearch(numberOfSlaves, numberOfSearchSlaves) {

}

NeighborhoodSearchFactory::~NeighborhoodSearchFactory() {
	// TODO Auto-generated destructor stub
}

} /* namespace grasp */
} /* namespace resolution */
