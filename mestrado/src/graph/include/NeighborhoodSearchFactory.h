/*
 * NeighborhoodSearchFactory.h
 *
 *  Created on: 09/07/2013
 *      Author: mlevorato
 */

#ifndef NEIGHBORHOODSEARCHFACTORY_H_
#define NEIGHBORHOODSEARCHFACTORY_H_

#include "SequentialNeighborhoodSearch.h"
#include "ParallelNeighborhoodSearch.h"
#include "Neighborhood.h"

namespace resolution {
namespace grasp {

class NeighborhoodSearchFactory {
public:
	static const int SEQUENTIAL = 0, PARALLEL = 1;

	NeighborhoodSearchFactory(unsigned long numberOfSlaves, unsigned long numberOfSearchSlaves);
	virtual ~NeighborhoodSearchFactory();

	SequentialNeighborhoodSearch sequentialNeighborhoodSearch;
	ParallelNeighborhoodSearch parallelNeighborhoodSearch;

	NeighborhoodSearch* build(int neighborhoodType) {
		if(neighborhoodType == NeighborhoodSearchFactory::SEQUENTIAL) {
			return &sequentialNeighborhoodSearch;
		} else {
			return &parallelNeighborhoodSearch;
		}
	}
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* NEIGHBORHOODSEARCHFACTORY_H_ */
