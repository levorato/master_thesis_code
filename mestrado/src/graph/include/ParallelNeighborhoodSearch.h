/*
 * ParallelNeighborhoodSearch.h
 *
 *  Created on: Jul 2, 2013
 *      Author: mario
 */

#ifndef PARALLELNEIGHBORHOODSEARCH_H_
#define PARALLELNEIGHBORHOODSEARCH_H_

#include "Neighborhood.h"

namespace resolution {
namespace grasp {

class ParallelNeighborhoodSearch: public clusteringgraph::NeighborhoodSearch {
public:
	ParallelNeighborhoodSearch();
	virtual ~ParallelNeighborhoodSearch();
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* PARALLELNEIGHBORHOODSEARCH_H_ */
