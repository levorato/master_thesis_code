/*
 * CCProblem.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef CCPROBLEM_H_
#define CCPROBLEM_H_

#include "ClusteringProblem.h"

namespace problem {

class CCProblem: public problem::ClusteringProblem {
public:
	CCProblem();
	virtual ~CCProblem();

	int gainFunction();
};

} /* namespace problem */
#endif /* CCPROBLEM_H_ */
