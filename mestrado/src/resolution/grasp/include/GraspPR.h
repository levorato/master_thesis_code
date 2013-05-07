/*
 * GraspPR.h
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#ifndef GRASPPR_H_
#define GRASPPR_H_

#include "../../include/ResolutionStrategy.h"
#include "../../../graph/include/Graph.h"
#include "../../../graph/include/Clustering.h"

using namespace clusteringgraph;

namespace resolution {
namespace grasp {

class GraspPR: public ResolutionStrategy {
public:
	GraspPR();
	virtual ~GraspPR();

	void pathRelinking(SignedGraph* g);
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* GRASPPR_H_ */
