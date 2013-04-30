/*
 * Grasp.h
 *
 *  Created on: 30/04/2013
 *      Author: czt0
 */

#ifndef GRASP_H_
#define GRASP_H_

#include "../../include/ResolutionStrategy.h"

namespace resolution {
namespace grasp {

class Grasp: public ResolutionStrategy {
public:
	Grasp();
	virtual ~Grasp();
};

} /* namespace grasp */
} /* namespace resolution */

#endif /* GRASP_H_ */
