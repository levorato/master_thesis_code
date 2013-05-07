/*
 * VertexSet.cpp
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#include "include/VertexSet.h"

namespace resolution {
namespace grasp {

VertexSet::VertexSet(int n) : vertexSetPtr(new GainFunctionVertexSet) {
	for(int i = 0; i < n; i++) {
		vertexSetPtr->insert(i);
	}
}

VertexSet::~VertexSet() {
	// TODO Auto-generated destructor stub
}

int VertexSet::size() {
	return vertexSetPtr->size();
}

void VertexSet::removeVertex(int i) {
	vertexSetPtr->erase(i);
}

} /* namespace grasp */
} /* namespace resolution */
