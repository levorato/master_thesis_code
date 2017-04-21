/*
 * CCImbalanceCalculator.h
 *
 *  Created on: 21 de abr de 2017
 *      Author: mlevorato
 */

#ifndef SRC_RESOLUTION_VALIDATION_INCLUDE_CCIMBALANCECALCULATOR_H_
#define SRC_RESOLUTION_VALIDATION_INCLUDE_CCIMBALANCECALCULATOR_H_

#include "./Graph.h"
#include "./Clustering.h"
#include "./CCProblem.h"
#include "./SimpleTextGraphFileReader.h"

namespace clusteringgraph {
namespace validation {

class CCImbalanceCalculator {
private:
	CCImbalanceCalculator() : g(), built(false) {

	}

public:
	void build(const std::string &filepath) {
		SimpleTextGraphFileReader reader;
		g = reader.readGraphFromFile(filepath);
	}

	static CCImbalanceCalculator& instance(const std::string &filepath)
	{
	  static CCImbalanceCalculator INSTANCE;
	  if(not INSTANCE.built) {
		  INSTANCE.build(filepath);
	  }
	  return INSTANCE;
	}

	Imbalance objectiveFunction(ClusterArray& cArray) {
		clusteringgraph::validation::CCProblem problem;
		clusteringgraph::validation::Clustering validation(cArray, *g, problem);
		return validation.getImbalance();
	}

	SignedGraphPtr g;
	bool built;
};

} /* namespace validation */
} /* namespace clusteringgraph */

#endif /* SRC_RESOLUTION_VALIDATION_INCLUDE_CCIMBALANCECALCULATOR_H_ */
