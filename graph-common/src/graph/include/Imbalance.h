/*
 * Imbalance.h
 *
 *  Created on: 13/06/2013
 *      Author: Mario Levorato
 */

#ifndef IMBALANCE_H_
#define IMBALANCE_H_

#include <iostream>
#include <boost/serialization/access.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::numeric::ublas;

namespace clusteringgraph {

/**
 * Models the imbalance value of a given clustering associated to a graph.
 */
class Imbalance {
public:
	Imbalance();
	Imbalance(double positive, double negative);
	Imbalance(const Imbalance& i);
	virtual ~Imbalance();

	Imbalance& operator+=(const Imbalance &i);
	Imbalance& operator-=(const Imbalance &i);

    friend bool operator> (Imbalance &i1, Imbalance &i2);
    friend bool operator<= (Imbalance &i1, Imbalance &i2);
    friend bool operator< (Imbalance &i1, Imbalance &i2);
    friend bool operator>= (Imbalance &i1, Imbalance &i2);
    friend bool operator== (Imbalance &i1, Imbalance &i2);

    ostream &operator<<(ostream &out);

	double getValue() const {
		return positiveValue + negativeValue;
	}

	double getNegativeValue() const {
		return negativeValue;
	}

	void setNegativeValue(double negativeValue) {
		this->negativeValue = negativeValue;
	}

	double getPositiveValue() const {
		return positiveValue;
	}

	void setPositiveValue(double positiveValue) {
		this->positiveValue = positiveValue;
	}

private:
	double positiveValue;
	double negativeValue;

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & positiveValue;
		ar & negativeValue;
	}
};

class ImbalanceMatrix {
public:
	matrix<double> pos, neg;
	ImbalanceMatrix() : pos(), neg() { }
	ImbalanceMatrix(int nc) : pos(zero_matrix<double>(nc, nc)), neg(zero_matrix<double>(nc, nc)) { }

	ImbalanceMatrix& operator+=(const ImbalanceMatrix &i) {
		this->pos += (i.pos);
		this->neg += (i.neg);
		return *this;
	}

	// serialization-specific code
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & pos;
		ar & neg;
	}
};

} /* namespace clusteringgraph */
#endif /* IMBALANCE_H_ */
