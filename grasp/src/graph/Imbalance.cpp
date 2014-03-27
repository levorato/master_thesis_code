/*
 * Imbalance.cpp
 *
 *  Created on: 13/06/2013
 *      Author: czt0
 */

#include "include/Imbalance.h"
#include <iomanip>

namespace clusteringgraph {

Imbalance::Imbalance(double positive, double negative) : positiveValue(positive),
		negativeValue(negative) {
	// TODO Auto-generated constructor stub

}

Imbalance::Imbalance(const Imbalance& i) : positiveValue(i.positiveValue),
		negativeValue(i.negativeValue) {

}

Imbalance::~Imbalance() {
	// TODO Auto-generated destructor stub
}

Imbalance& Imbalance::operator+=(const Imbalance &i) {
	this->positiveValue += (i.positiveValue);
	this->negativeValue += (i.negativeValue);
	return *this;
}

Imbalance& Imbalance::operator-=(const Imbalance &i) {
	this->positiveValue -= (i.positiveValue);
	this->negativeValue -= (i.negativeValue);
	return *this;
}

bool operator> (Imbalance &i1, Imbalance &i2) {
	return i1.getValue() > i2.getValue();
}

bool operator<= (Imbalance &i1, Imbalance &i2) {
	return i1.getValue() <= i2.getValue();
}

bool operator< (Imbalance &i1, Imbalance &i2) {
	return i1.getValue() < i2.getValue();
}

bool operator>= (Imbalance &i1, Imbalance &i2) {
	return i1.getValue() >= i2.getValue();
}

ostream& Imbalance::operator<<(ostream &out) {     //output
	out << "I(P) = " << std::fixed << std::setprecision(2) << getValue()
			<< " (" << positiveValue << "+, " << negativeValue << "-)";
	return out;
}

} /* namespace clusteringgraph */
