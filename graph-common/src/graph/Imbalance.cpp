/*
 * Imbalance.cpp
 *
 *  Created on: 13/06/2013
 *      Author: Mario Levorato
 */

#include "include/Imbalance.h"
#include <iomanip>
#include "../util/include/Precision.h"

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

bool Imbalance::operator==(const Imbalance &other) {
    double lhs = this->getValue();
    double rhs = other.getValue();
    return Precision::nearlyEqual(lhs, rhs);
}

bool operator> (Imbalance &i1, Imbalance &i2) {
	return Precision::rough_gt(i1.getValue(), i2.getValue());
}

bool operator<= (Imbalance &i1, Imbalance &i2) {
	return Precision::rough_lte(i1.getValue(), i2.getValue());
}

bool operator< (Imbalance &i1, Imbalance &i2) {
	return Precision::rough_lt(i1.getValue(), i2.getValue());
}

bool operator>= (Imbalance &i1, Imbalance &i2) {
	return Precision::rough_gte(i1.getValue(), i2.getValue());
}

ostream& Imbalance::operator<<(ostream &out) {     //output
	out << "I(P) = " << std::fixed << std::setprecision(2) << getValue()
			<< " (" << positiveValue << "+, " << negativeValue << "-)";
	return out;
}

} /* namespace clusteringgraph */
