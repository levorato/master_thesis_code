//
// Created by mlevorato on 15/03/17.
//

#ifndef MESTRADO_PRECISION_H
#define MESTRADO_PRECISION_H

#include <cmath>
#include <limits>

class Precision {

public:
    static bool nearlyEqual(double a, double b, double epsilon = std::numeric_limits<double>::epsilon()) {
        double absA = fabs(a);
        double absB = fabs(b);
        double diff = fabs(a - b);

        if (a == b) { // shortcut, handles infinities
            return true;
        } else if (a == 0 || b == 0 || diff < std::numeric_limits<double>::min()) {
            // a or b is zero or both are extremely close to it
            // relative error is less meaningful here
            return diff < (epsilon * std::numeric_limits<double>::min());
        } else { // use relative error
            return diff / (absA + absB) < epsilon;
        }
    }

    static bool rough_lt(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon()) // operator<
    {
        return rhs - lhs >= epsilon;
        // tricky >= because if the difference is equal to epsilon
        // then they are not equal per the rough_eq method
    }

    static bool rough_lte(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon()) // operator<=
    {
        return rhs - lhs > -epsilon;
    }

    static bool rough_gt(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon()) // operator>
    {
        return lhs - rhs >= epsilon;
        // tricky >= because if the difference is equal to epsilon
        // then they are not equal per the rough_eq method
    }

    static bool rough_gte(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon()) // operator>=
    {
        return lhs - rhs > -epsilon;
    }

};


#endif //MESTRADO_PRECISION_H
