/*
 * RandomUtil.h
 *
 *  Created on: 06/05/2014
 *      Author: czt0
 */

#ifndef RANDOMUTIL_H_
#define RANDOMUTIL_H_

#include <boost/random/uniform_int.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>

namespace util {

// other options of generator: minstd_rand, mt19937
typedef boost::minstd_rand RandomGeneratorType;

class RandomUtil {
public:
	virtual ~RandomUtil();
	RandomUtil();

	static void setSeed(int seed) {
		// other option for seed: generator.seed(boost::random::random_device()());
		rg.seed(seed);
		isSeeded = true;
	}

	static int next(int lowerLimit, int upperLimit) {
		// we are supposing the random generator is already seeded here
		rg.seed(boost::random::random_device()());
		boost::uniform_int<> distribution(lowerLimit, upperLimit);
		boost::variate_generator<RandomGeneratorType&, boost::uniform_int<> > LimitedInt(
				rg, distribution);
		return LimitedInt();
	}

private:
	static RandomGeneratorType rg;
	static bool isSeeded;
};

} /* namespace util */
#endif /* RANDOMUTIL_H_ */
