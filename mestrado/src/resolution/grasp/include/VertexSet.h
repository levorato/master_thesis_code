/*
 * VertexSet.h
 *
 *  Created on: May 7, 2013
 *      Author: mario
 */

#ifndef VERTEXSET_H_
#define VERTEXSET_H_

#include <set>
#include <boost/shared_ptr.hpp>

namespace resolution {
namespace grasp {

using namespace std;
using namespace boost;

// TODO implement the gain function according to the CC problem.
class GainFunctionComparison
{
  bool reverse;
public:
  GainFunctionComparison(const bool& revparam=false)
    {reverse=revparam;}
  bool operator() (const int& lhs, const int&rhs) const
  {
    if (reverse) return (lhs>rhs);
    else return (lhs<rhs);
  }
};

typedef set<int, GainFunctionComparison> GainFunctionVertexSet;
typedef boost::shared_ptr<GainFunctionVertexSet> GainFunctionVertexSetPtr;

class VertexSet {
public:
	/**
	 * Creates a vertex set containing nodes 0..n-1.
	 */
	VertexSet(int n);
	virtual ~VertexSet();

	/**
	 * Returns the size of the vertex list.
	 */
	int size();
	/**
	 * Removes the specified vertex from the list.
	 */
	void removeVertex(int i);
private:
	GainFunctionVertexSetPtr vertexSetPtr;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* VERTEXSET_H_ */
