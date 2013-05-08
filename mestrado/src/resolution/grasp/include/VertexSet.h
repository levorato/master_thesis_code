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

// TODO implement the gain function according to the gain function
// gc(i) specified in the article.
// See Class ClusteringProblem.
class GainFunctionComparison : std::binary_function <int, int, bool>
{
  bool reverse;
public:
  GainFunctionComparison(const bool& revparam=false)
    {reverse=revparam;}
    bool operator () ( const int& a, const int& b ) const
    {
      return a < b;
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
	/**
	 * Picks a random vertex among the first x elements.
	 */
	int chooseRandomVertex(int x);
private:
	GainFunctionVertexSetPtr vertexSetPtr;
};

} /* namespace grasp */
} /* namespace resolution */
#endif /* VERTEXSET_H_ */
