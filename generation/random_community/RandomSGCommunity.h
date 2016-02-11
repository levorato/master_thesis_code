/*
 * RandomSGCommunity.h
 *
 *  Created on: Feb 1, 2016
 *      Author: mlevorato
 */

#ifndef SRC_RandomSGCommunity_H_
#define SRC_RandomSGCommunity_H_

#include <map>
#include <list>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

using namespace boost::numeric::ublas;
using namespace std;

// generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > >
typedef generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > > Matrix;

namespace generation {

class Edge {
public:
	long x, y;

	Edge() : x(-1), y(-1) {  }
	Edge(long a, long b) : x(a), y(b) {  }
};

class RandomSGCommunity {
public:
	RandomSGCommunity();
	~RandomSGCommunity();
	
	bool SG(const long& c, const long& n, const long& k,
				const double& p_in, const double& p_minus, const double& p_plus);
	bool generateRandomSG(const long& c, const long& n, const long& k,
			const double& p_in, const double& p_minus, const double& p_plus,
			const unsigned int &myRank, const unsigned int &numProcessors);
	
private:
	std::vector<long> mycluster;
	std::vector< std::vector<long> >cluster_node_list;
	std::list<long> vertex_list;

	std::string get_file_contents(const char *filename);
	void find_and_replace(string& source, string const& find, string const& replace);
	
	long pick_random_vertex(std::vector<long>& indegree, std::vector<long>& outdegree, const long& k);
	Edge pick_random_external_edge(const long& N, Matrix& matrix, std::vector<long>& indegree,
			std::vector<long>& outdegree, const long& k, const long& c);
	Edge pick_random_internal_edge(const long& N, Matrix& matrix, std::vector<long>& indegree,
			std::vector<long>& outdegree, const long& k);
	double CCObjectiveFunction(const long& N, Matrix& matrix);
	
};

class InputMessage {
public:
	// the identifier of the graph
	unsigned int id;
	static const int TERMINATE_MSG_TAG;
	static const int DONE_MSG_TAG;

	InputMessage() : id(0) {

	}

	InputMessage(unsigned int i) :
		id(i) {
	}

	~InputMessage(){};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & id;
	}
};

} /* namespace generation */

#endif /* SRC_RandomSGCommunity_H_ */
