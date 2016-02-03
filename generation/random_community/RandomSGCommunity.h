/*
 * RandomSGCommunity.h
 *
 *  Created on: Feb 1, 2016
 *      Author: mlevorato
 */

#ifndef SRC_RandomSGCommunity_H_
#define SRC_RandomSGCommunity_H_

#include <map>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <list>
#include <vector>

using namespace boost::numeric::ublas;
using namespace std;

namespace generation {

class RandomSGCommunity {
public:
	RandomSGCommunity();
	~RandomSGCommunity();
	
	bool generateRandomSG(const long& c, const long& n, const long& k,
			const double& p_in, const double& p_minus, const double& p_plus,
			const unsigned int &myRank, const unsigned int &numProcessors);
	bool generateSGFromMovieRatings(const long& max_user_id, const long& max_movie_id,
			const string& outputFileName, const unsigned int &myRank,
			const unsigned int &numProcessors);
	
private:
	std::vector<long> mycluster;
	std::vector< std::vector<long> >cluster_node_list;
	std::list<long> vertex_list;
	// original directed graph from file
	std::vector< std::vector< std::pair<long, char> > > matrix;

	std::string get_file_contents(const char *filename);
	void find_and_replace(string& source, string const& find, string const& replace);
	
	long pick_random_vertex(vector<long>& indegree, vector<long>& outdegree, const long& k);
	Edge pick_random_external_edge(const long& N, const long& indegree, const long& outdegree, 
		const long& k, const long& c);
	Edge pick_random_internal_edge(const long& N, const long& indegree, const long& outdegree, const long& k);
	double CCObjectiveFunction(const long& N);
	
};

class Edge {
	long x, y;
	
	Edge() {  }
	Edge(long a, long b) : x(a), y(b) {  }
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
