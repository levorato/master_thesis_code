/*
 * MovieLensSGConverter.h
 *
 *  Created on: Feb 1, 2016
 *      Author: mlevorato
 */

#ifndef SRC_MOVIELENSSGCONVERTER_H_
#define SRC_MOVIELENSSGCONVERTER_H_

#include <map>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::numeric::ublas;
using namespace std;
using boost::multiprecision::cpp_dec_float_50;

namespace generation {

class MovieLensSGConverter {
public:
	MovieLensSGConverter();
	~MovieLensSGConverter();
	
	/**
	 * pos_edge_perc -> edge percentual when comparing 2 users and assuming their relation is positive (e.g. 80%)
	 * neg_edge_perc -> edge percentual when comparing 2 users and assuming their relation is negative (e.g. 20%)
	 * number_chunks -> the number of chunks the matrix will be split for processing (saves memory)
	 */
	bool processMovieLensFolder(const string& folder, const string& filter,
			const unsigned int &myRank, const unsigned int &numProcessors, const double& minimum_edge_weight,
			const int& number_chunks);
	bool readMovieLensCSVFile(const string& filename, long& max_user_id, long& max_movie_id);
	bool generateSGFromMovieRatings(const long& max_user_id, const long& max_movie_id,
			const string& outputFileName, const unsigned int &myRank,
			const unsigned int &numProcessors, const double& minimum_edge_weight,
			int number_chunks);
	
private:
	// the movie_users structure maps a movie_id (long) to a vector of <user_id, rating> pairs
	std::vector< std::vector< std::pair<long, char> > > movie_users;
	// the sparse matrix with movie ratings
	compressed_matrix<double> star;

	std::string get_file_contents(const char *filename);
	void find_and_replace(string& source, string const& find, string const& replace);
	cpp_dec_float_50 cosine_similarity(std::vector<double>& votesFromUserA,
			std::vector<double>& votesFromUserB);
	cpp_dec_float_50 pearson_correlation_coefficient(std::vector<double>& votesFromUserA,
			std::vector<double>& votesFromUserB);
	// http://projekter.aau.dk/projekter/files/32181941/Report.pdf - pages 24/25
	cpp_dec_float_50 pearson_correlation_coefficient2(std::vector<double>& votesFromUserA,
				std::vector<double>& votesFromUserB);
	cpp_dec_float_50 spearman_correlation_coefficient(std::vector<int>& votesFromUserA,
			std::vector<int>& votesFromUserB);
	bool is_zero_array(std::vector<double>& array);

	bool writeSignedGraphToFile(const string& file_prefix, const string& partialFilename,
			const std::vector<string>& edgeList, const long& max_user_id, const int& myRank);
	bool writeHistogramToFile(const string& file_prefix, const string& partialFilename,
			string fileContent, const int& myRank);

	bool mergeSignedGraphToFile(const string& file_prefix, const string& partialFilename,
			const std::vector<string>& edgeList, const long& max_user_id, const int& numProcessors);
	bool mergeHistogramToFile(const string& file_prefix, const string& partialFilename,
			const int& numProcessors, std::map<string, long>& histogram);
};

class InputMessage {
public:
	// the maximum common rating count
	long max;
	static const int MAX_COUNT_MSG_TAG = 10;
	static const int TERMINATE_MSG_TAG = 20;
	static const int DONE_MSG_TAG = 30;

	InputMessage() : max(0) {

	}

	InputMessage(long i) :
		max(i) {
	}

	~InputMessage(){};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & max;
	}
};

} /* namespace generation */

#endif /* SRC_MOVIELENSSGCONVERTER_H_ */
