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

using namespace boost::numeric::ublas;
using namespace std;

namespace generation {

class MovieLensSGConverter {
public:
	MovieLensSGConverter();
	~MovieLensSGConverter();
	
	bool processMovieLensFolder(const string& folder, const string& filter);
	bool readMovieLensCSVFile(const string& filename, long& max_user_id, long& max_movie_id);
	bool generateSGFromMovieRatings(const long& max_user_id, const long& max_movie_id, const string& outputFileName);
	
private:
	// the movie_users structure maps a movie_id (long) to a vector of <user_id, rating> pairs
	std::vector< std::vector< std::pair<long, double> > > movie_users;
	// the sparse matrix with movie ratings
	compressed_matrix<double> star;

	std::string get_file_contents(const char *filename);
	void find_and_replace(string& source, string const& find, string const& replace);
};

} /* namespace generation */

#endif /* SRC_MOVIELENSSGCONVERTER_H_ */
