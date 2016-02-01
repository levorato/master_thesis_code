/*
 * MovieLensSGConverter.h
 *
 *  Created on: Feb 1, 2016
 *      Author: mlevorato
 */

#ifndef SRC_MOVIELENSSGCONVERTER_H_
#define SRC_MOVIELENSSGCONVERTER_H_

#include "../../include/ParallelILS.h"

#include <map>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;

namespace generation {

class MovieLensSGConverter {
public:
	MovieLensSGConverter();
	
	void processMovieLensFolder(const string& folder, const string& filter);
	void readMovieLensCSVFile(const string& filename);
	void generateSGFromMovieRatings();
	
private:
	// the movie_users structure maps a movie_id (long) to a vector of <user_id, rating> pairs
	std::map< long, std::vector< std::pair<long, double> > > movie_users;
	// the sparse matrix with movie ratings
	ublas::compressed_matrix<double> star;

	std::string get_file_contents(const char *filename);
	void find_and_replace(string& source, string const& find, string const& replace);
};

} /* namespace generation */

#endif /* SRC_MOVIELENSSGCONVERTER_H_ */
