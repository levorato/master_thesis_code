/*  MovieLens dataset: http://grouplens.org/datasets/movielens/
 *  All ratings are contained in the file "ratings.dat" and are in the
 *  following format:
 * 
 *  UserID::MovieID::Rating::Timestamp
 * 
 *  - UserIDs range between 1 and 6040
 *  - MovieIDs range between 1 and 3952
 *  - Ratings are made on a 5-star scale (whole-star ratings only)
 *  - Timestamp is represented in seconds since the epoch as returned by time(2)
 *  - Each user has at least 20 ratings
 */
#include <fstream>      // fstream
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include <iterator>     // ostream_operator

#include <boost/tokenizer.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Global variables
#define EPS = 0.0001
// Be sure to update these numbers as the movielens dataset grows!
#define MAX_MOVIES = 100000
#define MAX_USERS = 250000
// BAD_MOVIE -> maximum number of stars so that a rated movie is judged as a bad movie (e.g. 2 stars)
#define BAD_MOVIE = 2
// GOOD_MOVIE -> minimum number of stars so that a rated movie is judged as a good movie (e.g. 4 stars)
#define GOOD_MOVIE = 4
// REGULAR_MOVIE -> rating equals to 3
#define REGULAR_MOVIE = 3
// %POS -> edge percentual when comparing 2 users and assuming their relation is positive (e.g. 80%)
#define POS_EDGE_PERC = 0.8
// %NEG -> edge percentual when comparing 2 users and assuming their relation is negative (e.g. 20%)
#define NEG_EDGE_PERC = 0.2

namespace ublas = boost::numeric::ublas;

namespace generation {

using namespace resolution::construction;
using namespace boost::algorithm;
using namespace std;
using namespace util;
using namespace util::parallel;

MovieLensSGConverter::MovieLensSGConverter() {
	
}

MovieLensSGConverter::~MovieLensSGConverter() {
	// TODO Auto-generated destructor stub
}

MovieLensSGConverter::processMovieLensFolder(const string& folder, const string& filter) {
	cout << "Processing folder " << folder << endl;
	boost::timer::cpu_timer timer;
	boost::timer::cpu_times start_time, end_time;
	double timeSpent = 0.0;
	// measure conversion execution time
	timer.start();
	start_time = timer.elapsed();
	
	// read the files in the specified folder
	fs::path inputDir(folder);
	fs::directory_iterator end_iter;
	if (!fs::exists(inputDir) || !fs::is_directory(inputDir)) {
		cerr << "Input file directory not found. Exiting." << endl;
		return 1;
	}
	for( fs::directory_iterator dir_iter(inputDir) ; dir_iter != end_iter ; ++dir_iter) {
		string ext = dir_iter->path().extension().string();
		boost::algorithm::to_lower(ext);
		if ((fs::is_regular_file(dir_iter->status())) &&
				ext == filter) {
			fs::path filePath = *dir_iter;
			fileList.push_back(filePath);
		}
	}
	// output folder
	boost::filesystem::path outputPath(inputDir + boost::filesystem::path::preferred_separator + "unweightedSG")
	boost::system::error_code returnedError;
	boost::filesystem::create_directories( outputPath, returnedError );
	
	// read the movielens dataset voting file from CSV
	for(unsigned int i = 0; i < fileList.size(); i++) {
		fs::path filePath = fileList.at(i);
		readMovieLensCSVFile(filePath.string());
	}
	// process the dataset file and generate signed graph
	generateSGFromMovieRatings();
	
	cout << "Created output signed graph file " << outputFile << endl;
	// Stops the timer and stores the elapsed time
	timer.stop();
	end_time = timer.elapsed();
	timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
	cout << "Graph generation took " << timeSpent << " seconds." << endl;
}

void MovieLensSGConverter::readMovieLensCSVFile(const string& filename) {
	ifstream in(filename.c_str());
    if (!in.is_open()) return 1;

    typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	char_separator<char> sep("\r\n");
	Tokenizer< char_separator<char> > tokens(graphContents, sep);
	vector< string > lines;
	lines.assign(tokens.begin(),tokens.end());
	
	star = ublas::compressed_matrix<double>(MAX_USERS, MAX_MOVIES);
    
	// reads all votes for movies
	long max_user_id = 0;
	long max_movie_id = 0;
    for(std::vector<string>::iterator line_iter = lines.begin(); line_iter != lines.end(); line_iter++) {
		string line = *line_iter;
		find_and_replace(text, "::", ":");
        Tokenizer tok(line);
		vector< string > column;
        column.assign(tok.begin(),tok.end());
        if (column.size() < 3)  continue;

		long user_id = boost::lexical_cast< long >( column[0] );
        long movie_id = boost::lexical_cast< long >( column[1] );
		double rating = boost::lexical_cast< double >( column[2] );
		// user_id and movie_id both begin with 1
		// print "{0} {1}".format(user_id, movie_id)
		star(user_id - 1, movie_id - 1) = rating;
		movie_users[movie_id - 1].append(std::make_pair(user_id - 1, rating));
		// detects the number of users and movies in the dataset
		if(user_id > max_user_id)
			max_user_id = user_id;
		if(movie_id > max_movie_id)
			max_movie_id = movie_id;
    }
	cout << "Successfully read input file, generating signed graph file." << endl;
}

void MovieLensSGConverter::generateSGFromMovieRatings() {
	// the signed graph representing the voting in movie lens dataset
	int previous_total = -1;
	int percentage = 0;

	ublas::compressed_matrix<long> common_rating_count(max_user_id, max_user_id);
	ublas::compressed_matrix<long> common_similar_rating_count(max_user_id, max_user_id);

	cout << "Begin movie list traversal..." << endl;
	for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
		for( user_a, rating_a in movie_users[movie_id]) {
			for(user_b, rating_b in movie_users[movie_id]) {
				if(user_a != user_b) {
					common_rating_count(user_a, user_b) += 1;
					if((rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE) 
							or (rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE) 
							or (rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE)) {
						// users A and B have the same opinion about the movie
						common_similar_rating_count(user_a, user_b) += 1;
						// print "agree"
					}
				}
			}
		}
		// display status of processing done
		threshold = int(std::floor(max_movie_id / 10.0));
		percentage = int(std::ceil(100 * (double(movie_id) / max_movie_id)));
		if((movie_id % threshold < EPS) and (percentage != previous_total)) {
			cout << str(percentage) << " % ";
			previous_total = percentage;
		}
	}

	cout << "\nBegin edge generation...\n";
	long count = max_user_id * max_user_id;
	previous_total = -1;
	percentage = 0;
	std::ostringstream out;
	long edgeCount = 0;
	for(long user_a = 0; user_a < max_user_id; user_a++) {
		for(long user_b = 0; user_b < max_user_id; user_b++) {
			if(user_a != user_b) {
				if(common_rating_count(user_a, user_b) > 0) {
					common_similar_rating_ratio = double(common_similar_rating_count(user_a, user_b)) / common_rating_count(user_a, user_b);
					if(common_similar_rating_ratio >= POS_EDGE_PERC) {
						// SG[user_a, user_b] = 1;
						out << user_a << " " << user_b << " 1\n";
						edgeCount++;
					}
					else if(common_similar_rating_ratio <= NEG_EDGE_PERC) {
						// SG[user_a, user_b] = -1;
						out << user_a << " " << user_b << " -1\n";
						edgeCount++;
					}
				}
			}
			// display status of processing done
			total_done = user_a * max_user_id + user_b;
			threshold = int(std::floor(count / 10.0));
			percentage = int(std::ceil(100 * (double(total_done) / count)));
			// print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
			if((total_done % threshold < EPS) and (percentage != previous_total)) {
				cout << percentage << " % ";
				previous_total = percentage;
			}
		}
	}
	
	output_file << max_user_id << "\t" << edgeCount << "\r\n";
	output_file << out.str();
}

std::string MovieLensSGConverter::get_file_contents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return(contents.str());
  }
  throw(errno);
}

void MovieLensSGConverter::find_and_replace(string& source, string const& find, string const& replace)
{
    for(string::size_type i = 0; (i = source.find(find, i)) != string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length();
    }
}

} /* namespace generation */
