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
#include <iomanip>

#include <boost/tokenizer.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "MovieLensSGConverter.h"

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Global variables
#define EPS 0.0000001
// Be sure to update these numbers as the movielens dataset grows!
#define MAX_MOVIES 150000
#define MAX_USERS 250000
// BAD_MOVIE -> maximum number of stars so that a rated movie is judged as a bad movie (e.g. 2 stars)
#define BAD_MOVIE 2
// GOOD_MOVIE -> minimum number of stars so that a rated movie is judged as a good movie (e.g. 4 stars)
#define GOOD_MOVIE 4
// REGULAR_MOVIE -> rating equals to 3
#define REGULAR_MOVIE 3

// const int generation::InputMessage::TERMINATE_MSG_TAG = 0;
// const int generation::InputMessage::DONE_MSG_TAG = 1;

using namespace boost::numeric::ublas;
namespace fs = boost::filesystem;
namespace mpi = boost::mpi;

namespace generation {

using namespace boost::algorithm;
using namespace boost;
using namespace std;
using boost::multiprecision::cpp_dec_float_50;

MovieLensSGConverter::MovieLensSGConverter() {
	
}

MovieLensSGConverter::~MovieLensSGConverter() {
	// TODO Auto-generated destructor stub
}

bool MovieLensSGConverter::processMovieLensFolder(const string& folder, const string& filter,
		const unsigned int &myRank, const unsigned int &numProcessors, const double& minimum_edge_weight,
		const int& number_chunks) {
	BOOST_LOG_TRIVIAL(info) << "Processing folder " << folder;
	
	// read the files in the specified folder
	fs::path inputDir(folder);
	fs::directory_iterator end_iter;
	if (!fs::exists(inputDir) || !fs::is_directory(inputDir)) {
		BOOST_LOG_TRIVIAL(fatal) << "Input file directory not found. Exiting." << endl;
		return false;
	}
	std::vector<fs::path> fileList;
	for( fs::directory_iterator dir_iter(inputDir) ; dir_iter != end_iter ; ++dir_iter) {
		string ext = dir_iter->path().extension().string();
		string filename = dir_iter->path().filename().string();
		boost::algorithm::to_lower(filename);
		if ((fs::is_regular_file(dir_iter->status())) &&
				filename == filter) {
			fs::path filePath = *dir_iter;
			fileList.push_back(filePath);
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Found " << fileList.size() << " file(s).";
	// output folder
	stringstream path_ss;
	path_ss << folder << boost::filesystem::path::preferred_separator << "unweightedSG";
	boost::filesystem::path outputPath(path_ss.str());
	boost::system::error_code returnedError;
	boost::filesystem::create_directories( outputPath, returnedError );
	
	// read the movielens dataset voting file from CSV
	for(unsigned int i = 0; i < fileList.size(); i++) {
		boost::timer::cpu_timer timer;
		boost::timer::cpu_times start_time, end_time;
		double timeSpent = 0.0;
		// measure conversion execution time
		timer.start();
		start_time = timer.elapsed();

		fs::path filePath = fileList.at(i);
		string filename = filePath.parent_path().filename().string();
		stringstream file_ss;
		file_ss << path_ss.str() << boost::filesystem::path::preferred_separator << filename <<
				"_mw" << std::fixed << std::setprecision(4) << minimum_edge_weight << ".g";

		long max_user_id = 0, max_movie_id = 0;
		readMovieLensCSVFile(filePath.string(), max_user_id, max_movie_id);
		// process the dataset file and generate signed graph
		if(generateSGFromMovieRatings(max_user_id, max_movie_id, file_ss.str(), myRank, numProcessors,
				minimum_edge_weight, number_chunks)) {
			if(myRank == 0) {
				BOOST_LOG_TRIVIAL(info) << "\nCreated output signed graph file " << file_ss.str();
			}
			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			BOOST_LOG_TRIVIAL(info) << "Graph generation took " << timeSpent << " seconds.";
		} else {
			BOOST_LOG_TRIVIAL(fatal) << "Error generating signed graph file " << file_ss.str() << endl;
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Done.";
	return true;
}

bool MovieLensSGConverter::readMovieLensCSVFile(const string& filename, long& max_user_id, long& max_movie_id) {
	ifstream in(filename.c_str());
    if (!in.is_open()) return false;
    string graphContents = get_file_contents(filename.c_str());

    typedef tokenizer< char_separator<char> > Tokenizer;
	char_separator<char> sep("\r\n");
	Tokenizer tokens(graphContents, sep);
	std::vector< string > lines;
	lines.assign(tokens.begin(),tokens.end());
	
	star = compressed_matrix<double>(MAX_USERS, MAX_MOVIES);
	movie_users = std::vector< std::vector< std::pair<long, char> > >(MAX_MOVIES);
    
	// reads all votes for movies
	for(std::vector<string>::iterator line_iter = lines.begin(); line_iter != lines.end(); line_iter++) {
		string line = *line_iter;
		// BOOST_LOG_TRIVIAL(info) << "Reading line => " << line;
		find_and_replace(line, "::", "\t");
		find_and_replace(line, ",", "\t");
        Tokenizer tok(line);
		std::vector< string > column;
        column.assign(tok.begin(),tok.end());
        if (column.size() < 3)  continue;
        std::size_t found = column[0].find("user");
        if (found != std::string::npos)  continue;  // skips header line

		long user_id = boost::lexical_cast< long >( column[0] );
        long movie_id = boost::lexical_cast< long >( column[1] );
		char rating = boost::lexical_cast< int >( column[2] );
		// user_id and movie_id both begin with 1
		// print "{0} {1}".format(user_id, movie_id)
		star(user_id - 1, movie_id - 1) = rating;
		movie_users[movie_id - 1].push_back(std::make_pair(user_id - 1, rating));
		// detects the number of users and movies in the dataset
		if(user_id > max_user_id)
			max_user_id = user_id;
		if(movie_id > max_movie_id)
			max_movie_id = movie_id;
    }
	BOOST_LOG_TRIVIAL(info) << "Successfully read input file, generating signed graph file.";
	return true;
}

bool MovieLensSGConverter::generateSGFromMovieRatings(const long& max_user_id, const long& max_movie_id,
		const string& outputFileName, const unsigned int &myRank, const unsigned int &numProcessors,
		const double& minimum_edge_weight, int number_chunks) {
	// the signed graph representing the voting in movie lens dataset
	int previous_total = -1;
	int percentage = 0;

	// compressed_matrix<long> common_rating_count(max_user_id, max_user_id);
	// compressed_matrix<long> common_similar_rating_count(max_user_id, max_user_id);
	BOOST_LOG_TRIVIAL(info) << "Processing graph with " << max_user_id << " users and " << max_movie_id << " movies...";

	// Calculates the average rating of each user
	boost::numeric::ublas::vector<double> avg_user_rating(max_user_id, (double)0.0);
	boost::numeric::ublas::vector<int> user_rating_count(max_user_id, 0);
	boost::numeric::ublas::vector<long> user_rating_sum(max_user_id, 0);
	for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
		std::vector< std::pair<long, char> > ratingList = movie_users[movie_id];
		for(long x = 0; x < ratingList.size(); x++) {
			std::pair<long, char> user_rating_x = ratingList[x];
			long user = user_rating_x.first;
			char rating = user_rating_x.second;
			user_rating_count(user)++;
			user_rating_sum(user) += rating;
		}
	}
	for(long u = 0; u < max_user_id; u++) {
		avg_user_rating(u) = user_rating_sum(u) / (double) user_rating_count(u);
	}

	// Parallel processing with MPI *************************************
	// Splits the processing in (n / numProcessors) chunks,
	// to be consumed by numProcessors processes
	long MPIChunkSize = long(std::floor(double(max_user_id) / numProcessors));
	long MPIRemainingVertices = max_user_id % numProcessors;
	long initialProcessorUserIndex = myRank * MPIChunkSize;
	long finalProcessorUserIndex = (myRank + 1) * MPIChunkSize - 1;
	if(MPIRemainingVertices > 0 and myRank == numProcessors - 1) {
		finalProcessorUserIndex = max_user_id - 1;
		MPIChunkSize = finalProcessorUserIndex - initialProcessorUserIndex + 1;
	}
	BOOST_LOG_TRIVIAL(info) << "This is process " << myRank << ". Will process user range [" <<
			initialProcessorUserIndex << ", " << finalProcessorUserIndex << "]";

	// Begin dividing the MPI Chunk into smaller chunks for local processing
	long chunkSize = long(std::floor(double(MPIChunkSize) / number_chunks));
	long remainingVertices = MPIChunkSize % number_chunks;
	std::vector<string> outSG_cosine, outSG_pearson;
	long edgeCount = 0;
	std::map<string, long> histogram_cosine, histogram_pearson;
	if(chunkSize == 0) {
		number_chunks = 1;
	}
	boost::mpi::communicator world;

	// For each movie, calculates the variance of the ratings it has received
	std::vector<double> movie_rating_variance;
	for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
		std::vector< std::pair<long, char> > ratingList = movie_users[movie_id];
		if(ratingList.size() > 0) {
			// compute the average rating of this movie
			double average = double(0.0);
			for(long x = 0; x < ratingList.size(); x++) {
				std::pair<long, char>& user_rating_x = ratingList[x];
				char rating_a = user_rating_x.second;
				average += double(rating_a);
			}
			average /= ratingList.size();
			// compute the variance of this movie
			double sq_sum = double(0.0);
			for(long x = 0; x < ratingList.size(); x++) {
				std::pair<long, char>& user_rating_x = ratingList[x];
				char rating_a = user_rating_x.second;
				sq_sum += pow(double(rating_a) - average, 2);
			}
			movie_rating_variance.push_back(sq_sum / ratingList.size());
		} else {
			movie_rating_variance.push_back(double(0.0));
		}
	}

	double avg_rating_variance = double(0.0);
	for(long x = 0; x < movie_rating_variance.size(); x++) {
		avg_rating_variance += movie_rating_variance[x];
	}
	if(movie_rating_variance.size() > 0) {
		avg_rating_variance /= movie_rating_variance.size();
	}

	BOOST_LOG_TRIVIAL(info) << "2. Calculating the similarity measure between user votes for this process...";
	chunkSize = long(std::floor(double(MPIChunkSize) / number_chunks));
	remainingVertices = MPIChunkSize % number_chunks;
	if(chunkSize == 0) {
		number_chunks = 1;
	}
	for(int i = 0; i < number_chunks; i++) {
		// will process only the range (initialUserIndex <= user_a <= finalUserIndex)
		long initialUserIndex = initialProcessorUserIndex + i * chunkSize;
		long finalUserIndex = initialProcessorUserIndex + (i + 1) * chunkSize - 1;
		if(remainingVertices > 0 and i == number_chunks - 1) {
			finalUserIndex = finalProcessorUserIndex;
			chunkSize = finalUserIndex - initialUserIndex + 1;
		}
		BOOST_LOG_TRIVIAL(info) << "\nProcessing user range [" << initialUserIndex << ", " << finalUserIndex << "]";
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > >
				common_rating_count(chunkSize, max_user_id);;
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > >
				common_similar_rating_count(chunkSize, max_user_id);
		// for each pair of users, stores a list containing the common votes each user gave to a set of movies
		std::vector< std::map < long, CommonVoteList > > common_rating_list(chunkSize, std::map < long, CommonVoteList >());

		BOOST_LOG_TRIVIAL(info) << "1. Begin movie list traversal (step " << (i+1) << " of " << number_chunks << ")...";
		cout << "\n1. Begin movie list traversal (step " << (i+1) << " of " << number_chunks << ")..." << endl;
		for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
			std::vector< std::pair<long, char> > ratingList = movie_users[movie_id];
			for(long x = 0; x < ratingList.size(); x++) {
				std::pair<long, char> user_rating_x = ratingList[x];
				long user_a = user_rating_x.first;
				char rating_a = user_rating_x.second;
				if(initialUserIndex <= user_a and user_a <= finalUserIndex) {
					for(long y = 0; y < ratingList.size(); y++) {
						std::pair<long, char> user_rating_y = ratingList[y];
						long user_b = user_rating_y.first;
						char rating_b = user_rating_y.second;
						if(user_a != user_b) {
							if(common_rating_list[user_a - initialUserIndex].count(user_b) == 0) {
								common_rating_list[user_a - initialUserIndex][user_b] = CommonVoteList();
							}
							common_rating_list[user_a - initialUserIndex][user_b].push_back(CommonVote(rating_a, rating_b, movie_id));

							common_rating_count(user_a - initialUserIndex, user_b) += 1;
							/* if((rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE)
									or (rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE)
									or (rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE)) { */
							if(rating_a == rating_b) {
								// users A and B have the same opinion about the movie
								common_similar_rating_count(user_a - initialUserIndex, user_b) += 1;
							}
						}
					}
				}
			}
		}

		BOOST_LOG_TRIVIAL(info) << "Begin edge generation (step " << (i+1) << " of " << number_chunks << ")...";
		cout << "\nBegin edge generation (step " << (i+1) << " of " << number_chunks << ")..." << endl;
		long count = (finalUserIndex - initialUserIndex + 1) * max_user_id;
		previous_total = -1;
		percentage = 0;

		// determine which pairs of users have the ratio (common_rating_count(a, b) / max_common_rating_count) bigger than a certain threshold
		for(long user_a = initialUserIndex; user_a <= finalUserIndex; user_a++) {
			for(long user_b = 0; user_b < max_user_id; user_b++) {
				if(user_a < user_b) {
					if(common_rating_list[user_a - initialUserIndex].count(user_b) > 0) {
						CommonVoteList& voteList = common_rating_list[user_a - initialUserIndex][user_b];
						if(voteList.size() < 20)  continue;  // TODO considerar numero de votos em comum > x
						std::vector<double> normalizedVotesFromUserA, normalizedVotesFromUserB;
						std::vector<int> votesFromUserA, votesFromUserB;
						std::vector<long> common_movie_ids;
						for(CommonVoteList::iterator it = voteList.begin(); it != voteList.end(); it++) {
							votesFromUserA.push_back(it->rating_a);
							// normalizes the votes of each user by subtracting the rating from the average rating of the user
							double vote_a = it->rating_a - avg_user_rating(user_a - initialUserIndex);
							normalizedVotesFromUserA.push_back(vote_a);

							votesFromUserB.push_back(it->rating_b);
							double vote_b = it->rating_b - avg_user_rating(user_b - initialUserIndex);
							normalizedVotesFromUserB.push_back(vote_b);

							common_movie_ids.push_back(it->movie_id);
						}
						// calculates the cosine similarity and Pearson correlation coefficient of the vector
						// Warning: this operation can only be executed on non-zero arrays (v != (0, ..., 0)) !!!
						if(not (is_zero_array(normalizedVotesFromUserA) or is_zero_array(normalizedVotesFromUserB))) {
							double cosine_edge_weight = cosine_similarity(normalizedVotesFromUserA, normalizedVotesFromUserB);
							// double pearson_edge_weight = pearson_correlation_coefficient2(normalizedVotesFromUserA, normalizedVotesFromUserB);
							double pearson_edge_weight = pearson_correlation_coefficient_with_variance(normalizedVotesFromUserA, normalizedVotesFromUserB,
									common_movie_ids, movie_rating_variance, avg_rating_variance);
							// double spearman_edge_weight = spearman_correlation_coefficient(votesFromUserA, votesFromUserB);
							/*
							if(pearson_edge_weight > 1) {
								BOOST_LOG_TRIVIAL(error) << "Edge_weight = " << pearson_edge_weight.str(20);
								stringstream ss, ss2;
								for(int x = 0; x < normalizedVotesFromUserA.size(); x++) {
									ss << normalizedVotesFromUserA[x] << " (" << (int)voteList[x].rating_a << "); ";
									ss2 << normalizedVotesFromUserB[x] << "( " << (int)voteList[x].rating_b <<  "); ";
								}
								BOOST_LOG_TRIVIAL(error) << "Votes from user_a = " << ss.str();
								BOOST_LOG_TRIVIAL(error) << "Votes from user_b = " << ss2.str();
							} */
							// converts the weight (floating point number) to a string with 5-digit precision
							stringstream ss, ss2;
							ss << std::setprecision(3) << std::fixed << cosine_edge_weight;
							ss2 << std::setprecision(3) << std::fixed << pearson_edge_weight;
							string key = ss.str();
							if(histogram_cosine.find(key) != histogram_cosine.end()) {
								histogram_cosine[key]++;
							} else {
								histogram_cosine[key] = 1;
							}
							string key2 = ss2.str();
							if(histogram_pearson.find(key2) != histogram_pearson.end()) {
								histogram_pearson[key2]++;
							} else {
								histogram_pearson[key2] = 1;
							}
							// the matrix was sparsified by zeroing each entry smaller than minimum_edge_weight
							double eps = std::numeric_limits<double>::epsilon();
							double abs_cosine_edge_weight = fabs(cosine_edge_weight);
							std::ostringstream ss_edge;
							if((abs_cosine_edge_weight - minimum_edge_weight > eps) and (fabs(abs_cosine_edge_weight - minimum_edge_weight) > eps)) {
								// (abs_cosine_edge_weight > minimum_edge_weight)
								// SG[user_a, user_b] = 1;
								ss_edge << user_a << " " << user_b << " " << std::setprecision(5) << std::fixed << cosine_edge_weight << endl;
								outSG_cosine.push_back(ss_edge.str());
							}
							double abs_pearson_edge_weight = fabs(pearson_edge_weight);
							if((abs_pearson_edge_weight - minimum_edge_weight > eps) and (fabs(abs_pearson_edge_weight - minimum_edge_weight) > eps)) {
								// (abs_pearson_edge_weight > minimum_edge_weight)
								// SG[user_a, user_b] = 1;
								ss_edge << user_a << " " << user_b << " " << std::setprecision(5) << std::fixed << pearson_edge_weight << endl;
								outSG_pearson.push_back(ss_edge.str());
							}
						}
					}
				}
				// display status of processing done
				int total_done = (user_a - initialUserIndex) * max_user_id + user_b;
				int threshold = int(std::floor(count / 10.0));
				int percentage = int(std::ceil(100 * (double(total_done) / count)));
				// print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
				if((total_done % threshold < EPS) and (percentage != previous_total)) {
					cout << percentage << " % ";
					cout.flush();
					BOOST_LOG_TRIVIAL(info) << percentage << " % ";
					previous_total = percentage;
				}
			}
		}
		BOOST_LOG_TRIVIAL(info) << "Writing partial set of edges to text file...";
		// write the signed graphs to output file
		bool append = true;
		if(i == 0)  append = false;
		writeSignedGraphToFile("cosine", outputFileName, outSG_cosine, max_user_id, myRank, append);
		writeSignedGraphToFile("pearson", outputFileName, outSG_pearson, max_user_id, myRank, append);
		outSG_cosine.clear();
		outSG_pearson.clear();
	}

	BOOST_LOG_TRIVIAL(info) << "Writing histogram to text file...";
	// process the histogram values and output to text file
	stringstream ss_histogram_cosine;
	for(map<string,long>::const_iterator it = histogram_cosine.begin(); it != histogram_cosine.end(); ++it) {
	    ss_histogram_cosine << it->first << " " << it->second << "\n";
	}
	stringstream ss_histogram_pearson;
	for(map<string,long>::const_iterator it = histogram_pearson.begin(); it != histogram_pearson.end(); ++it) {
		ss_histogram_pearson << it->first << " " << it->second << "\n";
	}

	// write the histogram to a secondary output file
	writeHistogramToFile("cosine", outputFileName, ss_histogram_cosine.str(), myRank);
	writeHistogramToFile("pearson", outputFileName, ss_histogram_pearson.str(), myRank);
	stringstream outputFileNameHist_ss;
	outputFileNameHist_ss << outputFileName << "-histogram.txt";

	if(myRank == 0) {
		// wait for the message 'done' from all workers
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg;
			boost::mpi::status msg = world.recv(boost::mpi::any_source, InputMessage::DONE_MSG_TAG, imsg);
		}
		// merge all the signed graph output files into a full graph output file
		mergeSignedGraphToFile("cosine", outputFileName, max_user_id, numProcessors);
		mergeSignedGraphToFile("pearson", outputFileName, max_user_id, numProcessors);

		mergeHistogramToFile("cosine", outputFileName, numProcessors, histogram_cosine);
		mergeHistogramToFile("pearson", outputFileName, numProcessors, histogram_pearson);

	} else {
		// send a message to the leader process to inform that the task is done
		BOOST_LOG_TRIVIAL(info) << "Sending done message to leader...";
		InputMessage imsg;
		world.send(0, InputMessage::DONE_MSG_TAG, imsg);
	}

	return true;
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

double MovieLensSGConverter::cosine_similarity(std::vector<double>& votesFromUserA,
		std::vector<double>& votesFromUserB)
{
	double dot(0.0);
	double denom_a(0.0);
	double denom_b(0.0);
	unsigned long len = votesFromUserA.size();
	assert(len == votesFromUserB.size());
	if(len == 0)  return double(0.0);
	for(unsigned long i = 0u; i < len; i++) {
		dot += double(votesFromUserA[i]) * votesFromUserB[i];
		denom_a += double(votesFromUserA[i]) * votesFromUserA[i];
		denom_b += double(votesFromUserB[i]) * votesFromUserB[i];
	}
	return dot / (sqrt(denom_a * denom_b));
}

// the following function is not working properly when array values of the same user are equal FIXME
double MovieLensSGConverter::pearson_correlation_coefficient(std::vector<double>& votesFromUserA,
		std::vector<double>& votesFromUserB)
{
	// summing up the product
	double prod_sum(0.0);
	// summing the array items
	double sum_a(0.0);
	double sum_b(0.0);
	// summing up the squares
	double sum_a_sq(0.0);
	double sum_b_sq(0.0);
	unsigned long len = votesFromUserA.size();
	assert(len == votesFromUserB.size());
	if(len == 0)  return double(0.0);

	for(unsigned long i = 0u; i < len; i++) {
		sum_a += double(votesFromUserA[i]);
		sum_b += double(votesFromUserB[i]);
		sum_a_sq += double(votesFromUserA[i]) * votesFromUserA[i];
		sum_b_sq += double(votesFromUserB[i]) * votesFromUserB[i];
		prod_sum += double(votesFromUserA[i]) * votesFromUserB[i];
	}
	double num(0.0), den(0.0);
	num = prod_sum - (sum_a * sum_b / len);
	den = sqrt((sum_a_sq - (sum_a * sum_a) / len) * (sum_b_sq - (sum_b * sum_b) / len));
	if(den == 0)  return 0;
	return num / den;
}

// http://projekter.aau.dk/projekter/files/32181941/Report.pdf - pages 24/25
double MovieLensSGConverter::pearson_correlation_coefficient2(std::vector<double>& votesFromUserA,
		std::vector<double>& votesFromUserB)
{
	// summing up the product
	double prod_sum(0.0);
	// summing the array items
	double sum_a(0.0);
	double sum_b(0.0);
	// summing up the squares
	double sum_a_sq(0.0);
	double sum_b_sq(0.0);
	unsigned long len = votesFromUserA.size();
	assert(len == votesFromUserB.size());
	if(len == 0)  return double(0.0);

	for(unsigned long i = 0u; i < len; i++) {
		sum_a += double(votesFromUserA[i]);
		sum_b += double(votesFromUserB[i]);
		sum_a_sq += double(votesFromUserA[i]) * votesFromUserA[i];
		sum_b_sq += double(votesFromUserB[i]) * votesFromUserB[i];
		prod_sum += double(votesFromUserA[i]) * votesFromUserB[i];
	}
	double num(0.0), den(0.0);
	num = prod_sum;
	den = sqrt((sum_a_sq) * (sum_b_sq));
	if(den == 0)  return 0;
	return num / den;
}

// http://projekter.aau.dk/projekter/files/32181941/Report.pdf - pages 29/30/31
double MovieLensSGConverter::pearson_correlation_coefficient_with_variance(std::vector<double>& votesFromUserA,
		std::vector<double>& votesFromUserB, std::vector<long>& common_movie_ids, std::vector<double>& movie_rating_variance,
		const double& avg_rating_variance)
{
	// summing up the product
	double prod_sum(0.0);
	// summing up the squares
	double sum_a_sq(0.0);
	double sum_b_sq(0.0);
	// average variance of the ratings between users a and b
	double avg_variance(0.0);
	unsigned long len = votesFromUserA.size();
	assert(len == votesFromUserB.size());
	assert(len == common_movie_ids.size());
	if(len == 0)  return double(0.0);

	for(unsigned long i = 0u; i < len; i++) {
		long movie_id = common_movie_ids[i];
		double movie_i_variance(movie_rating_variance[movie_id]);
		// den
		sum_a_sq += double(votesFromUserA[i]) * votesFromUserA[i] * movie_i_variance;
		sum_b_sq += double(votesFromUserB[i]) * votesFromUserB[i] * movie_i_variance;
		// num
		prod_sum += double(votesFromUserA[i]) * votesFromUserB[i] * movie_i_variance;
	}
	double num(0.0), den(0.0);
	if(avg_rating_variance == 0)  return double(0.0);
	num = prod_sum / avg_rating_variance;
	sum_a_sq /= avg_rating_variance;
	sum_b_sq /= avg_rating_variance;
	den = sqrt((sum_a_sq) * (sum_b_sq));
	if(den == 0)  return double(0.0);

	double result = num / den;
	if(isnan(result)) {
		BOOST_LOG_TRIVIAL(error) << "avg_variance = " << avg_rating_variance;
		BOOST_LOG_TRIVIAL(error) << "num = " << num;
		BOOST_LOG_TRIVIAL(error) << "den = " << den;
	}
	return result;
}

double MovieLensSGConverter::spearman_correlation_coefficient(std::vector<int>& votesFromUserA,
		std::vector<int>& votesFromUserB)
{
	unsigned long len = votesFromUserA.size();
	assert(len == votesFromUserB.size());
	double n(len);
	double sum(0.0);
	if(len == 0)  return double(0.0);

	for(unsigned long i = 0u; i < len; i++) {
		// computes the sum of squared differences between the pairs of ranks (values)
		double dif = double(votesFromUserA[i]) - double(votesFromUserB[i]);
		sum += dif * dif;
	}
	return double(1) - (double(6) * sum) / (n * (n * n - 1));
}

bool MovieLensSGConverter::is_zero_array(std::vector<double>& array)
{
	for(std::vector<double>::iterator it = array.begin(); it != array.end(); it++) {
		if(fabs(*it) > EPS) {
			return false;
		}
	}
	return true;
}

bool MovieLensSGConverter::writeSignedGraphToFile(const string& file_prefix, const string& partialFilename,
		const std::vector<string>& edgeList, const long& max_user_id, const int& myRank, const bool& append) {
	stringstream ss;
	ss << partialFilename << "-" << file_prefix << ".part" << myRank;
	ofstream output_file(ss.str().c_str(), ios::out | (append ? ios::app : ios::trunc));
	if(!output_file) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << partialFilename;
		return false;
	}
	for(std::vector<string>::const_iterator it = edgeList.begin(); it != edgeList.end(); it++) {
		output_file << (*it);
	}
	// Close the SG file
	output_file.close();
	return true;
}

bool MovieLensSGConverter::writeHistogramToFile(const string& file_prefix, const string& partialFilename,
		string fileContent, const int& myRank) {
	stringstream ss;
	ss << partialFilename << "-" << file_prefix << "-histogram.part" << myRank;
	ofstream output_file_hist(ss.str().c_str(), ios::out | ios::trunc);
	if(!output_file_hist) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << ss.str();
		return false;
	}
	output_file_hist << fileContent;
	output_file_hist.close();
	return true;
}

bool MovieLensSGConverter::mergeSignedGraphToFile(const string& file_prefix, const string& partialFilename,
		const long& max_user_id, const int& numProcessors) {
	stringstream ss;
	ss << partialFilename << "-" << file_prefix << ".g";
	ofstream output_full_file_sg(ss.str().c_str(), ios::out | ios::trunc);
	if(!output_full_file_sg) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open full SG output result file to: " << ss.str();
		return false;
	}
	int edgeCount = 0;
	for(int i = 0; i < numProcessors; i++) {
		stringstream partialSGFilenameN_ss;
		partialSGFilenameN_ss << partialFilename << "-" << file_prefix << ".part" << i;
		edgeCount += countLinesInFile(partialSGFilenameN_ss.str());
	}
	output_full_file_sg << max_user_id << "\t" << edgeCount << "\r\n";
	for(int i = 0; i < numProcessors; i++) {
		stringstream partialSGFilenameN_ss;
		partialSGFilenameN_ss << partialFilename << "-" << file_prefix << ".part" << i;
		ifstream inSG(partialSGFilenameN_ss.str().c_str());
		if (!inSG.is_open()){
			BOOST_LOG_TRIVIAL(fatal) << "Error opening output file number " << i;
			return false;
		}
		string SGfileContents = get_file_contents(partialSGFilenameN_ss.str().c_str());
		output_full_file_sg << SGfileContents;
	}
	output_full_file_sg.close();
	return true;
}

bool MovieLensSGConverter::mergeHistogramToFile(const string& file_prefix, const string& partialFilename,
		const int& numProcessors, std::map<string, long>& histogram) {
	stringstream ss, ss2;
	ss << partialFilename << "-" << file_prefix << "-histogram.txt";
	ss2 << partialFilename << "-" << file_prefix << "-abs_histogram.txt";
	ofstream output_full_file_hist(ss.str().c_str(), ios::out | ios::trunc);
	ofstream output_full_file_abs_hist(ss2.str().c_str(), ios::out | ios::trunc);
	if(!output_full_file_hist or !output_full_file_abs_hist) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open full histogram output result file to: " << ss.str();
		return false;
	}
	// histogram with absolute edge weights (unsigned)
	std::map<string, long> abs_histogram;
	for(int i = 0; i < numProcessors; i++) {
		stringstream partialHistFilenameN_ss;
		partialHistFilenameN_ss << partialFilename << "-" << file_prefix << "-histogram.part" << i;
		ifstream inHist(partialHistFilenameN_ss.str().c_str());
		if (!inHist.is_open()){
			BOOST_LOG_TRIVIAL(fatal) << "Error opening output file number " << i;
			return false;
		}
		// merges the histogram data using a full map
		string line;
		getline(inHist, line);  // Get the first line from the file, if any.
		while ( inHist ) {  // Continue if the line was successfully read.
			// Process the line.
			istringstream iss(line);
			string key;
			long value;
			iss >> key >> value;
			if(histogram.find(key) != histogram.end()) {
				histogram[key] += value;
			} else {
				histogram[key] = value;
			}
			string key2 = key;
			if(key.substr(0, 1) == "-") {  // remove trailing minus from edge weights
				key2 = key.substr(1);
			}
			if(abs_histogram.find(key2) != abs_histogram.end()) {
				abs_histogram[key2] += value;
			} else {
				abs_histogram[key2] = value;
			}

			getline(inHist, line);   // Try to get another line.
		}
	}
	// write the consolidated histogram to text file
	stringstream ss_histogram_full, ss_abs_histogram_full;
	for(map<string,long>::const_iterator it = histogram.begin(); it != histogram.end(); ++it) {
		ss_histogram_full << it->first << " " << it->second << "\n";
	}
	for(map<string,long>::const_iterator it = abs_histogram.begin(); it != abs_histogram.end(); ++it) {
		output_full_file_abs_hist << it->first << " " << it->second << "\n";
	}
	output_full_file_hist << ss_histogram_full.str();
	output_full_file_hist.close();
	output_full_file_abs_hist << ss_abs_histogram_full.str();
	output_full_file_abs_hist.close();
	return true;
}

unsigned int MovieLensSGConverter::FileRead( istream & is, std::vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

unsigned int MovieLensSGConverter::CountLines( const std::vector <char> & buff, int sz ) {
    int newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
    	if ( p[i] == '\n' ) {
    		newlines++;
    	}
    }
    return newlines;
}

unsigned int MovieLensSGConverter::countLinesInFile(const string& filename) {
	const int SZ = 1024 * 1024;
	std::vector <char> buff( SZ );
	ifstream ifs( filename.c_str() );
	int n = 0;
	while( int cc = FileRead( ifs, buff ) ) {
		n += CountLines( buff, cc );
	}
	return n;
}

} /* namespace generation */
