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

#include "MovieLensSGConverter.h"

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Global variables
#define EPS 0.0001
// Be sure to update these numbers as the movielens dataset grows!
#define MAX_MOVIES 100000
#define MAX_USERS 250000
// BAD_MOVIE -> maximum number of stars so that a rated movie is judged as a bad movie (e.g. 2 stars)
#define BAD_MOVIE 2
// GOOD_MOVIE -> minimum number of stars so that a rated movie is judged as a good movie (e.g. 4 stars)
#define GOOD_MOVIE 4
// REGULAR_MOVIE -> rating equals to 3
#define REGULAR_MOVIE 3
// %POS -> edge percentual when comparing 2 users and assuming their relation is positive (e.g. 80%)
#define POS_EDGE_PERC 0.8
// %NEG -> edge percentual when comparing 2 users and assuming their relation is negative (e.g. 20%)
#define NEG_EDGE_PERC 0.2
// NUMBER_OF_CHUNKS -> the number of chunks the matrix will be split for processing (saves memory)
#define NUMBER_OF_CHUNKS 256

const int generation::InputMessage::TERMINATE_MSG_TAG = 0;
const int generation::InputMessage::DONE_MSG_TAG = 1;

using namespace boost::numeric::ublas;
namespace fs = boost::filesystem;
namespace mpi = boost::mpi;

namespace generation {

using namespace boost::algorithm;
using namespace boost;
using namespace std;

MovieLensSGConverter::MovieLensSGConverter() {
	
}

MovieLensSGConverter::~MovieLensSGConverter() {
	// TODO Auto-generated destructor stub
}

bool MovieLensSGConverter::processMovieLensFolder(const string& folder, const string& filter,
		const unsigned int &myRank, const unsigned int &numProcessors) {
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
		file_ss << path_ss.str() << boost::filesystem::path::preferred_separator << filename << ".g";

		long max_user_id = 0, max_movie_id = 0;
		readMovieLensCSVFile(filePath.string(), max_user_id, max_movie_id);
		// process the dataset file and generate signed graph
		if(generateSGFromMovieRatings(max_user_id, max_movie_id, file_ss.str(), myRank, numProcessors)) {
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
		find_and_replace(line, "::", "\t");
        Tokenizer tok(line);
		std::vector< string > column;
        column.assign(tok.begin(),tok.end());
        if (column.size() < 3)  continue;

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
		const string& outputFileName, const unsigned int &myRank, const unsigned int &numProcessors) {
	// the signed graph representing the voting in movie lens dataset
	int previous_total = -1;
	int percentage = 0;

	// compressed_matrix<long> common_rating_count(max_user_id, max_user_id);
	// compressed_matrix<long> common_similar_rating_count(max_user_id, max_user_id);
	BOOST_LOG_TRIVIAL(info) << "Processing graph with " << max_user_id << " users and " << max_movie_id << " movies...";

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
	BOOST_LOG_TRIVIAL(info) << "This is process " << myRank << ". Will process user range [" << initialProcessorUserIndex << ", " << finalProcessorUserIndex << "]";

	// Begin dividing the MPI Chunk into smaller chunks for local processing
	long chunkSize = long(std::floor(double(MPIChunkSize) / NUMBER_OF_CHUNKS));
	long remainingVertices = MPIChunkSize % NUMBER_OF_CHUNKS;
	std::ostringstream out;
	long edgeCount = 0;
	for(int i = 0; i < NUMBER_OF_CHUNKS; i++) {
		// will process only the range (initialUserIndex <= user_a <= finalUserIndex)
		long initialUserIndex = initialProcessorUserIndex + i * chunkSize;
		long finalUserIndex = initialProcessorUserIndex + (i + 1) * chunkSize - 1;
		if(remainingVertices > 0 and i == NUMBER_OF_CHUNKS - 1) {
			finalUserIndex = finalProcessorUserIndex;
			chunkSize = finalUserIndex - initialUserIndex + 1;
		}
		BOOST_LOG_TRIVIAL(info) << "\nProcessing user range [" << initialUserIndex << ", " << finalUserIndex << "]";
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > > common_rating_count(chunkSize, max_user_id);;
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > > common_similar_rating_count(chunkSize, max_user_id);

		BOOST_LOG_TRIVIAL(info) << "Begin movie list traversal (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")...";
		cout << "\nBegin movie list traversal (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")..." << endl;
		for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
			std::vector< std::pair<long, char> > ratingList = movie_users[movie_id];
			// cout << movie_id << " with size " << ratingList.size() << endl;
			for(long x = 0; x < ratingList.size(); x++) {
				std::pair<long, char> user_rating_x = ratingList[x];
				long user_a = user_rating_x.first;
				char rating_a = user_rating_x.second;
				if(initialUserIndex <= user_a and user_a <= finalUserIndex) {
					for(long y = x + 1; y < ratingList.size(); y++) {
						std::pair<long, char> user_rating_y = ratingList[y];
						long user_b = user_rating_y.first;
						char rating_b = user_rating_y.second;
						if(user_a != user_b) {
							common_rating_count(user_a - initialUserIndex, user_b) += 1;
							if((rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE)
									or (rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE)
									or (rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE)) {
								// users A and B have the same opinion about the movie
								common_similar_rating_count(user_a - initialUserIndex, user_b) += 1;
							}
						}
					}
				}
			}
			// display status of processing done
			int threshold = int(std::floor(max_movie_id / 10.0));
			int percentage = int(std::ceil(100 * (double(movie_id) / max_movie_id)));
			if((movie_id % threshold < EPS) and (percentage != previous_total)) {
				cout << percentage << " % ";
				cout.flush();
				BOOST_LOG_TRIVIAL(info) << percentage << " % ";
				previous_total = percentage;
			}
		}

		BOOST_LOG_TRIVIAL(info) << "\nBegin edge generation (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")...";
		cout << "\nBegin edge generation (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")..." << endl;
		long count = (finalUserIndex - initialUserIndex + 1) * max_user_id;
		previous_total = -1;
		percentage = 0;
		for(long user_a = initialUserIndex; user_a <= finalUserIndex; user_a++) {
			for(long user_b = 0; user_b < max_user_id; user_b++) {
				if(user_a != user_b) {
					if(common_rating_count(user_a - initialUserIndex, user_b) > 0) {
						double common_similar_rating_ratio = double(common_similar_rating_count(user_a - initialUserIndex, user_b)) /
								common_rating_count(user_a - initialUserIndex, user_b);
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
	}

	// write the signed graph to output file
	stringstream partialFilename_ss;
	partialFilename_ss << outputFileName << ".part" << myRank;
	ofstream output_file(partialFilename_ss.str().c_str(), ios::out | ios::trunc);
	if(!output_file) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << outputFileName;
		return false;
	}
	if(myRank == 0) {
		output_file << max_user_id << "\t" << edgeCount << "\r\n";
	}
	output_file << out.str();
	// Close the file
	output_file.close();

	boost::mpi::communicator world;
	if(myRank == 0) {
		// wait for the message 'done' from all workers
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg;
			boost::mpi::status msg = world.recv(boost::mpi::any_source, InputMessage::DONE_MSG_TAG, imsg);
		}
		// merge all the output files into a full graph output file
		ofstream output_full_file(outputFileName.c_str(), ios::out | ios::trunc);
		if(!output_full_file) {
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open full output result file to: " << outputFileName;
			return false;
		}
		for(int i = 0; i < numProcessors; i++) {
			stringstream partialFilenameN_ss;
			partialFilenameN_ss << outputFileName << ".part" << i;
			ifstream in(partialFilenameN_ss.str().c_str());
			if (!in.is_open()){
				BOOST_LOG_TRIVIAL(fatal) << "Error opening output file number " << i;
				return false;
			}
			string fileContents = get_file_contents(partialFilenameN_ss.str().c_str());
			output_full_file << fileContents;
		}
		output_full_file.close();

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

} /* namespace generation */
