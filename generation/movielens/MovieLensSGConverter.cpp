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

#include "MovieLensSGConverter.h"

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Global variables
#define EPS 0.0000001
// Be sure to update these numbers as the movielens dataset grows!
#define MAX_MOVIES 100000
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
		const unsigned int &myRank, const unsigned int &numProcessors, const double& pos_edge_perc,
		const double& neg_edge_perc, const int& number_chunks) {
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
				"_p+" << std::fixed << std::setprecision(6) << pos_edge_perc << "_p-" << neg_edge_perc << ".g";

		long max_user_id = 0, max_movie_id = 0;
		readMovieLensCSVFile(filePath.string(), max_user_id, max_movie_id);
		// process the dataset file and generate signed graph
		if(generateSGFromMovieRatings(max_user_id, max_movie_id, file_ss.str(), myRank, numProcessors,
				pos_edge_perc, neg_edge_perc, number_chunks)) {
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
		const string& outputFileName, const unsigned int &myRank, const unsigned int &numProcessors,
		const double& pos_edge_perc, const double& neg_edge_perc, int number_chunks) {
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
	BOOST_LOG_TRIVIAL(info) << "This is process " << myRank << ". Will process user range [" <<
			initialProcessorUserIndex << ", " << finalProcessorUserIndex << "]";

	// Begin dividing the MPI Chunk into smaller chunks for local processing
	long chunkSize = long(std::floor(double(MPIChunkSize) / number_chunks));
	long remainingVertices = MPIChunkSize % number_chunks;
	std::ostringstream out;
	long edgeCount = 0;
	std::map<string, long> histogram;
	if(chunkSize == 0) {
		number_chunks = 1;
	}

	// determine the maximum common user rating count
	BOOST_LOG_TRIVIAL(info) << "1. Determining the maximum common_user_rating_count...";
	long max_common_rating_count = 0;
	for(int i = 0; i < number_chunks; i++) {
		// will process only the range (initialUserIndex <= user_a <= finalUserIndex)
		long initialUserIndex = initialProcessorUserIndex + i * chunkSize;
		long finalUserIndex = initialUserIndex + chunkSize - 1;
		if(remainingVertices > 0 and i == number_chunks - 1) {
			finalUserIndex = finalProcessorUserIndex;
			chunkSize = finalUserIndex - initialUserIndex + 1;
		}
		BOOST_LOG_TRIVIAL(info) << "1. Processing user range [" << initialUserIndex << ", " << finalUserIndex << "]";
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > >
				common_rating_count(chunkSize, max_user_id);;
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > >
				common_similar_rating_count(chunkSize, max_user_id);

		BOOST_LOG_TRIVIAL(info) << "1. Begin movie list traversal (step " << (i+1) << " of " << number_chunks << ")...";
		cout << "\n1. Begin movie list traversal (step " << (i+1) << " of " << number_chunks << ")..." << endl;
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
		// Discover the maximum common_rating_count between all pairs of users
		for(long user_a = initialUserIndex; user_a <= finalUserIndex; user_a++) {
			for(long user_b = 0; user_b < max_user_id; user_b++) {
				if(user_a != user_b) {
					long crc = common_rating_count(user_a - initialUserIndex, user_b);
					if(crc > 0 and crc > max_common_rating_count) {
						max_common_rating_count = crc;
					}
				}
			}
		}
	}
	boost::mpi::communicator world;
	if(myRank != 0) {  // every worker process sends its maximum count to the leader
		BOOST_LOG_TRIVIAL(info) << "1. Sending max_count to leader process. Local value = " << max_common_rating_count;
		InputMessage imsg(max_common_rating_count);
		world.send(0, InputMessage::MAX_COUNT_MSG_TAG, imsg);
		BOOST_LOG_TRIVIAL(info) << "1. Max count sent to leader process.";
		// receive the overall maximum common count
		boost::mpi::status msg = world.recv(0, InputMessage::MAX_COUNT_MSG_TAG, imsg);
		max_common_rating_count = imsg.max;
		BOOST_LOG_TRIVIAL(info) << "1. Overall max count received from leader. Step 1 done.";
	} else if(numProcessors > 1) {
		// wait for the maximum count from all workers
		BOOST_LOG_TRIVIAL(info) << "1. Waiting for common_user_rating_count from all worker processes...";
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg;
			boost::mpi::status msg = world.recv(boost::mpi::any_source, InputMessage::MAX_COUNT_MSG_TAG, imsg);
			if(imsg.max > max_common_rating_count) {
				max_common_rating_count = imsg.max;
			}
			BOOST_LOG_TRIVIAL(info) << "1. maximum common_user_rating_count received from process " << msg.source();
		}
		// send the overall maximum count back to all workers
		BOOST_LOG_TRIVIAL(info) << "1. (Leader) Sending the overall maximum common_user_rating_count to all worker processes...";
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg(max_common_rating_count);
			world.send(i, InputMessage::MAX_COUNT_MSG_TAG, imsg);
		}
		BOOST_LOG_TRIVIAL(info) << "1. (Leader) maximum common_user_rating_count sent to all worker processes. Step 1 done.";
	}
	BOOST_LOG_TRIVIAL(info) << "1. The maximum common count is " << max_common_rating_count;
	cout << "1. The maximum common count is " << max_common_rating_count;


	BOOST_LOG_TRIVIAL(info) << "2. Calculating all common_similar_rating_count for this process...";
	chunkSize = long(std::floor(double(MPIChunkSize) / number_chunks));
	remainingVertices = MPIChunkSize % number_chunks;
	if(chunkSize == 0) {
		number_chunks = 1;
	}
	edgeCount = 0;
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

		BOOST_LOG_TRIVIAL(info) << "Begin movie list traversal (step " << (i+1) << " of " << number_chunks << ")...";
		cout << "\nBegin movie list traversal (step " << (i+1) << " of " << number_chunks << ")..." << endl;
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

		BOOST_LOG_TRIVIAL(info) << "\nBegin edge generation (step " << (i+1) << " of " << number_chunks << ")...";
		cout << "\nBegin edge generation (step " << (i+1) << " of " << number_chunks << ")..." << endl;
		long count = (finalUserIndex - initialUserIndex + 1) * max_user_id;
		previous_total = -1;
		percentage = 0;

		// determine which pairs of users have the ratio (common_rating_count(a, b) / max_common_rating_count) bigger than a certain threshold
		for(long user_a = initialUserIndex; user_a <= finalUserIndex; user_a++) {
			for(long user_b = 0; user_b < max_user_id; user_b++) {
				if(user_a != user_b) {
					if(common_rating_count(user_a - initialUserIndex, user_b) > 0) {
						long n_a = common_similar_rating_count(user_a - initialUserIndex, user_b);
						// long n_b = common_rating_count(user_a - initialUserIndex, user_b);
						long n_b = max_common_rating_count;
						cpp_dec_float_50 common_similar_rating_ratio = cpp_dec_float_50(n_a);
						common_similar_rating_ratio /= cpp_dec_float_50(n_b);

						// converts the ratio (floating point number) to a string with 5-digit precision
						stringstream ss;
						ss << std::setprecision(5) << std::fixed << common_similar_rating_ratio;
						string key = ss.str();
						if(histogram.find(key) != histogram.end()) {
							histogram[key]++;
						} else {
							histogram[key] = 1;
						}

						/*
						double common_similar_rating_ratio = double(100 * common_similar_rating_count(user_a - initialUserIndex, user_b)) /
								common_rating_count(user_a - initialUserIndex, user_b); */
						cpp_dec_float_50 eps = std::numeric_limits<cpp_dec_float_50>::epsilon();
						/* TODO uncomment this
						if(((common_similar_rating_ratio - pos_edge_perc > eps) and (boost::multiprecision::abs(common_similar_rating_ratio - pos_edge_perc) > eps))  // (common_similar_rating_ratio > pos_edge_perc)
								or (boost::multiprecision::abs(common_similar_rating_ratio - pos_edge_perc) < eps)) {  // (common_similar_rating_ratio == pos_edge_perc)
							// SG[user_a, user_b] = 1;
							out << user_a << " " << user_b << " 1\n";
							edgeCount++;
						}
						else if(((common_similar_rating_ratio - neg_edge_perc < eps) and (boost::multiprecision::abs(common_similar_rating_ratio - neg_edge_perc) > eps))  // (common_similar_rating_ratio < neg_edge_perc)
								or (boost::multiprecision::abs(common_similar_rating_ratio - neg_edge_perc) < eps)) {  // (common_similar_rating_ratio == neg_edge_perc)
							// SG[user_a, user_b] = -1;
							out << user_a << " " << user_b << " -1\n";
							edgeCount++;
						} */
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

	// process the histogram values and output to text file
	stringstream ss_histogram;
	for(map<string,long>::const_iterator it = histogram.begin(); it != histogram.end(); ++it) {
	    ss_histogram << it->first << " " << it->second << "\n";
	}

	// write the signed graph to output file
	stringstream partialFilename_ss;
	partialFilename_ss << outputFileName << ".part" << myRank;
	ofstream output_file(partialFilename_ss.str().c_str(), ios::out | ios::trunc);
	if(!output_file) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << partialFilename_ss.str();
		return false;
	}
	if(myRank == 0) {
		output_file << max_user_id << "\t" << edgeCount << "\r\n";
	}
	output_file << out.str();
	// Close the SG file
	output_file.close();

	// write the histogram to a secondary output file
	stringstream partialFilenameHist_ss, outputFileNameHist_ss;
	partialFilenameHist_ss << outputFileName << "-histogram.part" << myRank;
	outputFileNameHist_ss << outputFileName << "-histogram.txt";
	ofstream output_file_hist(partialFilenameHist_ss.str().c_str(), ios::out | ios::trunc);
	if(!output_file_hist) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << partialFilenameHist_ss.str();
		return false;
	}
	output_file_hist << ss_histogram.str();
	output_file_hist.close();

	if(myRank == 0) {
		// wait for the message 'done' from all workers
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg;
			boost::mpi::status msg = world.recv(boost::mpi::any_source, InputMessage::DONE_MSG_TAG, imsg);
		}
		// merge all the signed graph output files into a full graph output file
		ofstream output_full_file_sg(outputFileName.c_str(), ios::out | ios::trunc);
		ofstream output_full_file_hist(outputFileNameHist_ss.str().c_str(), ios::out | ios::trunc);
		if(!output_full_file_sg) {
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open full SG output result file to: " << outputFileName;
			return false;
		}
		if(!output_full_file_hist) {
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open full histogram output result file to: " << outputFileNameHist_ss.str();
			return false;
		}
		for(int i = 0; i < numProcessors; i++) {
			stringstream partialSGFilenameN_ss, partialHistFilenameN_ss;
			partialSGFilenameN_ss << outputFileName << ".part" << i;
			partialHistFilenameN_ss << outputFileName << "-histogram.part" << i;
			ifstream inSG(partialSGFilenameN_ss.str().c_str());
			ifstream inHist(partialHistFilenameN_ss.str().c_str());
			if (!inSG.is_open() or !inHist.is_open()){
				BOOST_LOG_TRIVIAL(fatal) << "Error opening output file number " << i;
				return false;
			}
			string SGfileContents = get_file_contents(partialSGFilenameN_ss.str().c_str());
			output_full_file_sg << SGfileContents;

			if(i > 0) {  // the leader's histogram is already in the histogram map
				// merges the histogram data using a full map
				string line;
				getline(inHist, line);  // Get the frist line from the file, if any.
				while ( inHist ) {  // Continue if the line was sucessfully read.
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

					getline(inHist, line);   // Try to get another line.
				}
			}
		}
		// write the consolidated histogram to text file
		stringstream ss_histogram_full;
		for(map<string,long>::const_iterator it = histogram.begin(); it != histogram.end(); ++it) {
			ss_histogram_full << it->first << " " << it->second << "\n";
		}
		output_full_file_hist << ss_histogram_full.str();

		output_full_file_sg.close();
		output_full_file_hist.close();
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
