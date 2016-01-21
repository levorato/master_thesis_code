# MovieLens dataset: http://grouplens.org/datasets/movielens/
# All ratings are contained in the file "ratings.dat" and are in the
# following format:
#
# UserID::MovieID::Rating::Timestamp
#
# - UserIDs range between 1 and 6040
# - MovieIDs range between 1 and 3952
# - Ratings are made on a 5-star scale (whole-star ratings only)
# - Timestamp is represented in seconds since the epoch as returned by time(2)
# - Each user has at least 20 ratings

# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import argparse
import time
import math
from scipy.sparse import dok_matrix
import pandas as pd

# Global variables
EPS = 0.0001
# Be sure to update these numbers as the movielens dataset grows!
MAX_MOVIES = 35000
MAX_USERS = 250000
# BAD_MOVIE -> maximum number of stars so that a rated movie is judged as a bad movie (e.g. 2 stars)
BAD_MOVIE = 2
# GOOD_MOVIE -> minimum number of stars so that a rated movie is judged as a good movie (e.g. 4 stars)
GOOD_MOVIE = 4
# REGULAR_MOVIE -> rating equals to 3
REGULAR_MOVIE = 3
# %POS -> edge percentual when comparing 2 users and assuming their relation is positive (e.g. 80%)
POS_EDGE_PERC = 0.8
# %NEG -> edge percentual when comparing 2 users and assuming their relation is negative (e.g. 20%)
NEG_EDGE_PERC = 0.2


def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Convert MovieLens dataset file (ratings.dat) to unweighted signed graph.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the ratings.dat files')
    parser.add_argument('--filefilter', default='ratings.dat', required=False,
                        help='the filename for MovieLens ratings dataset files (default: ratings.dat)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter

    processMovieLensFiles(folders, filter)


def processMovieLensFiles(folders, filter):
    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = [f for f in files if filter in f]

                for file in file_list:
                    directory = root + os.sep + 'unweightedSG'
                    if not os.path.exists(directory):
                        os.makedirs(directory)

                    # measure elapsed time during graph generation
                    start = time.time()

                    file = os.path.join(root, file)
                    filename = file
                    print "Processing file " + filename
                    filename = filename[:filename.rfind(os.sep)]
                    filename = filename[filename.rfind(os.sep) + 1:]
                    content_file = open(file, 'r')
                    filename = "movielens-" + filename + ".g"
                    output_file = open(os.path.join(directory, filename), 'w')

                    # STAR(A,F) -> the number of stars (rating) user A gave to movie F
                    star = dok_matrix((MAX_USERS, MAX_MOVIES))
                    movie_users = [[] for x in xrange(0, MAX_MOVIES)]

                    try:
                        # detect the dialect (separator used in the csv file
                        dialect = csv.Sniffer().sniff(content_file.read(1024)) # , delimiters=";,"
                        content_file.seek(0)
                        reader = csv.reader(content_file, dialect)
                        max_user_id = 0
                        max_movie_id = 0

                        # reads all the votes for movies
                        for row in reader:
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if len(row) == 0:
                                continue

                            user_id = long(column[0])
                            movie_id = long(column[1])
                            rating = float(column[2])
                            # user_id and movie_id both begin with 1
                            star[user_id - 1, movie_id - 1] = rating
                            movie_users[movie_id - 1].append((user_id - 1, rating))

                            # detects the number of users and movies in the dataset
                            if user_id > max_user_id:
                                max_user_id = user_id
                            if movie_id > max_movie_id:
                                max_movie_id = movie_id
                        # end reading voting file
                        print "Successfully read input file, generating signed graph file."

                        # re-reads the contents of csv back to pandas dataframe
                        content_file2 = open(file, 'r')
                        df = pd.read_csv(content_file2, dialect=dialect, encoding="utf-8-sig",
                                         header=None, names=['user_id', 'movie_id', 'rating', 'timestamp'])  # index_col='Instance',
                        content_file2.close()

                        # the signed graph representing the voting in movie lens dataset
                        SG = dok_matrix((max_user_id, max_user_id))
                        edge_list = []
                        previous_total = -1
                        percentage = 0

                        common_rating_count = dok_matrix((max_user_id, max_user_id))
                        common_similar_rating_count = dok_matrix((max_user_id, max_user_id))

                        print "Begin movie list traversal..."
                        for movie_id in xrange(0, max_movie_id):
                            for user_a, rating_a in movie_users[movie_id]:
                                for user_b, rating_b in movie_users[movie_id]:
                                    if user_a != user_b:
                                        common_rating_count[user_a,user_b] += 1
                                        if(rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE) \
                                                or (rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE) \
                                                or (rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE):
                                            # users A and B have the same opinion about the movie
                                            common_similar_rating_count[user_a,user_b] += 1
                                            #print "agree"
                            # display status of processing done
                            threshold = int(math.floor(max_movie_id / 10.0))
                            percentage = int(math.ceil(100 * (float(movie_id) / max_movie_id)))
                            if movie_id % threshold < EPS and percentage != previous_total:
                                print str(percentage) + " % ",
                                previous_total = percentage

                        print "\nBegin edge generation..."
                        count = max_user_id * max_user_id
                        previous_total = -1
                        percentage = 0
                        for user_a in xrange(0, max_user_id):
                            for user_b in xrange(0, max_user_id):
                                if user_a != user_b:
                                    if common_rating_count[user_a,user_b] > 0:
                                        common_similar_rating_ratio = float(common_similar_rating_count[user_a,user_b]) / common_rating_count[user_a,user_b]
                                    else:
                                        common_similar_rating_ratio = 0

                                    if common_similar_rating_ratio >= POS_EDGE_PERC:
                                        #SG[user_a, user_b] = 1
                                        edge_list.append("{0} {1} 1\n".format(user_a, user_b))
                                    if common_similar_rating_ratio <= NEG_EDGE_PERC:
                                        #SG[user_a, user_b] = -1
                                        edge_list.append("{0} {1} -1\n".format(user_a, user_b))
                                # display status of processing done
                                total_done = user_a * max_user_id + user_b
                                threshold = int(math.floor(count / 10.0))
                                percentage = int(math.ceil(100 * (float(total_done) / count)))
                                #print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
                                if total_done % threshold < EPS and percentage != previous_total:
                                    print str(percentage) + " % ",
                                    previous_total = percentage

                        # for user_a in xrange(0, max_user_id):
                        #     for user_b in xrange(0, user_a):
                        #         if user_a != user_b:
                        #             common_rating_count = 0
                        #             common_similar_rating_count = 0
                        #             # print "common movies for users {0} and {1}: ".format(user_a, user_b)
                        #             # print common_movies_df
                        #
                        #             for movie_id in xrange(0, max_movie_id):
                        #                 if star[user_a, movie_id] != 0 and star[user_b, movie_id] != 0:
                        #                     #print row
                        #                     # users A and B have voted on the same movie
                        #                     common_rating_count += 1
                        #                     if(star[user_a, movie_id] <= BAD_MOVIE and star[user_b, movie_id] <= BAD_MOVIE) \
                        #                             or star[user_a, movie_id] >= GOOD_MOVIE and star[user_b, movie_id] >= GOOD_MOVIE \
                        #                             or star[user_a, movie_id] == REGULAR_MOVIE and star[user_b, movie_id] == REGULAR_MOVIE:
                        #                         # users A and B have the same opinion about the movie
                        #                         common_similar_rating_count += 1
                        #                         #print "agree"
                        #             # end for all movies
                        #             if common_rating_count > 0:
                        #                 common_similar_rating_ratio = float(common_similar_rating_count) / common_rating_count
                        #             else:
                        #                 common_similar_rating_ratio = 0
                        #
                        #             if common_similar_rating_ratio >= POS_EDGE_PERC:
                        #                 #SG[user_a, user_b] = 1
                        #                 edge_list.append("{0} {1} 1\n".format(user_a, user_b))
                        #             if common_similar_rating_ratio <= NEG_EDGE_PERC:
                        #                 #SG[user_a, user_b] = -1
                        #                 edge_list.append("{0} {1} -1\n".format(user_a, user_b))
                        #         # display status of processing done
                        #         total_done = user_a * user_b
                        #         threshold = int(math.floor(count / 10.0))
                        #         percentage = int(math.ceil(100 * (float(total_done) / count)))
                        #         #print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
                        #         if total_done % threshold < EPS and percentage != previous_total:
                        #             print str(percentage) + " % ",
                        #             previous_total = percentage

                        # writes the header info
                        output_file.write("{0}\t{1}\r\n".format(max_user_id, len(edge_list)))
                        for edge in edge_list:
                            output_file.write(edge)

                    finally:
                        content_file.close()
                        output_file.close()
                    print "\nCreated output signed graph file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Graph generation took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file


        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
