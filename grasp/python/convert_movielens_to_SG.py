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
                        count = max_user_id * max_user_id
                        previous_total = -1
                        percentage = 0
                        join_df = pd.merge(df, df, on='movie_id', how='inner')
                        #ab_df = join_df.groupby(['user_id_x', 'user_id_y']).count()

                        for user_a in xrange(0, max_user_id):
                            for user_b in xrange(0, user_a):
                                #for idx, row2 in ab_df.iterrows():
                                #user_a = row2['user_id_x']
                                #user_b = row2['user_id_y']
                                if user_a != user_b:
                                    common_rating_count = 0
                                    common_similar_rating_count = 0
                                    #df_a = df[df.user_id == user_a]
                                    #df_b = df[df.user_id == user_b]
                                    #common_movies_df = pd.merge(df_a, df_b, on='movie_id', how='inner')
                                    #movies = common_movies_df['movie_id'] #, 'rating_x', 'rating_y']
                                    movies = join_df[join_df.user_id_x == user_a]
                                    movies = movies[movies.user_id_y == user_b]
                                    common_rating_count = movies['movie_id'].count()

                                    bad = movies[movies.rating_x <= BAD_MOVIE]
                                    bad = bad[movies.rating_y <= BAD_MOVIE]

                                    good = movies[movies.rating_x >= GOOD_MOVIE]
                                    good = good[movies.rating_y >= GOOD_MOVIE]

                                    regular = movies[movies.rating_x == REGULAR_MOVIE]
                                    regular = regular[movies.rating_y == REGULAR_MOVIE]
                                    common_similar_rating_count = bad['movie_id'].count() + good['movie_id'].count() + regular['movie_id'].count()

                                    # if len(movies) > 0:
                                    #     # print "common movies for users {0} and {1}: ".format(user_a, user_b)
                                    #     # print common_movies_df
                                    #
                                    #     #for movie_id in xrange(0, max_movie_id):
                                    #     #    if star[user_a, movie_id] != 0 and star[user_b, movie_id] != 0:
                                    #     for index, row in movies.iterrows(): #common_movies_df.iterrows():
                                    #         #print row
                                    #         # columns: user_id_x  movie_id  rating_x  timestamp_x  user_id_y  rating_y
                                    #         movie_id = row['movie_id']
                                    #         rating_a = row['rating_x']
                                    #         rating_b = row['rating_y']
                                    #         #print "{0} {1} {2}".format(movie_id, rating_a, rating_b),
                                    #         # users A and B have voted on the same movie
                                    #         common_rating_count += 1
                                    #         #if(star[user_a, movie_id] <= BAD_MOVIE and star[user_b, movie_id] <= BAD_MOVIE) \
                                    #         #        or star[user_a, movie_id] >= GOOD_MOVIE and star[user_b, movie_id] >= GOOD_MOVIE:
                                    #         if(rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE) \
                                    #                 or rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE \
                                    #                 or rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE:
                                    #             # users A and B have the same opinion about the movie
                                    #             common_similar_rating_count += 1
                                    #             #print "agree"
                                    #     # end for all movies
                                    if common_rating_count > 0:
                                        common_similar_rating_ratio = float(common_similar_rating_count) / common_rating_count
                                    else:
                                        common_similar_rating_ratio = 0

                                    if common_similar_rating_ratio >= POS_EDGE_PERC:
                                        #SG[user_a, user_b] = 1
                                        edge_list.append("{0} {1} 1\n".format(user_a, user_b))
                                    if common_similar_rating_ratio <= NEG_EDGE_PERC:
                                        #SG[user_a, user_b] = -1
                                        edge_list.append("{0} {1} -1\n".format(user_a, user_b))
                                # display status of processing done
                                total_done = user_a * user_b
                                threshold = int(math.floor(count / 5.0))
                                percentage = int(math.ceil(100 * (float(total_done) / count)))
                                #print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
                                if total_done % threshold < EPS and percentage != previous_total:
                                    print str(percentage) + " % ",
                                    previous_total = percentage

                        # writes the header info
                        output_file.write("{0}\t{1}\r\n".format(max_user_id, len(edge_list)))
                        for edge in edge_list:
                            output_file.write(edge)

                    finally:
                        content_file.close()
                        output_file.close()
                    print "Created output signed graph file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Graph generation took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file


        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
