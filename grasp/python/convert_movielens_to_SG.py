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
                    #star = dok_matrix((MAX_USERS, MAX_MOVIES))

                    try:
                        # detect the dialect (separator used in the csv file
                        dialect = csv.Sniffer().sniff(content_file.read(1024)) #, delimiters=":: ")
                        max_user_id = 0

                        # reads the contents of csv to pandas dataframe
                        content_file2 = open(file, 'r')
                        text_content = content_file2.read()
                        text_content = text_content.replace("::", ":")
                        output = StringIO.StringIO(text_content)
                        df = pd.read_csv(output, dialect=dialect, encoding="utf-8-sig",
                                         header=None, names=['user_id', 'movie_id', 'rating', 'timestamp'])  # index_col='Instance',
                        content_file2.close()
                        print "Successfully read input file, generating signed graph file."

                        # the signed graph representing the voting in movie lens dataset
                        edge_list = []
                        join_df = pd.merge(df, df, on='movie_id', how='inner')
                        join_df = join_df[join_df.user_id_x != join_df.user_id_y]
                        join_df_g = join_df.groupby(['user_id_x', 'user_id_y']).size()
                        join_df_g = join_df_g.reset_index()
                        join_df_g.columns = join_df_g.columns.map(lambda x: 'count_cr' if "0" == str(x) else str(x))
                        #print join_df_g.columns
                        #print join_df_g

                        bad_df = join_df[join_df.rating_x <= BAD_MOVIE]
                        bad_df = bad_df[bad_df.rating_y <= BAD_MOVIE]
                        bad_df_g = bad_df.groupby(['user_id_x', 'user_id_y']).size()
                        bad_df_g = bad_df_g.reset_index()
                        bad_df_g.columns = bad_df_g.columns.map(lambda x: 'count_bad' if "0" == str(x) else str(x))

                        good_df = join_df[join_df.rating_x >= GOOD_MOVIE]
                        good_df = good_df[good_df.rating_y >= GOOD_MOVIE]
                        good_df_g = good_df.groupby(['user_id_x', 'user_id_y']).size()
                        good_df_g = good_df_g.reset_index()
                        good_df_g.columns = good_df_g.columns.map(lambda x: 'count_good' if "0" == str(x) else str(x))

                        regular_df = join_df[join_df.rating_x == REGULAR_MOVIE]
                        regular_df = regular_df[regular_df.rating_y == REGULAR_MOVIE]
                        regular_df_g = regular_df.groupby(['user_id_x', 'user_id_y']).size()
                        regular_df_g = regular_df_g.reset_index()
                        regular_df_g.columns = regular_df_g.columns.map(lambda x: 'count_reg' if "0" == str(x) else str(x))

                        result = pd.merge(join_df_g, bad_df_g, how='outer', on=['user_id_x', 'user_id_y'])
                        result = pd.merge(result, good_df_g, how='outer', on=['user_id_x', 'user_id_y'])
                        result = pd.merge(result, regular_df_g, how='outer', on=['user_id_x', 'user_id_y'])
                        #print result

                        count = result.count()['user_id_x']
                        previous_total = -1
                        percentage = 0
                        total_done = 0
                        for index, row in result.iterrows():
                            user_a = long(row['user_id_x'])
                            user_b = long(row['user_id_y'])
                            common_rating_count = row['count_cr']
                            common_similar_rating_count = row['count_bad'] + row['count_good'] + row['count_reg']

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
                            if user_a > max_user_id:
                                max_user_id = user_a
                            if user_b > max_user_id:
                                max_user_id = user_b
                            # display status of processing done
                            total_done += 1
                            threshold = int(math.floor(count / 20.0))
                            percentage = int(math.ceil(100 * (float(total_done) / count)))
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
                    print "\nCreated output signed graph file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Graph generation took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file


        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
