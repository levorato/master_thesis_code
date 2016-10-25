import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as st

# gambit para python 2.4
class _minmax(object):
    def __init__(self, func):
        self.func = func

    def __call__(self, *seq, **kwargs):
        key = kwargs.pop('key', None)
        if kwargs:
            raise TypeError("only 'key' accepted as a "
                            "keyword argument")
        if key is None:
            return self.func(*seq)
        if len(seq) == 1:
            seq, = seq
        return self.func((key(item), i, item)
                         for i, item in enumerate(seq))[-1]


min = _minmax(min)
max = _minmax(max)


def q(cond, on_true, on_false):
    return {True: on_true, False: on_false}[cond is True]


# ---------------------------------------------------------
# natsort.py: Natural string sorting.
# ---------------------------------------------------------
def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    print [tryint(c) for c in re.split('([0-9]+)', s)][1]
    return int([tryint(c) for c in re.split('([0-9]+)', s)][1])


def natsorted(l):
    """ Sort the given list in the way that humans expect.
    """
    temp = [(alphanum_key(x), x) for x in l]
    temp.sort()


def main(argv):

    csv.field_size_limit(sys.maxsize)

    # Make the graphs a bit prettier, and bigger
    pd.set_option('display.mpl_style', 'default')
    pd.set_option('display.line_width', 5000)
    pd.set_option('display.max_columns', 60)

    parser = argparse.ArgumentParser(description='Process GRASP / ILS - CC and SRCC result files.')
    parser.add_argument('folder',
                        help='the folder containing the result files')
    parser.add_argument('--filefilter', default='*.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    args = parser.parse_args()
    folder = args.folder
    filter = args.filefilter

    args = parser.parse_args()

    print 'Input dir is ', folder
    print 'File filter is ', filter

    processCCResult(folder)


def processCCResult(folder):
    # CC results
    all_files_summary = dict()
    prefix_name = folder[:folder.rfind('/')]
    experiment_name = prefix_name[prefix_name.rfind('/') + 1:] + '-' + folder[folder.rfind('/') + 1:]
    previous_filename = ""
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        print "Processing folder " + ''.join(root)
        if (len(files) and ''.join(root) != folder):
            file_list = []
            file_list.extend(glob.glob(root + "/CC*.csv"))
            file_list.extend(glob.glob(root + "/Node*.csv"))
            count = len(file_list) - 1

            # Process CC results
            if os.path.isfile(root + "/cc-result.txt"):
                input_file = open(root + "/cc-result.txt", "r")

                content = input_file.read()
                reader = csv.reader(StringIO.StringIO(content), delimiter='\n')
                for row in reader:
                    if ''.join(row).find("time spent:") >= 0:
                        line = ''.join(row)
                        # the global time is the total time spent by the algorithm (including the parallel version)
                        global_time = float(line[line.find("time spent:") + 11:])
                        break
                input_file.close

                filename = (root[:root.rfind("/")])
                datetime = root[root.rfind("/") + 1:]
                filename = filename[filename.rfind("/") + 1:]
                if not all_files_summary.has_key(filename):
                    all_files_summary[filename] = []

                best_value = sys.float_info.max
                best_time = 0
                best_param = ''
                value = 0
                pos_value = 0
                neg_value = 0
                K = 0
                iteration = 0
                total_iter = 0
                total_comb = 0

                while count >= 0:
                    # Process all files from all metaheuristic executions, including those in parallel
                    print "Processing file " + file_list[count] + "\n"
                    content_file = open(file_list[count], 'r')
                    try:
                        content = content_file.read()

                        reader = csv.reader(StringIO.StringIO(content), delimiter=',')
                        linecount = 0
                        best_value = sys.float_info.max
                        iter_without_improvement = 0
                        iter_bestsol = 0
                        time_bestsol = float(0.0)

                        for row in reader:
                            if linecount == 0:
                                linecount += 1
                                continue
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if len(column) < 3:
                                linecount += 1
                                continue
                            value = float(column[1])
                            if value < best_value:
                                best_value = value
                                iter_bestsol = linecount
                                time_bestsol = float(column[5])
                                iter_without_improvement = 0
                            else:
                                iter_without_improvement += 1

                            if iter_without_improvement == 400:
                                filepath = ''.join(file_list[count])
                                # obtains the best result found by a specific execution of a specific node (can be parallel)
                                value = float(column[1])
                                pos_value = float(column[2])
                                neg_value = float(column[3])
                                K = long(column[4])
                                iteration = long(iter_bestsol)
                                best_time = time_bestsol
                                global_time = float(column[5])
                                total_iter = long(column[0])
                                # totalizes the number of visited solutions of all nodes running in parallel
                                total_comb = 0
                                break
                            linecount += 1
                            # next line
                        count = count - 1
                        # next csv file
                    finally:
                        content_file.close()
                all_files_summary[filename].append(str(filename) + "; " + str(datetime) + "; " + str(best_value) + "; " + str(pos_value) + "; " + str(
                        neg_value) + "; " + str(K) + "; " + str(iteration) + "; " + str(best_time) + "; " + str(
                        global_time) + "; " + best_param + "; " + str(total_iter) + "; " + str(total_comb))

            # end loop
            # process last file
            # end process CC results
    print "\nExporting consolidated CC Results to csv file...\n"


    # Save CC results of all executions of all instances to csv file: Time to reach best solution and total time spent
    result_file = open(folder + "/summary-400wi.csv", "w")
    #print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
    result_file.write(
            "Instance; ExecutionID; I(P); I(P)+; I(P)-; k; IterBestSol; TimeToBestSol; Global time(s); Params; Total Iter; NumVisitedSolutions\n")
    for key in sorted(all_files_summary.iterkeys()):
        #print "%s, %s" % (key, all_files_summary[key])
        execution_list = all_files_summary[key]
        for execution in execution_list:
            result_file.write("%s\n" % (execution))
    result_file.close()

    # re-reads the contents of summary.csv back to pandas dataframe
    df = pd.read_csv(folder + "/summary-400wi.csv", sep='; ', encoding="utf-8-sig")  # index_col='Instance',

    grouped_results = df.groupby('Instance')
    avg_results = grouped_results.agg([np.mean, np.median, np.max, np.count_nonzero, lambda x: (np.std(x, ddof=1)/np.sqrt(x.count())) * 1.96])  # , np.std
    print avg_results

    # Obtain mean of each group
    grouped = grouped_results['TimeToBestSol']
    means = grouped.mean()
    # print means

    # Calculate 95% confidence interval for each group => margin of error
    # http://www.wikihow.com/Calculate-Confidence-Interval
    ci_time_best_sol = grouped.aggregate(lambda x: (np.std(x, ddof=1)/np.sqrt(x.count())) * 1.96)
    # confidence interval = means +/- ci
    #print ci_time_best_sol
    #print means + " +/- " + ci

    # http://pandas.pydata.org/pandas-docs/version/0.15.1/generated/pandas.DataFrame.to_latex.html
    # print avg_results.to_latex()
    # http://pandas.pydata.org/pandas-docs/version/0.15.1/generated/pandas.DataFrame.to_csv.html?highlight=to_csv#pandas.DataFrame.to_csv
    avg_results.to_csv(str(experiment_name) + '-400wi.csv')

    #print best_results.agg([np.sum, np.mean, np.std])


if __name__ == "__main__":
    main(sys.argv[1:])
