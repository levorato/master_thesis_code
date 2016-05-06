#!/usr/bin/python

# Find out if we're using the set built-in type (Python > 2.3)
# or the Set module (Python <= 2.3)
try:
    set
except NameError:
    from sets import Set as set
from itertools import *
from optparse import OptionParser
from ConfigParser import ConfigParser
import os
import os.path
import argparse
import time
import sys
import subprocess
import traceback
import csv
import StringIO
import time
import glob
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main(argv):
    #csv.field_size_limit(sys.maxsize)
    parser = argparse.ArgumentParser(description='Process Pajek results.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the Pajek report and solution files (txt) format')
    parser.add_argument('--instancepath', required=True,
                        help='the path for graph instance files in Pajek .net format')

    args = parser.parse_args()
    folders = args.folders
    instance_path = args.instancepath

    args = parser.parse_args()

    print 'Pajek result folders are ', folders
    print 'Graph instance path is ', instance_path

    process_result_files(instance_path, folders)


def read_graph_file_pajek_format(filename):
    print 'Reading graph file: {0}'.format(str(filename))
    with open(filename, 'rb') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=" \t")
        csvfile.seek(0)
        r = csv.reader(csvfile, dialect)
        line1=r.next()
        values = [col for col in line1]
        N = int(values[1].strip())
        while '*Vertices' not in values[0]:
            line1 = r.next()
            values = [col for col in line1]
            N = int(values[1].strip())
        print 'N = {0}'.format(str(N))
        while '*Arcs' not in values[0]:
            line1 = r.next()
            values = [col for col in line1]
        matrix = [[0 for x in xrange(N)] for x in xrange(N)]
        for line in r:
            i = int(line[0]) - 1
            j = int(line[1]) - 1
            w = float(line[2])
            matrix[i][j] = w
            #print '{0} {1} = {2}'.format(str(i), str(j), str(w))
        print 'Successfully read instance file.'
    return N, matrix


def CC_objective_function(N, matrix, mycluster):
    positive_sum = float(0)
    negative_sum = float(0)

    # For each vertex i
    for i in xrange(N):
        # Find out to which cluster vertex i belongs
        k = mycluster[i]
        # For each out edge of i
        for j in xrange(N):
            if matrix[i][j] <> 0:
                weight = float(matrix[i][j])
                same_cluster = (k == mycluster[j])
                if weight < 0 and same_cluster:  # negative edge
                    # i and j are in the same cluster
                    negative_sum += weight * (-1)
                elif weight > 0 and (not same_cluster):  # positive edge
                    # i and j are NOT in the same cluster
                    positive_sum += weight

    print 'CC calculated value. I(P) = {0}\n'.format(str(positive_sum + negative_sum))
    return positive_sum, negative_sum


def extract_solution_file(result_filepath):
    print 'Reading solution file: {0}'.format(str(result_filepath))
    with open(result_filepath, 'rb') as f:
        content = f.readlines()
        i = 0
        max_c = 0
        for line in content:
            if i == 0:
                N = int(line[line.rfind(' ') + 1:])
                mycluster = [0 for x in xrange(N)]
                print 'N = {0}'.format(str(N))
            else:
                c = int(line) - 1
                mycluster[i - 1] = c
                if c + 1 > max_c:
                    max_c = c + 1
                #print 'Vertex {0} is in cluster {1}'.format(str(i), str(c))
            i += 1
        print 'Successfully read solution file: {0} clusters found.'.format(str(max_c))
    return max_c, mycluster


def extract_report_file(report_filepath):
    # Example of Pajek solution report file
    # ==============================================================================
    # Partitioning Signed Networks according to Structural Balance
    # ==============================================================================
    #  Working...
    #  Number of clusters: 7, alpha: 0.500, min size of clusters: 1
    # -------  Starting partition -------
    # Errors:     368.00       Lines
    # -----------------------------------
    # -------     Improvements    -------
    #      1:     112.50
    #      2:      99.00
    #      5:      55.00
    #      6:      49.00
    #      7:      34.50
    #     14:      26.50
    #     54:      26.00
    #    228:      24.00
    #  133 solutions with 24.00 inconsistencies found.
    #  Time spent:  0:02:49
    print 'Reading report file: {0}'.format(str(report_filepath))
    with open(report_filepath, 'rb') as file:
        for line in reversed(file.readlines()):
            if line.find('Time spent') > 0:
                line = line.strip()
                time_string = line[line.find(':') + 1:].strip()
                print time_string
                time_spent_obj = datetime.strptime(time_string, '%H:%M:%S')
                #print time_spent_obj
                time_spent = (time_spent_obj - time_spent_obj.replace(
                        hour=0, minute=0, second=0, microsecond=0)).total_seconds()
                print 'Time spent: {0} s'.format(str(time_spent))
                break
        print 'Successfully read report file.'
    return time_spent


def extract_instance_name(instance_name_list, filename):
    for name in instance_name_list:
        if name in filename:
            return name
    return None


def process_result_files(instance_path, folders):
    # obtain the list of instance files in instance_path
    instance_file_list = []
    instance_file_list.extend(glob.glob(os.path.join(instance_path, "*.net")))
    instance_name_list = [s[s.rfind(os.path.sep)+1:s.rfind('.')] for s in instance_file_list]
    print "List of instance names: {0}".format(str(instance_name_list))
    
    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = []
                file_list.extend(glob.glob(os.path.join(root, "*.rep")))
                #file_list.extend(glob.glob(root + "/*.txt"))
                count = len(file_list) - 1
                # CC results
                all_files_summary = dict()
                prefix_name = folder[:folder.rfind(os.path.sep)]
                experiment_name = prefix_name[prefix_name.rfind(os.path.sep) + 1:] + '-' + folder[folder.rfind(os.path.sep) + 1:]
                previous_filename = ""
                count = 0

                for report_file in file_list:
                    print "\nProcessing file " + report_file
                    base_path = report_file[:report_file.rfind(os.path.sep) + 1]
                    filename = report_file[report_file.rfind(os.path.sep) + 1:]
                    solution_file = report_file
                    solution_file = solution_file.replace('.rep', '.txt')
                    instance_name = extract_instance_name(instance_name_list, filename)
                    if instance_name is None:
                        print "Error processing report file {0}. Skipping file.".format(report_file)
                        continue
                    print "Instance name is " + instance_name
                    instance_file = os.path.join(instance_path, instance_name + '.net')
                    print "Instance file path is " + instance_file

                    if not all_files_summary.has_key(filename):
                        all_files_summary[filename] = []

                    pos_value = 0
                    neg_value = 0
                    k = 0
                    time_spent = 0
                    try:
                        if previous_filename != instance_name:
                            N, matrix = read_graph_file_pajek_format(instance_file)
                        k, mycluster = extract_solution_file(solution_file)
                        pos_value, neg_value = CC_objective_function(N, matrix, mycluster)
                        time_spent = extract_report_file(report_file)

                        all_files_summary[filename].append(str(instance_name) + "; " + str(count) + "; " + str(pos_value+neg_value) + "; "
                                                           + str(pos_value) + "; " + str(neg_value) + "; " + str(k) + "; "
                                                           + str(time_spent))

                    except Exception,e:
                        print "Error processing report file {0}".format(report_file)
                        print str(e)
                        traceback.print_exc()
                        continue

                    previous_filename = instance_name
                    count += 1
                # end loop
                # process last file
                print "\nExporting consolidated Pajek results to csv file..."

                # Save CC results of all executions of all instances to csv file: Time to reach best solution and total time spent
                result_file = open(os.path.join(folder, "summary.csv"), "w")
                result_file.write("Instance; ExecutionID; I(P); I(P)+; I(P)-; k; Time(s)\n")
                for key in sorted(all_files_summary.iterkeys()):
                    #print "%s, %s" % (key, all_files_summary[key])
                    execution_list = all_files_summary[key]
                    for execution in execution_list:
                        result_file.write("%s\n" % (execution))
                result_file.close()
                print "Done.\n"

                # re-reads the contents of summary.csv back to pandas dataframe
                df = pd.read_csv(os.path.join(folder, "summary.csv"), sep='; ', encoding="utf-8-sig")  # index_col='Instance',

                grouped_results = df.groupby('Instance')
                avg_results = grouped_results.agg([np.mean, np.median, np.max, lambda x: (np.std(x, ddof=1)/np.sqrt(x.count())) * 1.96])  # , np.std
                print "Aggregated results:\n"
                print avg_results
                avg_results.to_csv('aggregated.csv')

        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])

