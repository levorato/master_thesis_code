# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import argparse
import pandas as pd
import scipy.stats as stats
import numpy as np
from matplotlib import pyplot as plt

def main(argv):

    csv.field_size_limit(sys.maxsize)

    parser = argparse.ArgumentParser(description='Process GRASP / ILS - CC log files.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the log files (one for each experiment / MH)')
    parser.add_argument('--filefilter', default='*.log', required=False,
                        help='the file extension for log files (default: *.log)')
    parser.add_argument('--labels', nargs='+',
                        help='the experiment labels (one for each experiment / MH)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    labels = args.labels

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter

    processCCResult(folders, labels, filter)

def convert_jobid(jobid):

    experiment_name = 'Unknown'
    if 19408 <= jobid <= 19421 or  18978 <= jobid <= 18982: experiment_name ='ils-seq--sdot-seq-l1-ils10s-a1.0-fi-imb-2h'
    if 19829 <= jobid <= 19853: experiment_name ='ils-seq-rcc-SetParA--unga-seq-l1-ils10s-a0.4-fi-imb-1h'
    if 19743 <= jobid <= 19757 or  19422 <= jobid <= 19436 or  18983 <= jobid <= 18987: experiment_name ='ils-seq--rbig-l1-ils10s-1x8-a1.0-8-2h'
    if 19809 <= jobid <= 19828: experiment_name ='ils-seq-SetParA--rbig-l1-ils10s-1x8-a0.4-8-2h'
    if 19687 <= jobid <= 19696 or  19457 <= jobid <= 19466 or  18988 <= jobid <= 19007: experiment_name ='grasp-SetParA-rcc--unga-seq-l2-grasp400s-a0.8-fi-imb-1h'
    if 19697 <= jobid <= 19711 or  19467 <= jobid <= 19481 or  19008 <= jobid <= 19012: experiment_name ='grasp-bullmpi--sdot-seq-l1-grasp400s-a1.0-fi-imb-2h'
    if 19712 <= jobid <= 19727 or  19482 <= jobid <= 19497 or  19013 <= jobid <= 19017: experiment_name ='grasp-bullmpi--rbig-seq-l2-grasp400s-a0.8-fi-imb-1h'
    if 19728 <= jobid <= 19742 or  19498 <= jobid <= 19512 or  19018 <= jobid <= 19022: experiment_name ='grasp-bullmpi--rbig-seq-l1-grasp400s-a1.0-fi-imb-2h'
    if 19854 <= jobid <= 19868: experiment_name ='grasp-bullmpi--agents-seq-l1-grasp400s-a1.0-fi-imb-2h'
    if 19377 <= jobid <= 19387 or  18964 <= jobid <= 18968: experiment_name ='ils--rbig-l1-ils10p-2x8-a1.0-16-2h'
    if 19388 <= jobid <= 19407 or  18973 <= jobid <= 18977: experiment_name ='ils-par-puro--sdot-l1-ils10p-a1.0-16-2h'
    if 19657 <= jobid <= 19671 or  19513 <= jobid <= 19527 or  19028 <= jobid <= 19032: experiment_name ='grasp-400wi-rcc--sdot-l1-grasp8p-1x8-a1.0-8-2h'
    if 19672 <= jobid <= 19686 or  19528 <= jobid <= 19542 or  19033 <= jobid <= 19037: experiment_name ='grasp-400wi-rcc--rbig-l1-grasp8p-1x8-a1.0-8-2h'
    if 19627 <= jobid <= 19641 or  19554 <= jobid <= 19563 or 19543 <= jobid <= 19547 or  19038 <= jobid <= 19042: experiment_name ='grasp-bullmpi--sdot-l1-grasp8p-pvnd-a1.0-32-2h'
    if 19642 <= jobid <= 19656 or  19564 <= jobid <= 19573 or 19548 <= jobid <= 19552 or  19043 <= jobid <= 19047: experiment_name ='ils--sdot-l1-ils10p-a1.0-pvnd2p-32-2h'
    if 20008 or  19869 <= jobid <= 19893: experiment_name ='grasp-bullmpi--sdot-l1-grasppar-1p-pvnd-3p-a1.0-4-2h'
    if 19919 <= jobid <= 19943: experiment_name ='ils--sdot-l1-ils1p-a1.0-pvnd-3p-4-2h'
    if 19894 <= jobid <= 19918: experiment_name ='grasp-bullmpi--sdot-l1-grasppar-2p-pvnd-3p-a1.0-8-2h'
    if 19944 <= jobid <= 19967: experiment_name ='ils--sdot-l1-ils2p-a1.0-pvnd-3p-8-2h'
    if 20016 <= jobid <= 20040: experiment_name ='grasp-bullmpi--rbig-l1-grasp8p-pvnd2p-a1.0-32-2h'
    if 20041 <= jobid <= 20060: experiment_name ='ils--rbig-l1-ils10p-pvnd2p-a1.0-32-2h'
    if 20061 <= jobid <= 20079: experiment_name ='grasp-bullmpi--sdot-l1-vote-boem-a0.0-1-2h'
    return experiment_name

def processCCResult(folders, labels, filter):

    # for each folder (algorithm), determines the best solution value
    # compares the solution values of each folder and uses the worst solution value as Target I(P)
    # Using this Target I(P), determines the average time spent by each algorithm and outputs a table with the comparison
    # Also outputs a sequence of all execution times of each algorithm, to allow the creation of TTTplots

    bestSolValues = []
    worstSolValues = []

    for folder in folders:
        print "Processing folder " + ''.join(folder)
        # CC results
        construction_times_for_experiment = dict()
        search_times_for_experiment = dict()

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            #print "Processing folder " + ''.join(root)
            if (len(files)):
                file_list = []
                file_list.extend(glob.glob(root + "/" + str(filter)))
                count = len(file_list) - 1
                print 'Found ' + str(count) + ' files'

                while count >= 0:
                    # Process all files from all algorithm executions, including those in parallel
                    print "Processing file " + file_list[count] + "\n"
                    filename = file_list[count]
                    jobid = (filename[filename.rfind('/') + 1: filename.rfind('.uff0')])
                    experiment_name = convert_jobid(long(jobid))
                    print 'Jobid is ' + str(jobid) + ' and experiment name is ' + str(experiment_name)

                    if not construction_times_for_experiment.has_key(experiment_name):
                        construction_times_for_experiment[experiment_name] = dict()
                    if not search_times_for_experiment.has_key(experiment_name):
                        search_times_for_experiment[experiment_name] = dict()

                    content_file = open(file_list[count], 'r')
                    try:
                        content = content_file.read()

                        reader = csv.reader(StringIO.StringIO(content), delimiter=']')
                        for row in reader:
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if 'Reading input file' in linestring:
                                # obtains the name of the instance file being processed
                                instance_name = linestring[linestring.rfind('/') + 1:linestring.rfind('\'')]

                                if not construction_times_for_experiment[experiment_name].has_key(instance_name):
                                    construction_times_for_experiment[experiment_name][instance_name] = []
                                if not search_times_for_experiment[experiment_name].has_key(instance_name):
                                    search_times_for_experiment[experiment_name][instance_name] = []

                                print 'Processing times for file ' + str(instance_name)
                            if 'Time spent on construction phase' in linestring:
                                construction_time = float(linestring[linestring.rfind(':') + 2: linestring.rfind('s')])
                                construction_times_for_experiment[experiment_name][instance_name].append(construction_time)
                                print 'Construction time: ' + str(construction_time)
                            if 'Time spent on local search' in linestring:
                                search_time = float(linestring[linestring.rfind(':') + 2: linestring.rfind('s')])
                                search_times_for_experiment[experiment_name][instance_name].append(search_time)
                                print 'Local search time: ' + str(search_time)
                        count = count - 1
                    finally:
                        content_file.close()

                # end loop
                # process last file
        print "\nProcessing CC Results...\n"

        # Save CC results of all executions of all instances to csv file: Time to reach best solution and total time spent
        #result_file = open(folder + "/summary.csv", "w")
        #print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
        #result_file.write("ExecutionID; I(P); I(P)+; I(P)-; k; IterBestSol; TimeToBestSol; Global time(s); Params; Total Iter; NumVisitedSolutions\n")
        #for key in sorted(all_files_summary.iterkeys()):
            #print "%s, %s" % (key, all_files_summary[key])
        #    result_file.write("%s; %s\n" % (key, all_files_summary[key]))
        #result_file.close()
        print "------ CC Best results:"
        print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
        for experiment_name in sorted(construction_times_for_experiment.iterkeys()):
            print 'Experiment ' + str(experiment_name)
            print 'Instance, Construction time, Local search time'
            for instance_name in sorted(construction_times_for_experiment[experiment_name]):
                print "%s, %s, %s" % (str(instance_name), str(np.mean(construction_times_for_experiment[experiment_name][instance_name])), str(np.mean(search_times_for_experiment[experiment_name][instance_name])))
            print '\n'

    # compare the best results of each algorithm, choosing the worst
    worse_sol = dict()
    best_sol = dict()

    for worst_result_dict in worstSolValues:
        for key in sorted(worst_result_dict.iterkeys()):
            instance_name = str(key)
            sol_value = float(worst_result_dict[key])
            if worse_sol.has_key(instance_name):
                current_value = float(worse_sol[instance_name])
                if sol_value > current_value:
                    worse_sol[instance_name] = str(sol_value)
            else:
                worse_sol[instance_name] = str(sol_value)

    for best_result_dict in bestSolValues:
        for key in sorted(best_result_dict.iterkeys()):
            instance_name = str(key)
            sol_value = float(best_result_dict[key])
            if best_sol.has_key(instance_name):
                current_value = float(best_sol[instance_name])
                if sol_value < current_value:
                    best_sol[instance_name] = str(sol_value)
            else:
                best_sol[instance_name] = str(sol_value)

    print "\n\nBest solution values for each instance, comparing all algorithms:"
    InstanceTargetSolDataSet = list(zip(best_file_summary.iterkeys(), best_file_summary.itervalues()))
    df = pd.DataFrame(data = InstanceTargetSolDataSet, columns=['Instance', 'Target I(P)'])
    print df

    print "\n\nWorst solution values for each instance, comparing all algorithms:"
    InstanceTargetSolDataSet2 = list(zip(worst_file_summary.iterkeys(), worst_file_summary.itervalues()))
    df2 = pd.DataFrame(data = InstanceTargetSolDataSet2, columns=['Instance', 'Target I(P)'])
    print df2

    # for each instance in worse_sol, discover the time each algorithm took to reach that solution
    ttt_for_algorithm = dict()  # contains the ttt dict for a specific algorithm / folder
    for folder in folders:
        ttt_for_algorithm[folder] = processTimeToTarget(folder, worse_sol)  # worse_sol, best_sol

    # now, for each instance, we need to calculate the columns: AvgTimeToTarget, Stddev(TTT), Speedup.
    # and also the Standard Error of the Mean (SEM)
    AvgTimeToTarget = dict()
    StddevTTT = dict()
    Speedup = []
    for folder in folders:
        AvgTimeToTarget[folder] = []
        StddevTTT[folder] = []
    for instance_name in best_file_summary.iterkeys():
        # for each algorithm,
        count = 0
        data_to_plot = []
        for folder in folders:
            tttdict = dict(ttt_for_algorithm[folder])
            tttlist = tttdict[instance_name]
            # obtains a list of execution times (time-to-target) of the algorithm for this specific instance
            # calculates the average, stddev and confidence internal for the ttt variable
            avg_ttt =  np.mean(tttlist)
            AvgTimeToTarget[folder].append(avg_ttt)
            stddev_ttt = (np.std(tttlist, ddof=1)/np.sqrt(len(tttlist)) * 1.96) # numpy.std(tttlist)
            StddevTTT[folder].append(stddev_ttt)

            # student's t-test for time-to-target
            #(tstat, pvalue) = stats.ttest_1samp(list(tttlist), avg_ttt)
            SEM = stats.sem(list(tttlist))
            data_to_plot.append(list(tttlist))

            # print 'For %s solved by %s : StdErrorOfMean = %6.4f' % (instance_name, labels[count], SEM)
            count += 1


        # box plot of time-to-target
        # Create a figure instance
        #fig = plt.figure(1, figsize=(9, 6))

        # Create an axes instance
        #ax = fig.add_subplot(111)
        #ax.set_xticklabels(list(labels))

        # Create the boxplot
        #bp = ax.boxplot(list(data_to_plot))
        #plt.show()

        # Save the figure
        #fig.savefig(str(instance_name) + '-boxplot.png', bbox_inches='tight')

    # TODO export TTT series to be used in TTT-plot-latex
    print '\nTTT data for use with TTT-Latex\n'
    for instance_name in best_file_summary.iterkeys():
        print str(instance_name) + ':\n'
        for folder in folders:
            tttdict = dict(ttt_for_algorithm[folder])
            tttlist = tttdict[instance_name]
            print str(folder[folder.rfind('/') + 1:]) + ': ' + str(tttlist).replace(',', '').replace('[', '').replace(']', '') + '\n'


    cols = ['Instance', 'Target I(P)']  # 'AvgTimeToTarget', 'Stddev-TTT', 'Speedup'
    CompleteInstanceDataSet = {'Instance' : list(worst_file_summary.iterkeys()), 'Target I(P)' : list(worst_file_summary.itervalues())}
    count = 0
    for folder in folders:
        colname = 'AvgTimeToTarget-' + str(labels[count])
        CompleteInstanceDataSet[colname] = AvgTimeToTarget[folder]
        cols.append(colname)
        colname = 'Stddev-TTT-' + str(labels[count])
        CompleteInstanceDataSet[colname] = StddevTTT[folder]
        cols.append(colname)
        count += 1

    df = pd.DataFrame(CompleteInstanceDataSet)
    print df[cols]

    experiment_name = ''
    for folder in folders:
        experiment_name = experiment_name + folder[folder.rfind('/') + 1:] + '-vs-'

    tdir = './compare_results'
    if not os.path.exists(tdir):
        os.makedirs(tdir)
    df.to_csv(tdir + '/' + str(experiment_name) + '.csv')


def processTimeToTarget(folder, worse_sol):

    print "\nProcessing time-to-target for result folder = " + str(folder) + "\n"
    instance_dict = dict()
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        #print "Processing folder " + ''.join(root)
        if (len(files) and ''.join(root) != folder):
            file_list = []
            file_list.extend(glob.glob(root + "/CC*-iterations.csv"))
            file_list.extend(glob.glob(root + "/Node*-iterations.csv"))
            count = len(file_list) - 1

            # Process CC results
            if os.path.isfile(root + "/cc-result.txt"):
                filename = (root[:root.rfind("/")])
                datetime = root[root.rfind("/") + 1:]
                filename = filename[filename.rfind("/") + 1:]
                myfolder = root[:root.rfind("/")]
                #print "Processing instance " + myfolder + ", " + filename
                if instance_dict.has_key(filename):
                    if instance_dict[filename] <> myfolder:
                        print "WARN: problem with result folder structure!\n"
                else:
                    instance_dict[filename] = root[:root.rfind("/")]
    # end result folder scan
    tttdict = dict()  # tttdict contains the time-to-target list for each instance for this algorithm
    for key in sorted(instance_dict.iterkeys()):
        instance_name = str(instance_dict[key])
        print "Processing instance " + str(instance_name) + ", " + str(key)
        tttlist = processTimeToTargetOnInstance(instance_name, str(key), worse_sol)
        tttdict[key] = tttlist

    return tttdict

def processTimeToTargetOnInstance(folder, filename, worse_sol):
    # CC results
    all_files_summary = dict()
    full_iteration_data = dict()
    target_iteration_history = dict()

    print "\n\nTime-to-target for instance " + str(filename)
    target_value = float(worse_sol[filename])
    print "\nDesired target I(P) = " + str(target_value)

    previous_filename = ""
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        #print "Processing folder " + ''.join(root)
        if (len(files) and ''.join(root) != folder):
            file_list = []
            file_list.extend(glob.glob(root + "/CC*-iterations.csv"))
            file_list.extend(glob.glob(root + "/Node*-iterations.csv"))
            count = len(file_list) - 1

            # Process CC results
            if os.path.isfile(root + "/cc-result.txt"):
                filename = (root[:root.rfind("/")])
                datetime = root[root.rfind("/") + 1:]
                filename = filename[filename.rfind("/") + 1:]
                # print "Processing instance " + filename

                value = 0
                pos_value = 0
                neg_value = 0
                K = 0
                iteration = 0
                time = 0

                while count >= 0:
                    # Process all files from all metaheuristic executions, including those in parallel
                    # print "Processing file " + file_list[count] + "\n"
                    content_file = open(file_list[count], 'r')
                    full_iteration_data.clear()
                    try:
                        content = content_file.read()

                        reader = csv.reader(StringIO.StringIO(content), delimiter=',')
                        linecount = 0
                        for row in reader:
                            if linecount == 0:  # skip the first line of csv file
                                linecount += 1
                                continue
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if linestring.startswith('Best value'):
                                break

                            # obtains each iteration result found by a specific execution of a specific node (can be parallel)
                            iteration = long(linecount)
                            value = float(column[1])
                            pos_value = float(column[2])
                            neg_value = float(column[3])
                            K = long(column[4])
                            time = float(column[5])

                            key = file_list[count]
                            iter_data = [iteration, value, pos_value, neg_value, K, time, key]

                            if not full_iteration_data.has_key(key):
                                full_iteration_data[key] = [iter_data]
                            full_iteration_data[key].append(iter_data)
                            linecount += 1
                        count = count - 1
                    finally:
                        content_file.close()

                    # for each execution of the algorithm
                    for executionId in sorted(full_iteration_data.iterkeys()):
                        executionIdIterations = full_iteration_data[executionId]
                        desired_iteration = -1
                        # for each iteration of an execution of the algorithm
                        for iter_data in executionIdIterations:
                            # find out in which iteration the algorithm reaches the target solution value
                            sol_value = float(iter_data[1])
                            if sol_value <= target_value:
                                desired_iteration = iter_data
                                #print "Found target at " + iter_data[5] + " seconds.\n"
                                break
                        # stores the info about the target iteration when the target solution was found
                        if desired_iteration < 0:
                            # if the code gets to this point, it means the algorithm was running in parallel and one
                            # of the nodes did not get to the best solution => can be safely ignored
                            continue
                        jobId = executionId[:executionId.rfind('/')]
                        jobId = jobId[jobId.rfind('/') + 1:]
                        if target_iteration_history.has_key(jobId):
                            current_iteration_info = target_iteration_history[jobId]
                            current_iteration_time = float(current_iteration_info[5])
                            new_iteration_time = float(desired_iteration[5])
                            if current_iteration_time > new_iteration_time:
                                # replaces current iteration info with another one, which is faster
                                # only happens in parallel algorithms
                                target_iteration_history[jobId] = desired_iteration
                                # print "found better result in the same job: " + str(jobId) + ", time = " + str(new_iteration_time) +  "\n"
                        else:
                            target_iteration_history[jobId] = desired_iteration
                # end process all result files for one execution of the algorithm
        # end process all executions of the algorithm

    # now we know the time-to-target of each algorithm execution for this instance
    print "Total: " + str(len(target_iteration_history)) + " results"
    print "Time, Iter, Info"
    tttlist = []
    for key in target_iteration_history.iterkeys():
        iteration_info = target_iteration_history[key]
        tttlist.append(iteration_info[5])
        print "%s %s %s" % (str(iteration_info[5]), str(long(iteration_info[0])), str(iteration_info))

    print "\nFinished TTT processing for instance " + str(filename) + ".\n"
    return tttlist


if __name__ == "__main__":
    main(sys.argv[1:])

