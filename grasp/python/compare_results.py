import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import argparse
#from pandas import *


def main(argv):

    csv.field_size_limit(sys.maxsize)

    parser = argparse.ArgumentParser(description='Compare GRASP / ILS - CC result files.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the result files (one for each experiment / MH)')
    parser.add_argument('--filefilter', default='*.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter

    processCCResult(folders)


def processCCResult(folders):

    # for each folder (algorithm), determines the best solution value
    # compares the solution values of each folder and uses the worst solution value as Target I(P)
    # Using this Target I(P), determines the average time spent by each algorithm and outputs a table with the comparison
    # Also outputs a sequence of all execution times of each algorithm, to allow the creation of TTTplots

    bestSolValues = []

    for folder in folders:
        print "Processing folder " + ''.join(folder)
        # CC results
        all_files_summary = dict()
        best_file_summary = dict()
        previous_filename = ""
        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            #print "Processing folder " + ''.join(root)
            if (len(files) and ''.join(root) != folder):
                file_list = []
                file_list.extend(glob.glob(root + "/CC*.csv"))
                file_list.extend(glob.glob(root + "/Node*.csv"))
                count = len(file_list) - 1

                # Process CC results
                if os.path.isfile(root + "/cc-result.txt"):

                    filename = (root[:root.rfind("/")])
                    filename = filename[filename.rfind("/")+1:]
                    datetime = root[root.rfind("/") + 1:]

                    best_value = 100000000L
                    value = 0

                    while count >= 0:
                        # Process all files from all metaheuristic executions, including those in parallel
                        #print "Processing file " + file_list[count] + "\n"
                        content_file = open(file_list[count], 'r')
                        try:
                            content = content_file.read()

                            reader = csv.reader(StringIO.StringIO(content), delimiter=',')
                            for row in reader:
                                linestring = ''.join(row)
                                column = []
                                for col in row:
                                    column.append(col)
                                if linestring.startswith('Best value'):
                                    # obtains the best result found by a specific execution of a specific node (can be parallel)
                                    value = float(column[1])
                                    if value < best_value:
                                        # if the results from the execution of this node are better, replace the data
                                        best_value = value
                            count = count - 1
                        finally:
                            content_file.close()

                    # captura o melhor resultado dadas todas as execucoes de um mesmo grafo
                    if best_file_summary.has_key(filename):
                        current_value = float(best_file_summary[filename])
                        if (best_value < current_value):
                            best_file_summary[filename] = str(best_value)
                    else:
                        best_file_summary[filename] = str(best_value)

                    previous_filename = filename
                # end loop
                # process last file
        print "\nProcessing CC Results...\n"
        bestSolValues.append(best_file_summary)


        # Group the results by instance (filename)
        for key in sorted(all_files_summary.iterkeys()):
            with open(folder + '/' + str(key) + '-all_executions.csv', 'wb') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)

                csvwriter.writerow(['','ExecutionID','BestSol','I(P)+','I(P)-','k','IterBestSol','TimeToBestSol','Total time','Params','Total iter','NumVisitedSolutions'])
                executions = all_files_summary[key]
                for item in executions:
                    csvwriter.writerow(str(item))



        # TODO export (to csv) the results of the time spent to get to a specific target I(P) for each instance



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
        for key in sorted(best_file_summary.iterkeys()):
            print "%s, %s" % (key, best_file_summary[key])

    # compare the best results of each algorithm, choosing the worst
    worse_sol = dict()

    for best_result_dict in bestSolValues:
        for key in sorted(best_result_dict.iterkeys()):
            instance_name = str(key)
            sol_value = float(best_result_dict[key])
            if worse_sol.has_key(instance_name):
                current_value = float(worse_sol[instance_name])
                if sol_value > current_value:
                    worse_sol[instance_name] = str(sol_value)
            else:
                worse_sol[instance_name] = str(sol_value)

    print "\n\nWorst solution values for each instance"
    for key in sorted(best_file_summary.iterkeys()):
        print "%s, %s" % (key, best_file_summary[key])

    # for each instance in worse_sol, discover the time each algorithm took to reach that solution
    for folder in folders:
        processTimeToTarget(folder, worse_sol)

def processTimeToTarget(folder, worse_sol):
    # CC results
    all_files_summary = dict()
    full_iteration_data = dict()
    target_iteration_history = dict()

    previous_filename = ""
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        print "Processing folder " + ''.join(root)
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

                value = 0
                pos_value = 0
                neg_value = 0
                K = 0
                iteration = 0
                time = 0

                while count >= 0:
                    # Process all files from all metaheuristic executions, including those in parallel
                    print "Processing file " + file_list[count] + "\n"
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
                            iteration = linecount
                            value = float(column[1])
                            pos_value = float(column[2])
                            neg_value = float(column[3])
                            K = long(column[4])
                            time = float(column[5])

                            key = file_list[count]
                            iter_data = [iter, value, pos_value, neg_value, K, time, key]

                            if not full_iteration_data.has_key(key):
                                full_iteration_data[key] = [iter_data]
                            full_iteration_data[key].append(iter_data)
                            linecount += 1
                        count = count - 1
                    finally:
                        content_file.close()

                    target_value = float(worse_sol[filename])
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
                                break
                        # stores the info about the target iteration when the target solution was found
                        if desired_iteration < 0:
                            print "Error retrieving target iteration info for executionId = " + str(executionId)
                        if target_iteration_history.has_key(filename):
                            target_iteration_history[filename].append(desired_iteration)
                        else:
                            target_iteration_history[filename] = [desired_iteration]

                # now we know the time-to-target of each algorithm execution for each instance file
                print "\n\nDesired target I(P) = " + str(target_value)
                print "\nTime-to-target for instance " + str(filename)
                print "Time, Iter, Info"
                for iteration_info in target_iteration_history[filename]:
                    print "%s %s %s" % (str(iteration_info[5]), str(iteration_info[0]), str(iteration_info))

            # end loop
            # process last file
    print "\nProcessing CC Results...\n"


if __name__ == "__main__":
    main(sys.argv[1:])

