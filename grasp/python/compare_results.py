import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import argparse
from pandas import *


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

    for folder in folders:
        # CC results
        all_files_summary = dict()
        best_file_summary = dict()
        avg_file_summary = dict()
        timeInterval = dict()
        timeCount = dict()
        previous_filename = ""
        avg_ip_const = 0
        avg_value = 0
        avg_k = 0
        avg_time = 0
        avg_count = 0
        avg_iter = 0
        avg_comb = 0
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

                    text_file = open(root + "/summary.txt", "w")
                    filename = (root[:root.rfind("/")])
                    datetime = root[root.rfind("/") + 1:]
                    filename = filename[filename.rfind("/") + 1:]
                    text_file.write("CC Summary for graph file: %s\n" % filename)
                    local_avg_ip_const = 0
                    local_avg_count = 0

                    best_value = 100000000L
                    best_pos_value = 0
                    best_neg_value = 0
                    best_K = 0
                    best_iteration = 0
                    best_time = 0
                    best_param = ''
                    value = 0
                    pos_value = 0
                    neg_value = 0
                    K = 0
                    iteration = 0
                    time = 0
                    total_iter = 0
                    total_comb = 0

                    while count >= 0:
                        # Process all files from all metaheuristic executions, including those in parallel
                        print "Processing file " + file_list[count] + "\n"
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
                                    filepath = ''.join(file_list[count])
                                    text_file.write(filepath[filepath.rfind("/") + 1:] + ' ' + linestring + '\n')
                                    # obtains the best result found by a specific execution of a specific node (can be parallel)
                                    value = float(column[1])
                                    pos_value = float(column[2])
                                    neg_value = float(column[3])
                                    K = long(column[4])
                                    iteration = long(column[5])
                                    time = float(column[6])
                                    total_iter = long(column[7])
                                    # totalizes the number of visited solutions of all nodes running in parallel
                                    total_comb += long(column[8])
                                    if value < best_value:
                                        # if the results from the execution of this node are better, replace the data
                                        best_value = value
                                        best_pos_value = pos_value
                                        best_neg_value = neg_value
                                        best_K = K
                                        best_iteration = iteration
                                        best_time = time
                                        best_param = filepath[filepath.rfind("/") + 1:]
                                    elif value == best_value and iteration < best_iteration:
                                        best_K = K
                                        best_iteration = iteration
                                        best_time = time
                                        best_param = filepath[filepath.rfind("/") + 1:]
                                        best_pos_value = pos_value
                                        best_neg_value = neg_value
                                elif linestring.startswith('Average initial I(P)'):
                                    # computa o valor medio da solucao inicial da fase construtiva para determinado no da execucao paralela
                                    ip_const = float(column[1])
                                    local_avg_ip_const = local_avg_ip_const + ip_const
                                    local_avg_count = local_avg_count + 1
                            count = count - 1
                        finally:
                            content_file.close()
                    text_file.close()
                    all_files_summary[filename + "/" + datetime] = str(best_value) + "; " + str(pos_value) + "; " + str(
                            neg_value) + "; " + str(best_K) + "; " + str(iteration) + "; " + str(best_time) + "; " + str(
                            global_time) + "; " + best_param + "; " + str(total_iter) + "; " + str(total_comb)

                    # calcula a media dos valores de I(P) da fase de construcao para este set de arquivos de resultado
                    if local_avg_count == 0:
                        local_avg_ip_const = 0
                    else:
                        local_avg_ip_const = local_avg_ip_const / local_avg_count
                    # armazena os valores de todas as execucoes de um mesmo grafo para calculo da media
                    if filename == previous_filename:
                        avg_comb = avg_comb + total_comb
                        avg_ip_const = avg_ip_const + local_avg_ip_const
                        avg_value = avg_value + best_value
                        avg_k = avg_k + best_K
                        avg_time = avg_time + global_time
                        avg_iter = avg_iter + total_iter
                        avg_count = avg_count + 1
                    else:
                        if avg_count > 0:
                            print "storing " + previous_filename + ", number of executions = " + str(avg_count)
                            avg_file_summary[previous_filename] = str(avg_ip_const / avg_count) + ", " + str(
                                    avg_value / avg_count) + ", " + str(avg_k / avg_count) + ", " + str(
                                    avg_time / avg_count) + ", " + str(avg_iter / avg_count) + ", " + str(
                                    avg_comb / avg_count) + ", " + str(avg_count)
                            print "average execution times for file " + previous_filename
                            tdir = "./times"
                            if not os.path.exists(tdir):
                                os.makedirs(tdir)
                            times_file = open(tdir + "/" + previous_filename + "-executionTimes.txt", "w")
                            for key, value in sorted(timeInterval.items()):
                                times_file.write(str(key) + "," + str(value / timeCount[key]) + "\n")
                            times_file.close()
                            timeInterval = dict()
                            timeCount = dict()
                        avg_comb = total_comb
                        avg_ip_const = local_avg_ip_const
                        avg_value = best_value
                        avg_k = best_K
                        avg_time = global_time
                        avg_iter = total_iter
                        avg_count = 1

                    # captura o melhor resultado dadas todas as execucoes de um mesmo grafo
                    if best_file_summary.has_key(filename):
                        element = best_file_summary[filename]
                        value = float(element[0:element.find(';') - 1])
                        if (best_value < value):
                            best_file_summary[filename] = str(all_files_summary[filename + "/" + datetime])
                    else:
                        best_file_summary[filename] = str(all_files_summary[filename + "/" + datetime])

                    previous_filename = filename
                # end loop
                # process last file
                print "storing " + previous_filename + ", number of executions = " + str(avg_count)
                avg_file_summary[previous_filename] = str(avg_ip_const / avg_count) + ", " + str(
                        avg_value / avg_count) + ", " + str(avg_k / avg_count) + ", " + str(
                    avg_time / avg_count) + ", " + str(
                        avg_iter / avg_count) + ", " + str(avg_comb / avg_count) + ", " + str(avg_count)
                print "average execution times for file " + previous_filename
                tdir = "./times"
                if not os.path.exists(tdir):
                    os.makedirs(tdir)
                    times_file = open(tdir + "/" + previous_filename + "-executionTimes.txt", "w")
                    for key, value in sorted(timeInterval.items()):
                        times_file.write(str(key) + "," + str(value / timeCount[key]) + "\n")
                    times_file.close()
                    timeInterval = dict()
                    timeCount = dict()
                    # end process CC results
        print "\nProcessing CC Results...\n"


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
        result_file = open(folder + "/summary.csv", "w")
        #print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
        result_file.write(
                "ExecutionID; I(P); I(P)+; I(P)-; k; IterBestSol; TimeToBestSol; Global time(s); Params; Total Iter; NumVisitedSolutions\n")
        for key in sorted(all_files_summary.iterkeys()):
            #print "%s, %s" % (key, all_files_summary[key])
            result_file.write("%s; %s\n" % (key, all_files_summary[key]))
        result_file.close()
        print "------ CC Best results:"
        print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
        for key in sorted(best_file_summary.iterkeys()):
            print "%s, %s" % (key, best_file_summary[key])
        print "------ CC Average results:"
        print "Instance, Avg I(P) const, Avg I(P), Avg K, Avg Time(s), Avg Iter, Avg combinations, Num executions"
        for key in sorted(avg_file_summary.iterkeys()):
            print "%s, %s" % (key, avg_file_summary[key])


if __name__ == "__main__":
    main(sys.argv[1:])

