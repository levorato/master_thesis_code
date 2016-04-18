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

def generate_gnuplot_gantt_from_solution(solution_file, output_folder, alphasort = False,
                                         colorfile = False, plottitle = "Default"):

    # Single argument: file with gantt data
    ganttfile = solution_file
    base_folder = ganttfile[:ganttfile.rfind(os.sep)]
    filename = base_folder[:base_folder.rfind(os.sep)]
    filename = filename[filename.rfind(os.sep)+1:]
    gnufilename = filename + '.gpl'
    outputfile = os.path.join(output_folder, gnufilename)

    # obtain the objective function value (z) and makespan (Termino) from the result.txt file
    z, Cmax, time_spent = extract_solution_value(os.path.join(base_folder, "result.txt"))

    options = gantt1.make_default_options()
    options.alphasort = alphasort
    options.plottitle = filename + " ; z = {:.2f} (Cmax = {:d}); t = {:.2f} s".format(z, Cmax, time_spent)
    file_content = gantt1.compute(options, ganttfile)

    options.outputfile = outputfile
    gantt1.compute(options, ganttfile)

    # runs gnuplot to generate a png file with the gantt chart
    proc = subprocess.Popen(['gnuplot','-p'],
                        shell=True,
                        stdin=subprocess.PIPE,
                        )
    proc.stdin.write('set terminal png\n')
    proc.stdin.write("set output '" + outputfile + ".png'\n")
    proc.stdin.write(file_content)
    proc.stdin.write('\nquit\n') #close the gnuplot window


def extract_solution_value(result_filepath):

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

    try:
        content_file2 = open(result_filepath, 'r')
        text_content = content_file2.read()
        output = StringIO.StringIO(text_content)
        reader = csv.reader(output, delimiter =',')

        row_count = 1
        Cmax = -1
        Termino = -1
        for row in reader:
            linestring = ''.join(row)
            column = []
            for col in row:
                column.append(col)
            if len(row) == 0:
                continue

            line = column[0]
            if row_count == 1:  # linha 1

                z = float(line[line.find('objective') + 10:])
            else:
                break
            row_count += 1
        # end read lines
        content_file2.close()
        print "Successfully read result file."
        return z, Cmax, time_spent
    except Exception,e:
        print "Error obtaining result file {0}".format(result_filepath)
        print str(e)
        traceback.print_exc()
        return -1, -1, -1


def main(argv):
    parser = argparse.ArgumentParser(description='Invoke Pajek via command-line to solve .net files.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the graph files in Pajek (.net) format')
    parser.add_argument('--filefilter', default='.net', required=False,
                        help='the filename extension for graph files (default: .net)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter

    args = parser.parse_args()

    print 'Graph file folders are ', folders
    print 'File filter is ', filter

    processInstanceFiles(folders, filter)


def processInstanceFiles(folders, filter):

# Template for Pajek macro file (commandline invocation):
# NETBEGIN 1
# CLUBEGIN 1
# PERBEGIN 1
# CLSBEGIN 1
# HIEBEGIN 1
# VECBEGIN 1
#
# Msg Reading Network   ---    C:\Users\czt0\Documents\Mestrado\git\grasp\tests\Instances\Social Media\slashdot-undirected\pajek\slashdot-undirected-size200-part0.g.net
# N 1 RDN "C:\Users\czt0\Documents\Mestrado\git\grasp\tests\Instances\Social Media\slashdot-undirected\pajek\slashdot-undirected-size200-part0.g.net" (200)
# C 1 RANDOMC 7 (200)
# Msg Partitioning Signed Networks according to Structural Balance
# C 134 BALANCE 1 1 [1000 0.500 0 1 0 0] (200)
# Msg Saving Partition(s)/Vector(s) to tab delimited file   ---    C:\Users\czt0\Documents\Mestrado\git\grasp\tests\Instances\Social Media\slashdot-undirected\pajek\solution_1.txt
# C 134 WTABCLUVEC "C:\Users\czt0\Documents\Mestrado\git\grasp\tests\Instances\Social Media\slashdot-undirected\pajek\solution_1.txt" 0 1 2 11 1
# Msg Partitioning Signed Networks according to Structural Balance
# C 269 BALANCE 1 134 [1000 0.500 0 1 0 0] (200)
# SAVEREPORT "C:\Users\czt0\Documents\Mestrado\git\grasp\tests\Instances\Social Media\slashdot-undirected\pajek\sdot-200.rep"


    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = [f for f in files if filter in f]

                for file in file_list:
                    directory = os.path.join(root, 'gantt')
                    if not os.path.exists(directory):
                        os.makedirs(directory)

                    # measure elapsed time during instance generation
                    start = time.time()

                    file = os.path.join(root, file)
                    filename = file
                    print "\nProcessing file " + filename

                    try:
                        generate_gnuplot_gantt_from_solution(filename, directory)
                    except Exception,e:
                        print "Error creating gantt chart for solution file {0}".format(filename)
                        print str(e)
                        traceback.print_exc()
                        continue

                    print "Created gantt chart for solution file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Generation took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file

        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])


# (200, 7), (300, 9), (400, 6), (600, 9), (800, 14), (1000, 18), (2000, 31)
