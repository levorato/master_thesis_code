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

    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = [f for f in files if filter in f]

                for file in file_list:
                    directory = os.path.join(root, 'dm_results')
                    if not os.path.exists(directory):
                        os.makedirs(directory)

                    # measure elapsed time during instance generation
                    start = time.time()

                    file = os.path.join(root, file)
                    filename = file
                    print "\nProcessing file " + filename

                    try:
                        invoke_pajek_solve_dm(filename, directory)
                    except Exception,e:
                        print "Error invoking Pajek for instance file {0}".format(filename)
                        print str(e)
                        traceback.print_exc()
                        continue

                    print "Invoked Pajek for instance file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Execution took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file

        print "\nDone.\n"


def invoke_pajek_solve_dm(instance_file, output_folder):

    pajek_path = 'C:\Users\czt0\Downloads\Pajek4'
    pajek_exe_path = pajek_path + '\Pajek.exe'
    # open the instance_file and retrieves the number of vertices n
    n = 0
    with open(instance_file, "r") as i_file:
        content = i_file.read()
        reader = csv.reader(StringIO.StringIO(content), delimiter=',')
        for row in reader:
            linestring = ''.join(row)
            print linestring
            linestring = linestring[linestring.find(' ') + 1:]
            n = long(linestring)
            break
    print "Instance with n = " + str(n)

    # Determine the number of clusters for each file, depending on n
    k_dict = {200 : 7, 300 : 9, 400 : 6, 600 : 9, 800 : 14, 1000 : 18, 2000 : 31}

    # Single argument: file with gantt data
    outputfile = os.path.join(output_folder, instance_file + ".rep")
    solutionfile = os.path.join(output_folder, instance_file + ".txt")

    if n in k_dict.keys():
        with open(os.path.join(output_folder, "Pajek.log"), "w") as t_file:
            try:
                t_file.write('NETBEGIN 1\n')
                t_file.write('CLUBEGIN 1\n')
                t_file.write('PERBEGIN 1\n')
                t_file.write('CLSBEGIN 1\n')
                t_file.write('HIEBEGIN 1\n')
                t_file.write('VECBEGIN 1\n')
                t_file.write('\n')
                t_file.write('Msg Reading Network   ---    ' + instance_file + '\n')
                t_file.write('N 1 RDN "' + instance_file + '" (' + str(n) + ')\n')
                t_file.write('C 1 RANDOMC ' + str(k_dict[n]) + ' (' + str(n) + ')\n')
                t_file.write('Msg Partitioning Signed Networks according to Structural Balance\n')
                t_file.write('C 2 BALANCE 1 1 [1000 0.500 0 1 0 0] (' + str(n) + ')\n')
                t_file.write('SAVEREPORT "' + outputfile + '"\n')
                t_file.write('Msg Saving partition to file   ---    ' + solutionfile + '\n')
                t_file.write('C 2 WC "' + solutionfile + '" (' + str(n) + ')\n')
                t_file.write('EXIT\n')
            except Exception,e:
                print 'Processing error, skipping file.'

        # runs pajek using the macro file create before
        proc = subprocess.Popen(pajek_exe_path, cwd=output_folder, shell=True, stdin=subprocess.PIPE, )
        proc.communicate()  # wait for process to finish
        #proc.stdin.write('\nquit\n') #close the gnuplot window
    else:
        print 'Skipping file, since k value was not found.'

if __name__ == "__main__":
    main(sys.argv[1:])



