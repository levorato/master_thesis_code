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
    parser = argparse.ArgumentParser(description='Invoke Pajek via command-line to solve .net files via Louvain Method.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the graph files in Pajek (.net) format')
    parser.add_argument('--filefilter', default='.net', required=False,
                        help='the filename extension for graph files (default: .net)')
    parser.add_argument('--pajekpath', default='C:\Users\czt0\Downloads\Pajek4', required=False,
                        help='the installation folder of Pajek')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    pajek_path = args.pajekpath

    print 'Graph file folders are ', folders
    print 'File filter is ', filter
    print 'Pajek installation folder is ', pajek_path

    processInstanceFiles(folders, filter, pajek_path)


def processInstanceFiles(folders, filter, pajek_path):

    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = [f for f in files if filter in f]

                for exec_num in xrange(0, 25):
                    for file in file_list:
                        directory = os.path.join(root, 'louvain_results')
                        if not os.path.exists(directory):
                            os.makedirs(directory)
                        print "Output folder is " + str(directory)

                        file = os.path.join(root, file)
                        filename = file
                        print "\nProcessing file " + filename

                        # measure elapsed time during instance generation
                        start = time.time()
                        try:
                            invoke_pajek_solve_louvain(pajek_path, filename, directory, exec_num+1, True)
                        except Exception,e:
                            print "Error invoking Pajek Louvain single-refinement for instance file {0}".format(filename)
                            print str(e)
                            traceback.print_exc()
                            #continue

                        print "Invoked Pajek Louvain single-refinement for instance file {0}".format(filename)
                        end = time.time()
                        elapsed = end - start
                        print "Execution took {0:.2f} seconds.".format(elapsed)

                        # measure elapsed time during instance generation
                        start = time.time()
                        try:
                            invoke_pajek_solve_louvain(pajek_path, filename, directory, exec_num+1, False)
                        except Exception, e:
                            print "Error invoking Pajek Louvain multiple-refinement for instance file {0}".format(filename)
                            print str(e)
                            traceback.print_exc()
                            #continue

                        print "Invoked Pajek Louvain multiple-refinement for instance file {0}".format(filename)
                        end = time.time()
                        elapsed = end - start
                        print "Execution took {0:.2f} seconds.".format(elapsed)
                    # end loop
                    # process last file

        print "\nDone.\n"


def invoke_pajek_solve_louvain(pajek_path, instance_file, output_folder, exec_num, single_refinement = True):

    print 'Execution number ' + str(exec_num)
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

    if single_refinement:
        suffix = 'single'
    else:
        suffix = 'ML'
    base_filename = instance_file[instance_file.rfind(os.path.sep)+1:instance_file.rfind('.')]
    outputfile = os.path.join(output_folder, base_filename + '-' + suffix + '-' + str(exec_num) + ".rep")
    solutionfile = os.path.join(output_folder, base_filename + '-' + suffix + '-' + str(exec_num) + ".txt")

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
            if single_refinement:
                t_file.write('Louvain Community Detection for Signed Network. Coarsening only.\n')
                t_file.write('C 1 MODULARITY 1 10000 10000 10000 1.0000000000 0.500000 (' + str(n) + ')\n')
            else:
                t_file.write('Louvain Community Detection for Signed Network. Coarsening and Refinement. \n')
                t_file.write('C 1 MODULARITYML 1 10000 10000 10000 10000 1.0000000000 0.500000 (' + str(n) + ')\n')
            t_file.write('SAVEREPORT "' + outputfile + '"\n')
            t_file.write('Msg Saving partition to file   ---    ' + solutionfile + '\n')
            t_file.write('C 1 WC "' + solutionfile + '" (' + str(n) + ')\n')
            t_file.write('EXIT\n')
        except Exception,e:
            print 'Processing error, skipping file.'

    # runs pajek using the macro file create before
    proc = subprocess.Popen(pajek_exe_path, cwd=output_folder, shell=True, stdin=subprocess.PIPE, )
    proc.communicate()  # wait for process to finish
    #proc.stdin.write('\nquit\n') #close the gnuplot window

if __name__ == "__main__":
    main(sys.argv[1:])



