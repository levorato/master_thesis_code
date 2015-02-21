# Script de processamento de resultados do RCC que aponta potenciais grupos de mediadores
# The existence of a group of individuals who have only positive relationships with everyone 
# in the network counts as imbalance in the CC Problem. Nonetheless, the individuals in this group 
# could be construed as mediators (i.e. their relations probably won't change over time) and 
# their relations should not be considered as a contribution to the imbalance of the network.

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re

import numpy as np
import multiprocessing as mp
from scipy.sparse import dok_matrix
from multiprocessing import Pool

# Define an output queue
output = mp.Queue()

#gambit para python 2.4
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
    print [ tryint(c) for c in re.split('([0-9]+)', s) ][1]
    return int([ tryint(c) for c in re.split('([0-9]+)', s) ][1])

def natsorted(l):
    """ Sort the given list in the way that humans expect.
    """
    temp = [(alphanum_key(x), x) for x in l]
    temp.sort()

def process_edges(matrix, x, n, output):
          edges = ""
          size = n / double(4)
          start = x * size
          end = (x + 1) * size - 1
          print "Processing from " + str(start) + " to " + str(end)
          for i in xrange(start, end):
              print "i = " + str(i)
              for j in xrange(N):
                  if (i < j):
                      mij = matrix[i,j]
                      mji = matrix[j,i]
                      mult = mij * mji
                      if (mult > 0):  # opposite edges with the same sign become one undirected edge
                          edges.append( i + '\t' + j + '\t' + str(mij + mji) + '\n' )
                      elif (mult < 0):  # opposite edges with different signs remain as 2 diff. edges
                          edges.append( i + '\t' + j + '\t' + str(mij) + '\n' )
                          edges.append( j + '\t' + i + '\t' + str(mji) + '\n' )
                      else:  # mult = 0
                          edges.append( i + '\t' + j + '\t' + str(mij + mji) + '\n' )

          output.put(edges)

def main(argv):

   folder = ''
   filter = '*'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["filter=","instancespath="])
   except getopt.GetoptError:
      print 'convert_to_symmetric.py --filter <filter> --instancespath <instances_path>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convert_to_symmetric.py --filter <filter> --instancespath <instances_path>'
         sys.exit()
      if opt in ("-f", "--filter"):
         filter = arg
      if opt in ("-p", "--instancespath"):
         instances_path = arg
   if(instances_path == ''):
      print 'Please specify the graph instances path'
      sys.exit()
   print 'Graph instances dir is ', instances_path

   # RCC results
   all_files_summary = dict()
   print "\nProcessing Graphs...\n"

   for root, subFolders, files in os.walk(instances_path):
      # sort dirs and files
      subFolders.sort()
      files.sort()      
      print 'Found {0} files\n'.format(str(len(files)))
      print "Processing folder " + ''.join(root)
      if(len(files)):
         file_list = []
         file_list.extend(glob.glob(root + "/" + filter))
         count = len(file_list) - 1
         print 'Found {0} files\n'.format(str(count))
         
         while count >= 0:
          print "Processing file " + file_list[count] + "\n"
          prefix = file_list[count]
          graphfile = prefix
          # reads graph file, filling the adjecency matrix
          print 'Reading graph file: {0}'.format(str(graphfile))
          with open(graphfile, 'rb') as csvfile:
             lines = [line.strip() for line in open(graphfile, 'rb')]
             if str(lines[0]).find("people") >= 0:  # xpress mrel format
                 print "Detected xpress mrel format"
                 line = str(lines[0])
                 N = int(line[line.find(':')+1:].strip())
                 print "n = {0}".format(str(N))
                 # uses scipy's dictionary of keys sparse matrix (http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.dok_matrix.html#scipy.sparse.dok_matrix)
                 matrix = dok_matrix((N, N), dtype=np.float32)   # original directed graph from file
                 i = 0
                 while str(lines[i]).find("Mrel") < 0:
                    i += 1
                 line = str(lines[i])
                 line = line[line.find('[')+1:].strip()
                 i = int(line[line.find('(')+1:line.find(',')]) - 1
                 j = int(line[line.find(',')+1:line.find(')')]) - 1
                 w = float(line[line.find(')')+1:])
                 matrix[i,j] = w
                 #print "{0} {1} {2}".format(i, j, w)
                 for line in lines:
                    line = str(line)
                    if len(line) == 0:
                        continue
                    if not line[0] == '(':
                        continue
                    #print line
                    i = int(line[line.find('(')+1:line.find(',')]) - 1
                    j = int(line[line.find(',')+1:line.find(')')]) - 1
                    w = float(line[line.find(')')+1:])
                    matrix[i,j] = w
                    i += 1
             else:  # traditional graph edge tabulated format
                 dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=" \t")
                 csvfile.seek(0)
                 r = csv.reader(csvfile, dialect)
                 line1=r.next()
                 print "Detected traditional graph edge tabulated format"
                 values = [col for col in line1]
                 if(values[0].find(' ') > 0):
                    N = int(values[0][:values[0].find(' ')])
                 else:
                    N = int(values[0])
                 matrix = dok_matrix((N, N), dtype=np.float32)   # original directed graph from file
                 for line in r:
                    i = int(line[0])
                    j = int(line[1])
                    w = float(line[2])
                    matrix[i,j] = w

          print "Successfully read graph file. Converting to symmetric...\n"
          # converts the graph from directed to undirected version
          # Setup a list of processes that we want to run
          processes = [mp.Process(target=process_edges, args=(matrix, x, N, output)) for x in range(4)]
          # results = [pool.apply(process_edges, args=(matrix,x,output)) for x in range(1,N)]
          # Get process results from the output queue
          results = [output.get() for p in processes]
          
          print "Graph converted. Generating output file...\n" + str(len(results))
          # writes all edges to output graph file
          output_file = open(graphfile[:graphfile.rfind('.')] + "-sym.g", 'w')
          try:
             content = output_file.write('Test')
          finally:
              content_file.close()
                
          count = count - 1
	 # end process results


   
if __name__ == "__main__":
   main(sys.argv[1:])
