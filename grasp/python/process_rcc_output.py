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
from scipy.sparse import dok_matrix

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

def main(argv):

   folder = ''
   filter = ''

   threshold = float(90.0)  # percentual threshold for the quantity of relashionships
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["folder=","filter=","instancespath="])
   except getopt.GetoptError:
      print 'process_rcc_output.py --folder <folder> --filter <filter> --instancespath <instances_path>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'process_rcc_output.py --folder <folder> --filter <filter> --instancespath <instances_path>'
         sys.exit()
      elif opt in ("-i", "--folder"):
         folder = arg
      if opt in ("-f", "--filter"):
         filter = arg
      if opt in ("-p", "--instancespath"):
         instances_path = arg
   if(folder == ''):
      print 'Please specify the results dir'
      sys.exit()
   if(instances_path == ''):
      print 'Please specify the graph instances path'
      sys.exit()
   print 'Input dir is ', folder
   print 'File filter is ', filter
   print 'Graph instances dir is ', instances_path
   print 'Mediator detection threshold is ', threshold

   # RCC results
   all_files_summary = dict()
   print "\nProcessing RCC Results...\n"

   for root, subFolders, files in os.walk(folder):
      # sort dirs and files
      subFolders.sort()
      files.sort()      
      print 'Found {0} files\n'.format(str(len(files)))
      print "Processing folder " + ''.join(root)
      if(len(files)):
         file_list = []
         file_list.extend(glob.glob(root + "/rcc-result.txt"))
         count = len(file_list) - 1
         print 'Found {0} files\n'.format(str(count))
         
         while count >= 0:
          print "Processing file " + file_list[count] + "\n"
          prefix = file_list[count]
          graphfile = prefix[0:prefix.rfind('/')]
          graphfile = graphfile[0:graphfile.rfind('/')]
          graphfile = graphfile[graphfile.rfind('/')+1:]
          # reads graph file, filling the adjecency matrix
          print 'Reading graph file: {0}'.format(str(instances_path + '/' + graphfile))
          with open(instances_path + '/' + graphfile, 'rb') as csvfile:
             lines = [line.strip() for line in open(instances_path + '/' + graphfile, 'rb')]
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
                 matrix = matrix = dok_matrix((N, N), dtype=np.float32)   # original directed graph from file
                 for line in r:
                    i = int(line[0])
                    j = int(line[1])
                    w = float(line[2])
                    matrix[i,j] = w

          print "Successfully read graph file."
          smatrix = matrix
          '''
          print "Converting to symmetric...\n"
          # converts the graph from directed to undirected version
          smatrix = dok_matrix((N, N), dtype=np.float32)  # symmetric version of the directed graph
          for i in xrange(N):
              print "i = " + str(i)
              for j in xrange(N):
                  if (i < j):
                      mij = matrix[i,j]
                      mji = matrix[j,i]
                      mult = mij * mji
                      if (mult > 0):  # opposite edges with the same sign become one undirected edge
                          smatrix[i,j] = mij + mji
                      elif (mult < 0):  # opposite edges with different signs remain as 2 diff. edges
                          smatrix[i,j] = mij
                          smatrix[j,i] = mji
                      else:  # mult = 0
                          smatrix[i,j] = mij + mji
          
          print "Graph converted."
          '''

          print "Reading SRCC results files..."
          # reads rcc results file, with cluster X node configuration
          content_file = open(file_list[count], 'r')
          cluster = [int(0) for x in xrange(N)]
          partitionNumber = 0
          try:
             content = content_file.read().splitlines()
             for row in content:
               linestring = ''.join(row).strip()
               # Captura a lista de particoes / clusters da solucao do RCC
               if linestring.startswith( 'Partition' ):          
                  strcluster = linestring[linestring.find("[")+1:linestring.rfind("]")]
                  # print "linha: " + strcluster.strip()
                  nodeList = strcluster.split()
                  for node in nodeList:
                    cluster[int(node.strip())] = partitionNumber
                  partitionNumber += 1
          finally:
              content_file.close()

          print "Results file OK. Processing percentage of edges..."
          # lists of interesting clusters (with mediation properties)
          PlusMediators = []
          PlusMutuallyHostileMediators = []
          InternalSubgroupHostility = []
          DifferentialPopularity = []
          # for each partition, try to find a partition where most of the external edges are positive
          # the search is done over the UNDIRECTED version of the graph (smatrix)
          print "Cluster,%IntPosEdges,%IntNegEdges,%ExtPosEdges,%ExtNegEdges"
          for c in xrange(partitionNumber):
            numberOfExtNegEdges = 0
            numberOfExtPosEdges = 0
            numberOfIntNegEdges = 0
            numberOfIntPosEdges = 0
            totalNumberOfIntEdges = 0
            totalNumberOfExtEdges = 0
            for i in xrange(N):
              if cluster[i] == c:
                for j in xrange(N):
                  if smatrix[i,j] != 0:
                    if cluster[j] == c:  # internal edges (within the same cluster)
                      totalNumberOfIntEdges += 1
                      if smatrix[i,j] < 0:
                        numberOfIntNegEdges += 1
                      else:
                        numberOfIntPosEdges += 1
                    else:  # external edges (between clusters)
                      totalNumberOfExtEdges += 1
                      if smatrix[i,j] > 0:
                        numberOfExtPosEdges += 1
                      else:
                        numberOfExtNegEdges += 1
            if totalNumberOfIntEdges > 0:
                PIntPosEdges = float(numberOfIntPosEdges)*100 / totalNumberOfIntEdges
                PIntNegEdges = float(numberOfIntNegEdges)*100 / totalNumberOfIntEdges
            else:
                PIntPosEdges = PIntNegEdges = 0
            if totalNumberOfExtEdges > 0:
                PExtPosEdges = float(numberOfExtPosEdges)*100 / totalNumberOfExtEdges
                PExtNegEdges = float(numberOfExtNegEdges)*100 / totalNumberOfExtEdges
            else:
                PExtPosEdges = PExtNegEdges = 0
            print str(c) + ",%.2f,%.2f,%.2f,%.2f" % (PIntPosEdges, PIntNegEdges, PExtPosEdges, PExtNegEdges)
            
            # internal pos + external pos : "plus mediators" fig 2, according to Doreian et. al
            # internal neg + external pos : "plus mutually hostile mediators" fig 3
            # maybe internal neg + external neg would be "internal subgroup hostility" ???
            if (PIntPosEdges > threshold and PExtPosEdges > threshold):
              PlusMediators.append(str(c) + str(" (PercIntPosEdges = %.2f" % (PIntPosEdges)) + str(" and PercExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
            if (PIntNegEdges > threshold and PExtPosEdges > threshold):
              PlusMutuallyHostileMediators.append(str(c) + str(" (PercIntNegEdges = %.2f" % (PIntNegEdges)) + str(" and PercExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
            if (PIntNegEdges > threshold and PExtNegEdges > threshold):
              InternalSubgroupHostility.append(str(c) + str(" (PercIntNegEdges = %.2f" % (PIntNegEdges)) + str(" and PercExtNegEdges = %.2f" % (PExtNegEdges)) + ")")

          ''' DISABLED
          # the processing below is done over the DIRECTED version of the graph (matrix)
          # detects differential popularity in the directed graph
          # (from a different cluster to the current one)
          for c in xrange(partitionNumber):
            numberOfExtNegEdges = 0
            numberOfExtPosEdges = 0
            totalNumberOfIntEdges = 0
            totalNumberOfExtEdges = 0
            for i in xrange(N):
              if cluster[i] <> c:
                for j in xrange(N):
                  if matrix[i,j] != 0:
                    if cluster[j] == c:  # external edges (between clusters)
                      totalNumberOfExtEdges += 1
                      if matrix[i,j] > 0:
                        numberOfExtPosEdges += 1
                      else:
                        numberOfExtNegEdges += 1
            if totalNumberOfExtEdges > 0:
                PExtPosEdges = float(numberOfExtPosEdges)*100 / totalNumberOfExtEdges
                PExtNegEdges = float(numberOfExtNegEdges)*100 / totalNumberOfExtEdges
            else:
                PExtPosEdges = PExtNegEdges = 0
            
            # external pos : Differential popularity
            if (PExtPosEdges > threshold):
              DifferentialPopularity.append(str(c) + str(" PercExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
          '''

          print "\nProcessing mediation properties...\n"
          # Print results on display
          #result_file = open(folder + "/xpress-summary.csv", "w")
          print "Clusters with plus mediators - undirected (90%+ internal + edges and 90%+ external - edges): "
          if not PlusMediators:
              print "None"
          for elem in (PlusMediators):
              print "%s" % (elem)
          print "\nClusters with plus mutually hostile mediators - undirected (90%+ internal - edges and 90%+ external + edges): "
          if not PlusMutuallyHostileMediators:
              print "None"
          for elem in (PlusMutuallyHostileMediators):
              print "%s" % (elem)
          print "\nClusters with differential popularity ( (any internal) and (90%+ directed external + edges) ): "
          if not DifferentialPopularity:
              print "None"
          for elem in (DifferentialPopularity):
              print "%s" % (elem)         

          print "\n\nClusters with internal subgroup hostility (maybe internal neg + external neg): "
          if not InternalSubgroupHostility:
              print "None   "
          for elem in (InternalSubgroupHostility):
              print "%s" % (elem)
                  
          count = count - 1
	 # end process results


   
if __name__ == "__main__":
   main(sys.argv[1:])
