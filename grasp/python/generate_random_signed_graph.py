#*****************************************************
# Instance generator for random symmetric signed graphs
# INPUT:
#       n number of vertices
#       d graph density ( 2m/n(n-1) )
#       dn negative percentual (|E-|/m)
#       
#       where m=|E+|+|E-|
# OUTPUT
#       undirected graph with NO parallel edges
# OUTPUT FORMAT
#       <n> <|E+ \cup E-|>
#       <i> <j> <s>
#       .
#       .
#       .
# ****************************************************/

import sys
import csv
import StringIO
import os
import math
import argparse
import collections
import random

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __str__(self):
        return "range(" + str(self.start) + ", " + str(self.end) + ")"

edge = collections.namedtuple('Edge', ['x', 'y'])

# Returns a random edge from the graph
def pick_random_edge(n):
    vertex_list = range(n)
    a = random.choice(vertex_list)
    b = random.choice(vertex_list)
    while(b == a):
        b = random.choice(vertex_list)
    return edge(a, b)

def generate_random_graph(n, d, dneg):

   M = [[0 for x in xrange(n)] for x in xrange(n)]
   
   m = 0     # number of edges
   mneg = 0  # number of neg edges
   mpos = 0  # number of pos edges
   
   m_min = int(math.ceil(d*(n*(n-1)/2.0)))
   mneg_min = int(math.ceil(dneg*m_min))
   
   while(mneg < mneg_min):
      e = pick_random_edge(n)
      i = e.x
      j = e.y
      if(M[i][j] == 0):
         mneg += 1
         m += 1
         M[i][j] = -1
         M[j][i] = -1
   
   while(m < m_min):
      e = pick_random_edge(n)
      i = e.x
      j = e.y
      if(M[i][j] == 0):
         mpos += 1
         m += 1
         M[i][j] = 1
         M[j][i] = 1
   
   # Writes output file in XPRESS Mosel format (.mos)
   # --------------------------------------------------
   # Stores graph file in output folder
   # Vertex numbers start at 1
   filename_prefix = "file_" + str(n) + "_" + str(d) + "_" + str(dneg)
   filename = filename_prefix + ".g"
   directory = "output-random-sg"
   if not os.path.exists(directory):
        os.makedirs(directory)
   with open(directory+"/"+filename, "w") as g_file:
        # file header
        g_file.write('people: {0}\r\n\r\n'.format(str(n)))
        g_file.write('VarErr: 0.5\r\n\r\n')
        vertex_names = ""
        for x in xrange(1,n+1):
             vertex_names += "V" + str(x) + " "
        g_file.write('Names: [{0}]\r\n\r\n'.format(vertex_names))
        # graph contents
        edge_list = ""
        for i in xrange(n):
             for j in xrange(n):
                   if(M[i][j] <> 0):
                       edge_list += '({0},{1}){2} \r\n'.format(str(i+1), str(j+1), str(M[i][j]))
        g_file.write('Mrel: [ {0}]'.format(edge_list))

   print "Number of negative edges: {0} ({1}%), number of positive edges: {2} ({3}%), total edges: {4} ({5}%).".format(str(mneg), str(100*mneg/m), str(mpos), str(100*mpos/m), str(m), str(100*m/(n*(n-1)/2)))
   print "Output file successfully generated."

def main(argv):

   parser = argparse.ArgumentParser()
   parser.add_argument('-n', type=int, required=True, help="number of nodes")
   parser.add_argument('-d', type=float, required=True, choices=[Range(0.0, 1.0)], help="density of edges of the graph")
   parser.add_argument('-nd', type=float, required=True, choices=[Range(0.0, 1.0)], help="negative density of edges of the graph")
   args = parser.parse_args()
     
   # Call random graph generation procedure
   generate_random_graph(args.n, args.d, args.nd)

if __name__ == "__main__":
   main(sys.argv[1:])


