import sys, getopt
import csv
import StringIO
import glob
import os
from array import *

def main(argv):

   filepath = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["file="])
   except getopt.GetoptError:
      print 'convert_directed_to_undirected_g.py -i <file>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convert_directed_to_undirected_g.py -i <file>'
         sys.exit()
      elif opt in ("-i", "--file"):
         filepath = arg
   if filepath == '':
      print "Please specify the input file."
      sys.exit(2)
   print 'Input .g file is ', filepath

   input_file = open(filepath, "r")
   g_file = open(filepath[0:filepath.rfind('.')]+"-undirected.g", "w")
   content = input_file.read()
   reader = csv.reader(StringIO.StringIO(content), delimiter='\t')
   vertex_a = 0
   vertex_b = 0
   value = 0
   vertex_count = 0
   edge_count = 0
   graph = {}
   out_matrix = {}
   additions = list()
   
   for row in reader:
      if row[0][0] == '#':
         line = ''.join(row)
         pos = line.find("Nodes:")
         if pos > 0:
            vertex_count = int(line[pos+6:line.find("Edges")-1].strip())
            print "Vertex count is " + str(vertex_count)
      elif len(row) == 1:  # file header with num of vertices and edges
         x = row[0]
         vertex_count = int(x[:x.find(" ")])
      else:
         vertex_a = int(row[0])
         vertex_b = int(row[1])
         value = int(row[2])
         key = (vertex_a,vertex_b)
         if not graph.has_key(key):
            graph[key] = value
         else:
            if graph[key] != value:
               print "WARNING: duplicated directed edge with different value!"
         rkey = (vertex_b,vertex_a)
         value2 = graph.get(rkey)
         if value2 != None:
            if key != rkey:  # ignore self-relationships
               if value == value2:
                  if vertex_a < vertex_b:
                     out_matrix[key] = value
                  else:
                     out_matrix[rkey] = value
               else:
                  additions.append(str(vertex_a)+"\t"+str(vertex_b)+"\t"+"*\n")
               del graph[key]
               del graph[rkey]
               edge_count = edge_count + 1
            else:
               del graph[key]

   edge_count = edge_count + len(graph)
   
   # prints file header
   g_file.write(str(vertex_count)+" "+str(edge_count)+"\n")
   # for each remaining edges...
   for key in sorted(graph.iterkeys()):
      skey = ''.join(str(key))
      v_a = skey[1:skey.find(",")]
      v_b = skey[skey.find(",")+2:skey.find(")")]
      g_file.write("%s\t%s\t%s\n" % (str(v_a).rstrip("L"), str(v_b).rstrip("L"), str(graph[key])))
      print "%s\t%s\t%s" % (str(v_a).rstrip("L"), str(v_b).rstrip("L"), str(graph[key]))

   for key in sorted(out_matrix.iterkeys()):
      skey = ''.join(str(key))
      v_a = skey[1:skey.find(",")]
      v_b = skey[skey.find(",")+2:skey.find(")")]
      g_file.write("%s\t%s\t%s\n" % (str(v_a).rstrip("L"), str(v_b).rstrip("L"), str(out_matrix[key])))
      print "%s\t%s\t%s" % (str(v_a).rstrip("L"), str(v_b).rstrip("L"), str(out_matrix[key]))

   #for item in additions:
   #   g_file.write(item)
   print "Add " + str(len(additions))

   input_file.close
   g_file.close
   print "Output file successfully generated."

if __name__ == "__main__":
   main(sys.argv[1:])

