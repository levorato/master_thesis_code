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
   out_matrix = {}
   additions = list()
   
   for row in reader:
      if row[0][0] == '#':
         line = ''.join(row)
         pos = line.find("Nodes:")
         if pos > 0:
            vertex_count = long(line[pos+6:line.find("Edges")-1].strip())
            print "Vertex count is " + str(vertex_count)
            matrix = array('i',(0,)*vertex_count*vertex_count)
      else:
         vertex_a = long(row[0])
         vertex_b = long(row[1])
         value = long(row[2])
         matrix[vertex_a+vertex_b*vertex_count] = value
   
   # for each pair os nodes, processes edges in both directions
   for i in range(0,vertex_count):
      for j in range(0 ,i-1):
         edge_in = matrix[i+j*vertex_count]
         edge_out = matrix[j+i*vertex_count]
         if edge_in == edge_out:
            if edge_in != None:
               out_matrix[(i,j)] = edge_in
         else:
            additions.append(str(i)+"\t"+str(j)+"\t"+"*")

   for key in sorted(out_matrix.iterkeys()):
      skey = ''.join(key)
      v_a = skey[1:skey.find(",")]
      v_b = skey[skey.find(",")+2:skey.find(")")]
      g_file.write("%s\t%s\t%s\n" % (v_a, v_b, str(out_matrix[key])))
   #for item in line_list:
   #   g_file.write(item)
   input_file.close
   g_file.close
   print "Output file successfully generated."

if __name__ == "__main__":
   main(sys.argv[1:])

