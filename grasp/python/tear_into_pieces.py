import sys, getopt
import csv
import StringIO
import glob
import os
import math
from array import *

def main(argv):

   filepath = ''
   size = 10000
   try:
      opts, args = getopt.getopt(argv,"hi:n:",["file=","size="])
   except getopt.GetoptError:
      print 'tear_into_pieces.py -i <file> -n <size>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'tear_into_pieces.py -i <file> -n <size>'
         sys.exit()
      elif opt in ("-i", "--file"):
         filepath = arg
      elif opt in ("-n", "--size"):
         size = int(arg)
   if filepath == '':
      print "Please specify the input file."
      sys.exit(2)
   print 'Input .g file is ', filepath

   input_file = open(filepath, "r")
   content = input_file.read()
   reader = csv.reader(StringIO.StringIO(content), delimiter='\t')
   vertex_a = 0
   vertex_b = 0
   value = 0
   vertex_count = 0
   edge_count = array('l')
   additions = array('l')
   g_file = []
   file_content = []
   directory = ""
      
   for row in reader:
      if len(row) == 1:  # file header with num of vertices and edges
         x = row[0]
         vertex_count = int(x[:x.find(" ")])
         num_pieces = int(math.ceil(vertex_count / size))
         print "Breaking the graph file into " + str(num_pieces) + " pieces of size " + str(size)
         edge_count = range(num_pieces+1)
         additions = range(num_pieces+1)
         for i in range(0, num_pieces+1):
            path = filepath[0:filepath.rfind('/')]
            filename = filepath[filepath.rfind('/')+1:filepath.rfind('.')]+"-size"+str(size)
            directory = path+"/"+filepath[filepath.rfind('/')+1:filepath.rfind('.')]
            if not os.path.exists(directory):
               os.makedirs(directory)
            g_file.append(open(directory+"/"+filename+"-part"+str(i)+".g", "w"))
            file_content.append("")
            edge_count[i] = 0
            additions[i] = 0 
      else:
         vertex_a = int(row[0])
         vertex_b = int(row[1])
         div_a = int(math.floor(vertex_a / size))
         div_b = int(math.floor(vertex_b / size))
         if div_a == div_b:  # ignore edges between different pieces
            if(row[2] == "*"):
               additions[div_a] = additions[div_a] + 1
            else:
               value = int(row[2])
               edge_count[div_a] = edge_count[div_a] + 1
               # print "div_a = " + str(div_a)
               # apply offset to vertex numbers (to always start from zero)
               file_content[div_a] += str(vertex_a-size*div_a)+"\t"+str(vertex_b-size*div_a)+"\t"+str(value)+"\n"
     
   additions_file = open(directory+"/additions.txt", "w")
   for i in range(0, num_pieces+1):
      v_count = 0
      if(i < num_pieces):
         v_count = size
      else:
         v_count = vertex_count - (size * num_pieces)
      # prints file header
      g_file[i].write(str(v_count)+" "+str(edge_count[i])+"\n")
      # prints file contents
      g_file[i].write(file_content[i])
      # writes the additions to file
      additions_file.write("Add "+str(additions[i])+" to part "+str(i)+"\n")

   additions_file.close()
   input_file.close()
   for i in range(0, num_pieces+1):
      g_file[i].close()
   print "Output files successfully generated."

if __name__ == "__main__":
   main(sys.argv[1:])

