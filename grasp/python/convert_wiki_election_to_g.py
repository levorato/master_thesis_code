import sys, getopt
import csv
import StringIO
import glob
import os

def main(argv):

   filepath = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["file="])
   except getopt.GetoptError:
      print 'convert_wiki_election_to_g.py -i <file>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convert_wiki_election_to_g.py -i <file>'
         sys.exit()
      elif opt in ("-i", "--file"):
         filepath = arg
   if filepath == '':
      print "Please specify the input file."
      sys.exit(2)
   print 'Input file is ', filepath

   input_file = open(filepath, "r")
   g_file = open(filepath[0:filepath.rfind('.')]+".g", "w")
   content = input_file.read()
   reader = csv.reader(StringIO.StringIO(content), delimiter='\t')
   vertex_a = 0
   vertex_b = 0
   value = 0
   vertex_count = 0
   edge_count = 0
   line_list = list()
   for row in reader:
      line = ''.join(row)
      if len(row) > 0 and row[0][0] != '#' and row[0][0] != 'E' and row[0][0] != 'T' and row[0][0] != 'N':  # ignore comments (#), E, T and N
         if(row[0][0] == 'U'):
            vertex_a = long(row[1])
            vertex_count = max(vertex_count, vertex_a)
         else:  # V
            value = long(row[1])
            vertex_b = long(row[2])
            vertex_count = max(vertex_count, vertex_b)
            line_list.append(str(vertex_a)+"\t"+str(vertex_b)+"\t"+str(value)+"\n")
            edge_count = edge_count + 1
   
   g_file.write(str(vertex_count)+" "+str(edge_count)+"\n")
   for item in line_list:
      g_file.write(item)
   input_file.close
   g_file.close
   print "Output file successfully generated."

if __name__ == "__main__":
   main(sys.argv[1:])

