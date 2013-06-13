import sys, getopt
import csv
import StringIO
import glob
import os

def main(argv):

   folder = ''
   filter = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["folder=","filter="])
   except getopt.GetoptError:
      print 'process_output.py -i <folder> -f <filter>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'process_output.py -i <folder> -f <filter>'
         sys.exit()
      elif opt in ("-i", "--folder"):
         folder = arg
      if opt in ("-f", "--filter"):
         filter = arg
   print 'Input dir is ', folder
   print 'File filter is ', filter

   all_files_summary = dict()

   for root, subFolders, files in os.walk(folder):
      print "Processing folder " + ''.join(root)
      if(len(files) and ''.join(root) != folder):
         file_list = []
         file_list.extend(glob.glob(root + "/" + filter))
         count = len(file_list) - 1
         
         text_file = open(root + "/summary.txt", "w")
         filename = (root[:root.rfind("/")])
         datetime = root[root.rfind("/")+1:]
         filename = filename[filename.rfind("/")+1:]
         text_file.write("Summary for graph file: %s\n"%filename)
         best_value = 1000000L
         best_pos_value = 0
         best_neg_value = 0
         best_K = 0
         best_iteration = 0
         best_time = 0
         best_param = ''

         while count >= 0:
            with open(file_list[count], 'r') as content_file:
               content = content_file.read()

            reader = csv.reader(StringIO.StringIO(content), delimiter=',')
            for row in reader:
               linestring = ''.join(row)
               if linestring.startswith( 'Best value' ):
                  column = []
                  for col in row:
                     column.append(col)
                  
                  filepath = ''.join(file_list[count])
                  text_file.write(filepath[filepath.rfind("/")+1:] + ' ' + linestring + '\n')
                  value = long(column[1])
                  pos_value = long(column[2])
                  neg_value = long(column[3])
                  K = long(column[4])
                  iteration = long(column[5]) 
                  time = float(column[6])
                  if value < best_value :
                     best_value = value
                     best_pos_value = pos_value
                     best_neg_value = neg_value
                     best_K = K
                     best_iteration = iteration
                     best_time = time
                     best_param = filepath[filepath.rfind("/")+1:]
                  elif value == best_value and iteration < best_iteration :
                     best_K = K
                     best_iteration = iteration
                     best_time = time
                     best_param = filepath[filepath.rfind("/")+1:]
                     best_pos_value = pos_value
                     best_neg_value = neg_value
            count = count - 1

         text_file.close()
         all_files_summary[filename+"/"+datetime] = str(best_value)+", "+str(pos_value)+", "+str(neg_value)+", "+str(best_K)+", "+str(iteration)+", "+str(best_time)+", "+best_param

   result_file = open(folder + "/summary.txt", "w")
   print "Filename, I(P), I(P)+, I(P)-, k, Iter, Time(s), Params"
   result_file.write("Filename, I(P), I(P)+, I(P)-, k, Iter, Time(s), Params")
   for key in sorted(all_files_summary.iterkeys()):
      print "%s, %s" % (key, all_files_summary[key])
      result_file.write("%s: %s\n" % (key, all_files_summary[key]))
   result_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
