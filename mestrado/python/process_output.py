import sys, getopt
import csv
import StringIO
import glob
import os

def main(argv):

   folder = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["folder=","ofile="])
   except getopt.GetoptError:
      print 'process_output.py -i <folder>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'process_output.py -i <folder>'
         sys.exit()
      elif opt in ("-i", "--folder"):
         folder = arg
   print 'Input dir is ', folder

   all_files_summary = dict()

   for root, subFolders, files in os.walk(folder):
      print "Processing folder " + ''.join(root)
      if(len(files) and ''.join(root) != folder):
         file_list = []
         file_list.extend(glob.glob(root + "/*.csv"))
         count = len(file_list) - 1
         
         text_file = open(root + "/summary.txt", "w")
         filename = root[root.rfind("file"):root.rfind("/")]
         text_file.write("Summary for graph file: %s\n"%filename)

         while count >= 0:
            with open(file_list[count], 'r') as content_file:
               content = content_file.read()
            reader = csv.reader(StringIO.StringIO(content), delimiter=',')
            best_value = long(1000000)
            best_iteration = long(0)
            best_time = long(0)
            best_param = ''
            for row in reader:
               linestring = ''.join(row)
               if linestring.startswith( 'Best value' ):
                  filepath = ''.join(file_list[count])
                  text_file.write(filepath[filepath.rfind("/")+1:] + ' ' + linestring + '\n')
                  value = long(linestring[linestring.find(":")+2:linestring.find("Iteration")-1])
                  iteration = long(linestring[linestring.find("Iteration:")+11:linestring.rfind("Time")-1]) 
                  time = long(linestring[linestring.rfind(":")+2:])
                  if value < best_value :
                     best_value = value
                     best_iteration = iteration
                     best_time = time
                     best_param = filepath[filepath.rfind("/")+1:]
                  elif value == best_value and iteration < best_iteration :
                     best_iteration = iteration
                     best_time = time
                     best_param = filepath[filepath.rfind("/")+1:]
            count = count - 1

         text_file.close()
         all_files_summary[filename] = (str(best_value), str(iteration), str(best_time), best_param)

   result_file = open(folder + "/summary.txt", "w")
   for key in sorted(all_files_summary.iterkeys()):
      print "%s: %s" % (key, all_files_summary[key])
      result_file.write("%s: %s\n" % (key, all_files_summary[key]))
   result_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
