import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re

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
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["folder=","filter="])
   except getopt.GetoptError:
      print 'extract_random_summary.py --folder <folder> --filter <filter>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'extract_random_summary.py --folder <folder> --filter <filter>'
         sys.exit()
      elif opt in ("-i", "--folder"):
         folder = arg
      if opt in ("-f", "--filter"):
         filter = arg
   print 'Input dir is ', folder
   print 'File filter is ', filter

   # CC results
   all_files_summary = dict()

   for root, subFolders, files in os.walk(folder):
      # sort dirs and files
      subFolders.sort()
      files.sort()      
      print 'Found {0} files\n'.format(str(len(files)))
      print "Processing folder " + ''.join(root)
      if(len(files)):
         file_list = []
         file_list.extend(glob.glob(root + "/*-info.txt"))
         count = len(file_list) - 1
         print 'Found {0} files\n'.format(str(count))
         
         while count >= 0:
		   print "Processing file " + file_list[count] + "\n"
		   content_file = open(file_list[count], 'r')
		   try:
		    content = content_file.read()
		    reader = csv.reader(StringIO.StringIO(content), delimiter=' ')
		    k = 0
		    for row in reader:
		       linestring = ''.join(row)
		       column = []
                       for col in row:
                          column.append(col)
		       if linestring.startswith( 'I(P)' ):          
		          filepath = ''.join(file_list[count])
		          filename = file_list[count]  
		          filename = filename[filename.rfind("/")+1:]
			  value = float(column[2])
		       elif linestring.startswith( 'c(clusters)' ):
			  k = int(column[2])
		    count = count - 1
		   finally:
		    content_file.close()
		   all_files_summary[filename] = str(value)+", "+str(k)
	 # end process results

   print "\nProcessing Random Graph Generation Results...\n"
   # Print results on display
   result_file = open(folder + "/summary.txt", "w")
   print "Instance, I(P), k"
   result_file.write("Filename, I(P), k\n")
   for key in sorted(all_files_summary.iterkeys()):
      print "%s, %s" % (key, all_files_summary[key])
      result_file.write("%s, %s\n" % (key, all_files_summary[key]))
   result_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
