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
   matrix = []

   # RCC results
   all_files_summary = dict()

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
          # reads graph file, filling adjecency matrix
          print 'Reading graph file: {0}'.format(str(instances_path + '/' + graphfile))
          with open(instances_path + '/' + graphfile, 'rb') as csvfile:
             dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=" \t")
             csvfile.seek(0)
             r = csv.reader(csvfile, dialect)
             line1=r.next()
             values = [col for col in line1]
             if(values[0].find(' ') > 0):
                N = int(values[0][:values[0].find(' ')])
             else:
                N = int(values[0])
             matrix = [[0 for x in xrange(N)] for x in xrange(N)]
             for line in r:
                i = int(line[0])
                j = int(line[1])
                w = float(line[2])
                matrix[i][j] = w

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
       			
               print str(matrix[1960][1988])
          finally:
              content_file.close()
          count = count - 1
	 # end process results

   print "\nProcessing RCC Results...\n"
   # Print results on display
   #result_file = open(folder + "/xpress-summary.csv", "w")
   #print "Instance, n, e, e+, e-, I(P), I(P)+, I(P)-, k, Time, Status"
   #result_file.write("Instance, n, e, e+, e-, I(P), I(P)+, I(P)-, k, Time, Status\n")
   #for key in sorted(all_files_summary.iterkeys()):
   #   print "%s, %s" % (key, all_files_summary[key])
   #   result_file.write("%s, %s\n" % (key, all_files_summary[key]))
   #result_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
