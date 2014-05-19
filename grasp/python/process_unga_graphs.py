# script para processamento do resultado das instancias UNGA
# associa os numeros dos vertices aos nomes dos paises e
# 

# .ccode files location
ccode_dir = '../tests/Instances/UNGA/UNGA_Weighted_SignedGraphs'

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

   global ccode_dir
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


   error_summary = []

   # lookup country full name and iso abbreviation from csv file
   country_full_name = dict()
   with open(ccode_dir + '/' + 'UN_country_names.csv', 'r') as ccdata:
      reader = csv.DictReader(ccdata, delimiter=',')
      for line in reader:
         country_full_name[line["<iso_code>"].strip()] = line["<country_name>"].strip()


   for root, subFolders, files in os.walk(folder):
      # sort dirs and files
      subFolders.sort()
      files.sort()      
      print "Processing folder " + ''.join(root)
      if(len(files)):
         file_list = []
         file_list.extend(glob.glob(root + "/*-result.txt"))
         count = len(file_list) - 1
         
         while count >= 0:
		   print "Processing file " + file_list[count] + "\n"
		   content_file = open(file_list[count], 'r')
                   prefix = file_list[count]
                   result_file = open(prefix[0:prefix.rfind('.txt')] + '-addinfo.txt', "w")
                   result_file.write("Additional results info\n")

                   # lookup of country code data from file
                   graphfile = prefix[0:prefix.rfind('/')]
                   graphfile = graphfile[0:graphfile.rfind('/')]
                   graphfile = graphfile[graphfile.rfind('/')+1:]
                   country = dict()
                   with open(ccode_dir + '/' + graphfile + '.ccode', 'r') as ccdata:
                      reader = csv.DictReader(ccdata, delimiter=' ')
                      for line in reader:
                         country[int(line["<vertex_label>"])] = line["<cname>"]

		   try:
		    content = content_file.read()
		    reader = csv.reader(StringIO.StringIO(content), delimiter='\n')
		    k = 0
                    processed = False
                    clustering_numbers = []
                    clustering_names = []
                    clustering_full_names = []
		    for row in reader:
		       linestring = ''.join(row)
                       if(linestring.startswith(' ')):
                          processed = True
                          vertex_list = linestring[linestring.find('[')+1:linestring.rfind(']')-1].strip()
                          reader2 = csv.reader(StringIO.StringIO(vertex_list), delimiter=' ')
                          for line in reader2:
                             line_out = ''
                             line_out2 = ''
                             line_out3 = ''
                             for vertex in line:
                                line_out += str(country[int(vertex)]) + ","
                                line_out2 += str(vertex) + ","
                                key = str(country[int(vertex)]).strip()
                                if(country_full_name.has_key(key)):
                                   line_out3 += str(country_full_name[key]) + ","
                                else:
                                   line_out3 += str(country[int(vertex)]) + ","
                                   error_summary.append(str(country[int(vertex)]))
                             clustering_names.append(line_out + '\r\n')
                             clustering_numbers.append(line_out2 + '\r\n')
                             clustering_full_names.append(line_out3 + '\r\n')
                       else:
                          if(not processed):
                             result_file.write(str(row[0]) + '\n')
                          else:
                             break

                    for line in clustering_full_names:
                       result_file.write(line)
                    result_file.write('\r\n')
                    for line in clustering_names:
                       result_file.write(line)
                    result_file.write('\r\n')
                    for line in clustering_numbers:
                       result_file.write(line)

		   finally:
		    content_file.close()
                    result_file.close()
                   count = count - 1
	 # end process results

   print "\nSuccessfully processed all graph files.\n"
   if(len(error_summary) > 0):
      print "WARNING: The following country full names could not be found: "
      items = ''
      for item in set(error_summary):
         items += item + ", "
      print items

if __name__ == "__main__":
   main(sys.argv[1:])
