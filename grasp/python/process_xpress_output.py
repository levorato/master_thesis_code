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
      print 'process_xpress_output.py --folder <folder> --filter <filter>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'process_xpress_output.py --folder <folder> --filter <filter>'
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
         file_list.extend(glob.glob(root + "/*result.dat"))
         count = len(file_list) - 1
         print 'Found {0} files\n'.format(str(count))
         
         while count >= 0:
		   print "Processing file " + file_list[count] + "\n"
		   filename = str(file_list[count])
		   filename = filename[filename.rfind('/') + 1:filename.rfind('-')]
		   content_file = open(file_list[count], 'r')
		   try:
		    content = content_file.read().splitlines()
		    status = ' - '
		    fo = ' - '
		    pos_imb = neg_imb = ' - '
		    k = ' - '
		    for row in content:
       			linestring = ''.join(row)
       			# Campos a serem capturados do arquivo de saida do xpress:
			    #x Total de vertices: 100
				#x Total de arestas: 990
				#x Total de arestas pos: 792
				#x Total de arestas neg: 198
				# Problem status: Optimum found
				#x Valor F.Objetivo: 198
				#x Total Erro positivo entre particoes: 0
				#x Total Erro negativo nas   particoes: 198
				#x Total de clusters gerados: 1
				#x Problem status: Unfinished
				#x *** Search unfinished ***    Time: 3817
				#x *** Search completed ***     Time:  2789
		
       			if linestring.startswith( 'Total de vertices:' ):          
		          n = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total de arestas:' ):
			  	  e = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total de arestas pos:' ):
			  	  e_pos = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total de arestas neg:' ):
			  	  e_neg = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Problem status:' ):
			  	  status = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Valor F.Objetivo:' ):
			  	  fo = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total Erro positivo entre particoes:' ):
			  	  pos_imb = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total Erro negativo nas   particoes:' ):
			  	  neg_imb = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Total de clusters gerados:' ):
			  	  k = linestring[linestring.rfind(":")+2:]
       			elif linestring.startswith( 'Tempo de Execu' ): #  '*** Search'
			  	  time_spent = linestring[linestring.rfind("Total:")+6:]
				  if 'Nodes:' in str(time_spent):
					time_spent = time_spent[:time_spent.rfind("Nodes:")]
		   finally:
		    content_file.close()
		   all_files_summary[filename] = str(n).strip()+", "+str(e).strip()+", "+str(e_pos).strip()+", "+str(e_neg).strip()+", "+str(fo).strip()+", "+str(pos_imb).strip()+", "+str(neg_imb).strip()+", "+str(k).strip()+", "+str(time_spent).strip()+", "+str(status).strip()
		   count = count - 1
	 # end process results

   print "\nProcessing Xpress Results...\n"
   # Print results on display
   result_file = open(folder + "/xpress-summary.csv", "w")
   print "Instance, n, e, e+, e-, I(P), I(P)+, I(P)-, k, Time, Status"
   result_file.write("Instance, n, e, e+, e-, I(P), I(P)+, I(P)-, k, Time, Status\n")
   for key in sorted(all_files_summary.iterkeys()):
      print "%s, %s" % (key, all_files_summary[key])
      result_file.write("%s, %s\n" % (key, all_files_summary[key]))
   result_file.close()

if __name__ == "__main__":
   main(sys.argv[1:])
