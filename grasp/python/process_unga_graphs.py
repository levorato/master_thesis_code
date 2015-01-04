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
import math

import HTML

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
   instances_path = ''
   threshold = float(90.0)  # percentual threshold for the quantity of relashionships
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["folder=","filter=","instancespath="])
   except getopt.GetoptError:
      print 'process_unga_graphs.py --folder <folder> --filter <filter> --instancespath <instances_path>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'process_unga_graphs.py --folder <folder> --filter <filter> --instancespath <instances_path>'
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
   print 'Mediation threshold is ', threshold

   error_summary = []
   cc_result_summary = dict()
   cc_imbalance_summary = dict()
   rcc_result_summary = dict()
   rcc_imbalance_summary = dict()
   cc_cluster_imb_matrix = dict()
   rcc_cluster_imb_matrix = dict()
   matrix = []
   neg_edge_sum = 0.0
   pos_edge_sum = 0.0
   # lists of interesting clusters (with mediation properties)
   PlusMediators = dict()
   PlusMutuallyHostileMediators = dict()
   InternalSubgroupHostility = dict()

   # lookup country full name and iso abbreviation from csv file
   country_full_name = dict()
   with open(ccode_dir + '/' + 'UN_country_names.csv', 'r') as ccdata:
      reader = csv.DictReader(ccdata, delimiter=',')
      for line in reader:
         country_full_name[line["<iso_code>"].strip()] = line["<country_name>"].strip()

   previous_filename = ''
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
                   result_file_path = prefix[0:prefix.rfind('.txt')]
                   result_file = open(result_file_path + '-addinfo.txt', "w")
                   result_file.write("Additional results info\n")
                   result_file_name = result_file_path[result_file_path.rfind('/')+1:]
                   
                   # lookup of country code data from file
                   graphfile = prefix[0:prefix.rfind('/')]
                   graphfile = graphfile[0:graphfile.rfind('/')]
                   graphfile = graphfile[graphfile.rfind('/')+1:]
                   country = dict()
                   with open(ccode_dir + '/' + graphfile + '.ccode', 'r') as ccdata:
                      reader = csv.DictReader(ccdata, delimiter=' ')
                      for line in reader:
                         country[int(line["<vertex_label>"])] = line["<cname>"]

		   # process graph file contents and obtains positive and negative edge weight sum
		   if(result_file_name <> previous_filename):
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
		       #print 'N is {0}'.format(str(N))
		       matrix = [[0 for x in xrange(N)] for x in xrange(N)]
		       pos_edge_sum = 0.0
		       neg_edge_sum = 0.0		    
		       for line in r:
		          i = int(line[0])
		          j = int(line[1])
		          w = float(line[2])
		          matrix[i][j] = w
		          #print '({0}, {1}) = {2}\n'.format(str(i), str(j), str(w))
		          if(w < 0):
		             neg_edge_sum += math.fabs(w)
		          else:
		             pos_edge_sum += w
		   print 'neg_edge_sum = {0}, pos_edge_sum = {1}'.format(str(neg_edge_sum), str(pos_edge_sum))

		   try:  # process *-result.txt file
		    content = content_file.read()
		    reader = csv.reader(StringIO.StringIO(content), delimiter='\n')
		    k = 0
                    processed = False
                    processed_imb_analysis = False
                    clustering_numbers = []
                    clustering_names = []
                    clustering_full_names = []
                    vertex_contrib = [[float(0.0), float(0.0)] for x in xrange(N)]
                    pos_contrib_sum = float(0.0)
                    neg_contrib_sum = float(0.0)
		    started_cluster_analysis = False
		    processed_imb_analysis = False
		    started_imb_analysis = False
		    partitionNumber = 0
		    cluster = [int(0) for x in xrange(N)]
		    for row in reader:
		       linestring = ''.join(row)
                       if(linestring == ''):
                          continue
                       if(linestring.startswith(' ')):
                          processed = True
                          vertex_list = linestring[linestring.find('[')+1:linestring.rfind(']')-1].strip()
                          reader2 = csv.reader(StringIO.StringIO(vertex_list), delimiter=' ')
                          for line in reader2:
                             line_out = ''
                             line_out2 = ''
                             line_out3 = ''
                             for vertex in line:
                                line_out += str(country[int(vertex)]) + ", "
                                line_out2 += str(vertex) + ", "
                                key = str(country[int(vertex)]).strip()
                                if(country_full_name.has_key(key)):
                                   line_out3 += str(country_full_name[key]) + ", "
                                else:
                                   line_out3 += str(country[int(vertex)]) + ", "
                                   error_summary.append(str(country[int(vertex)]))
                                cluster[int(vertex)] = partitionNumber
                             clustering_names.append(line_out + '\r\n')
                             clustering_numbers.append(line_out2 + '\r\n')
                             clustering_full_names.append(line_out3 + '\r\n')
                          partitionNumber += 1
                       else:
                          if(not processed):
                             result_file.write(str(row[0]) + '\n')
                          else:
                             # processes imbalance analysis data (out edges contribution)
                             #reader3 = csv.reader(StringIO.StringIO(linestring), delimiter=',')
                             #line = reader3.next()
                             if(linestring.startswith('Cluster')):
                                if(not started_cluster_analysis):
                                   started_cluster_analysis = True
				# initialize matrix with cluster contribution
				numberOfClusters = len(clustering_full_names)
				clusterImbMatrix = [[float(0.0) for x in xrange(numberOfClusters)] for x in xrange(numberOfClusters)]
				line = 0
                             else:
                                if(not processed_imb_analysis):
				   if(linestring.startswith('Vertex')):
				       continue
				   if(linestring.startswith('Imbalance')):
				       if(not started_imb_analysis):
				          started_imb_analysis = True
				       else:
					  processed_imb_analysis = True
				       continue
                                   # processes vertex contribution to imbalance
                                   #reader3 = csv.reader(StringIO.StringIO(linestring), delimiter=',')
                                   #line = reader3.next()
                                   # print 'line: {0}'.format(str(linestring))
                                   vertex = int(linestring[:linestring.find(',')])
                                   pos_contrib = float(linestring[linestring.find(',')+1:linestring.rfind(',')])
                                   neg_contrib = float(linestring[linestring.rfind(',')+1:])
                                   vertex_contrib[vertex][0] = pos_contrib
                                   vertex_contrib[vertex][1] =  neg_contrib
                                   pos_contrib_sum += pos_contrib
                                   neg_contrib_sum += neg_contrib
                                   # print '{0} {1}'.format(str(neg_contrib_sum), str(neg_contrib))
                                if(started_cluster_analysis):
                                   # reads matrix with cluster contribution to imbalance
				   reader3 = csv.reader(StringIO.StringIO(linestring), delimiter=',')
				   column = 0
                                   for elem in reader3.next():
					if(str(elem).strip() != ''):
						clusterImbMatrix[line][column] = float(str(elem))
						column += 1
				   line += 1

		    if(result_file_name.startswith('rcc')):  # process mediation info from rcc result
			  # for each partition, tries to find a partition where most of the external edges are positive
			  #print "Cluster,%IntPosEdges,%IntNegEdges,%ExtPosEdges,%ExtNegEdges"
			  for c in xrange(partitionNumber):
			    numberOfExtNegEdges = 0
			    numberOfExtPosEdges = 0
			    numberOfIntNegEdges = 0
			    numberOfIntPosEdges = 0
			    totalNumberOfIntEdges = 0
			    totalNumberOfExtEdges = 0
			    for i in xrange(N):
			      if cluster[i] == c:  # i is in cluster c
				for j in xrange(N):
				  if matrix[i][j] != 0:  # edge (i, j)
				    if cluster[j] == c:  # internal edges (within the same cluster)
				      totalNumberOfIntEdges += 1
				      if matrix[i][j] < 0:
				        numberOfIntNegEdges += 1
				      else:
				        numberOfIntPosEdges += 1
				    else:  # external edges (between clusters)
				      totalNumberOfExtEdges += 1
				      if matrix[i][j] > 0:
				        numberOfExtPosEdges += 1
				      else:
				        numberOfExtNegEdges += 1
			    if totalNumberOfIntEdges > 0:
				PIntPosEdges = float(numberOfIntPosEdges)*100 / totalNumberOfIntEdges
				PIntNegEdges = float(numberOfIntNegEdges)*100 / totalNumberOfIntEdges
			    else:
				PIntPosEdges = PIntNegEdges = 0
			    if totalNumberOfExtEdges > 0:
				PExtPosEdges = float(numberOfExtPosEdges)*100 / totalNumberOfExtEdges
				PExtNegEdges = float(numberOfExtNegEdges)*100 / totalNumberOfExtEdges
			    else:
				PExtPosEdges = PExtNegEdges = 0
			    #print str(c) + ",%.2f,%.2f,%.2f,%.2f" % (PIntPosEdges, PIntNegEdges, PExtPosEdges, PExtNegEdges)
			    
			    # internal pos + external pos : "plus mediators" fig 2, according to Doreian et. al
			    # internal neg + external pos : "plus mutually hostile mediators" fig 3
			    # maybe internal neg + external neg would be "internal subgroup hostility" ???
			    if (PIntPosEdges > threshold and PExtPosEdges > threshold):
			      PlusMediators[graphfile] = ("Cluster " + str(c+1) + str(" (%%IntPosEdges = %.2f" % (PIntPosEdges)) + str(" and %%ExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
			    if (PIntNegEdges > threshold and PExtPosEdges > threshold):
			      PlusMutuallyHostileMediators[graphfile] = ("Cluster " + str(c+1) + str(" (%%IntNegEdges = %.2f" % (PIntNegEdges)) + str(" and %%ExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
			    if (PIntNegEdges > threshold and PExtNegEdges > threshold):
			      InternalSubgroupHostility[graphfile] = ("Cluster " + str(c+1) + str(" (%%IntNegEdges = %.2f" % (PIntNegEdges)) + str(" and %%ExtNegEdges = %.2f" % (PExtNegEdges)) + ")")

		    # export cluster imb matrix to html
		    matrixline = ['Cluster']
		    for line in xrange(1, numberOfClusters+1):
			matrixline.append('%d' % line)
		    matrixline.append('Sum')
		    t = HTML.Table(header_row=matrixline)
		    for line in xrange(0, numberOfClusters):
			matrixline = ['<b>' + str(line+1) + '</b>']
			sum = float(0.0)
			for column in xrange(0, numberOfClusters):
				if(clusterImbMatrix[line][column] > 0):
					matrixline.append('%.4f' % clusterImbMatrix[line][column])
				else:
					matrixline.append('<font color=\"red\">%.4f</font>' % clusterImbMatrix[line][column])
				sum += abs(clusterImbMatrix[line][column])
			matrixline.append('%.4f' % sum)
	      		t.rows.append(matrixline)
		    if(result_file_name.startswith('cc')):
		    	cc_cluster_imb_matrix[graphfile] = str(t)
		    else:
		     	rcc_cluster_imb_matrix[graphfile] = str(t)
		    # end process *-result.txt files


                    # process *-Node0-*-iterations.csv file
                    csv_file_list = []
                    if(result_file_name.startswith('cc')):
                        prfx = 'CC'
                    else:
                        prfx = 'RCC'
                        pos_imbalance = pos_contrib_sum
                        neg_imbalance = neg_contrib_sum
                    csv_file_list.extend(glob.glob(root + '/' + prfx + '-Node0-*-iterations.csv'))
                    with open(csv_file_list[0], 'r') as r_csv_file:
                       r2 = csv.reader(r_csv_file)
		       line = r2.next()
                       while (not line[0].startswith('Best')):
                          line = r2.next()
                       # extracts total, positive and negative imbalance of the best result
                       imbalance = float(line[1])
                       if(result_file_name.startswith('cc')):
                          pos_imbalance = float(line[2])
                          neg_imbalance = float(line[3])
                          print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance), str(neg_imbalance))
                       else:  # for RCC problem the numbers are internal and external imbalance values
                          int_imbalance = float(line[2])
                          ext_imbalance = float(line[3])
                          print 'imbalances {0} (int){1} (ext){2}'.format(str(imbalance), str(int_imbalance), str(ext_imbalance))
                          print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance), str(neg_imbalance))

                    for line in clustering_full_names:
                       result_file.write(line)
                    result_file.write('\r\n')
                    for line in clustering_names:
                       result_file.write(line)
                    result_file.write('\r\n')
                    for line in clustering_numbers:
                       result_file.write(line)

                    if(result_file_name.startswith('cc')):
                       cc_result_summary[graphfile] = clustering_full_names
                       cc_imbalance_summary[graphfile] = [str(imbalance), str(pos_edge_sum+neg_edge_sum), str(round(100*imbalance/(pos_edge_sum+neg_edge_sum), 2)), str(pos_imbalance), str(pos_edge_sum), str(round(100*pos_imbalance/pos_edge_sum, 2)), str(neg_imbalance), str(neg_edge_sum), str(round(100*neg_imbalance/neg_edge_sum, 2))]
                    else:  # rcc
                       rcc_result_summary[graphfile] = clustering_full_names
                       rcc_imbalance_summary[graphfile] = [str(imbalance), str(pos_edge_sum+neg_edge_sum), str(round(100*imbalance/(pos_edge_sum+neg_edge_sum), 2)), str(pos_imbalance), str(pos_edge_sum), str(round(100*pos_imbalance/pos_edge_sum, 2)), str(neg_imbalance), str(neg_edge_sum), str(round(100*neg_imbalance/neg_edge_sum, 2)), str(int_imbalance), str(round(100*int_imbalance/(pos_edge_sum+neg_edge_sum), 2)), str(ext_imbalance), str(round(100*ext_imbalance/(pos_edge_sum+neg_edge_sum), 2))]
		   finally:
		    content_file.close()
                    result_file.close()
                   count = count - 1
                   previous_filename = graphfile
	 # end process results

   # TODO implementar tabela com os valores de contribuicao pos/neg entre clusters

   print "\nSuccessfully processed all graph files.\n"
   if(len(error_summary) > 0):
      print "WARNING: The following country full names could not be found: "
      items = ''
      for item in set(error_summary):
         items += item + ", "
      print items

   # exports summary to HTML file
   # Exporting CC and RCC results
   mhtml = '<!DOCTYPE html>\n<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://genshi.edgewall.org/">\n<head><title>UNGA CC and RCC Summary</title>'
   mhtml += '</head><body class="index">'
   mhtml += '<h1>UNGA Voting Data - CC and RCC Imbalance Analysis</h1>'
   startyear = 1946
   partial_html_cc = []
   partial_html_rcc = []
   
   for key, value in sorted(cc_result_summary.iteritems()):
      html = '<h2>CC - ' + str(key)  + ' (' + str(startyear) + ')</h2><p>'

      imbalance_array = cc_imbalance_summary[key]
      t = HTML.Table(header_row=['I(P)', 'Edge sum', 'I(P)%', 'I(P)+/I(P)ext', 'Pos edge sum', 'I(P)+%', 'I(P)-/I(P)int', 'Neg edge sum', 'I(P)-%'])
      t.rows.append(imbalance_array)
      html += str(t) + '</p><p>'

      t = HTML.Table(header_row=['Partition', '#', 'Countries'])
      count = 0
      partition_tuples = []
      for item in value:
         partition_tuples.append([str(count+1), item.count(','), str(item)])
         count += 1
      # sort partitions by partition size
      for item in sorted(partition_tuples, key=lambda partition: partition[1]):        
         t.rows.append(item)
      
      html += str(t)
      html += '</p><br/>Cluster to cluster imbalance contribution matrix'
      # append imbalance matrix
      html += cc_cluster_imb_matrix[key] + '<br/>'
      startyear += 1
      partial_html_cc.append(html)
   html = '<h2>CC - Frequency of countries per cluster per year</h2><p>'
   t = HTML.Table(header_row=['Year', 'Partition sizes'])
   year = 1946
   for key, value in sorted(cc_result_summary.iteritems()):   
      sizes = []
      for item in value:        
         sizes.append(item.count(','))
      pref = [year]
      for item in sorted(sizes):
         pref.append(str(item))
      t.rows.append(pref)
      year += 1
   html += str(t)
   html += '</p>'
   partial_html_cc.append(html)

   startyear = 1946
   
   for key, value in sorted(rcc_result_summary.iteritems()):
      html = '<h2>RCC - ' + str(key) + ' (' + str(startyear) + ')</h2><p>'

      imbalance_array = rcc_imbalance_summary[key]
      t = HTML.Table(header_row=['RI(P)', 'Edge sum', 'RI(P)%', 'RI(P)+', 'Pos edge sum', 'RI(P)+%', 'RI(P)-', 'Neg edge sum', 'RI(P)-%', 'RI(P)int', 'RI(P)int%', 'RI(P)ext', 'RI(P)ext%'])
      t.rows.append(imbalance_array)
      html += str(t) + '</p><p>'

      t = HTML.Table(header_row=['Partition', '#', 'Countries'])
      count = 0
      partition_tuples = []
      for item in value:
         partition_tuples.append([str(count+1), item.count(','), str(item)])
         count += 1
      
      for item in sorted(partition_tuples, key=lambda partition: partition[1]):        
         t.rows.append(item)
         
      html += str(t)
      html += '</p><br/><h4>Mediation and differential popularity analysis</h4>'
      
      html += "\n<br>Clusters with plus mediators (internal positive edges + external positive edges): <br>"
      if not key in PlusMediators:
         html += "None<br>\n"
      else:
         for elem in (PlusMediators[key]):
            html += (elem)
         html += "<br>\n"
      html += "\n<br>Clusters with plus mutually hostile mediators (internal negative edges + external positive edges): <br>"
      if not key in PlusMutuallyHostileMediators:
         html += "None<br>\n"
      else:
         for elem in (PlusMutuallyHostileMediators[key]):
            html += (elem)
         html += "<br>\n"

      html += '<br/>Cluster to cluster imbalance contribution matrix'
      # append imbalance matrix
      html += rcc_cluster_imb_matrix[key] + '<br/>'
      startyear += 1
      partial_html_rcc.append(html)
   html = '<h2>RCC - Frequency of countries per cluster per year</h2><p>'
   t = HTML.Table(header_row=['Year', 'Partition sizes'])
   year = 1946
   for key, value in sorted(rcc_result_summary.iteritems()):   
      sizes = []
      for item in value:        
         sizes.append(item.count(','))
      pref = [year]
      for item in sorted(sizes):
         pref.append(str(item))
      t.rows.append(pref)
      year += 1
   html += str(t)
   html += '</p>'
   partial_html_rcc.append(html)

   with open(folder + '/unga-cc-rcc-summary.html', 'w') as rccfile:
      rccfile.write(mhtml)
      for i in range(0, len(partial_html_cc)):
         rccfile.write(partial_html_cc[i])
         rccfile.write(partial_html_rcc[i])
      rccfile.write('</body>\n</html>')
      print 'Saved CC and RCC UNGA summary results.'
   

if __name__ == "__main__":
   main(sys.argv[1:])
