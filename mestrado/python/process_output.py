import sys, getopt
import csv
import StringIO
import glob
import os
import os.path

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
   best_file_summary = dict()
   avg_file_summary = dict()
   timeInterval = dict()
   timeCount = dict()
   previous_filename = ""
   avg_ip_const = 0
   avg_value = 0
   avg_k = 0
   avg_time = 0
   avg_count = 0
   avg_iter = 0

   for root, subFolders, files in os.walk(folder):
      # sort dirs and files
      subFolders.sort()
      files.sort()      

      print "Processing folder " + ''.join(root)
      if(len(files) and ''.join(root) != folder):
         file_list = []
         file_list.extend(glob.glob(root + "/" + filter))
         count = len(file_list) - 1
         
         if os.path.isfile(root + "/result.txt"):
		 input_file = open(root + "/result.txt", "r")
		 
		 content = input_file.read()
		 reader = csv.reader(StringIO.StringIO(content), delimiter='\n')
		 for row in reader:
		    if ''.join(row).find("time spent:") >= 0:
		       line = ''.join(row)
		       global_time = float(line[line.find("time spent:")+11:])
		       break
		 input_file.close

		 text_file = open(root + "/summary.txt", "w")
		 filename = (root[:root.rfind("/")])
		 datetime = root[root.rfind("/")+1:]
		 filename = filename[filename.rfind("/")+1:]
		 text_file.write("Summary for graph file: %s\n"%filename)
		 local_avg_ip_const = 0
		 local_avg_count = 0
		 best_value = 1000000L
		 best_pos_value = 0
		 best_neg_value = 0
		 best_K = 0
		 best_iteration = 0
		 best_time = 0
		 best_param = ''
		 
		 while count >= 0:
		    
		   content_file = open(file_list[count], 'r')
		   try:
		    content = content_file.read()

		    reader = csv.reader(StringIO.StringIO(content), delimiter=',')
		    for row in reader:
		       linestring = ''.join(row)
		       column = []
                       for col in row:
                          column.append(col)
		       if linestring.startswith( 'Best value' ):          
		          filepath = ''.join(file_list[count])
		          text_file.write(filepath[filepath.rfind("/")+1:] + ' ' + linestring + '\n')
		          value = float(column[1])
		          pos_value = float(column[2])
		          neg_value = float(column[3])
		          K = long(column[4])
		          iteration = long(column[5]) 
		          time = float(column[6])
		          total_iter = long(column[7])
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
		       elif linestring.startswith( 'Average initial I(P)' ):
			     # computa o valor medio da solucao inicial da fase construtiva para determinado no da execucao paralela
			     ip_const = float(column[1])
			     local_avg_ip_const = local_avg_ip_const + ip_const 
			     local_avg_count = local_avg_count + 1
		    count = count - 1
		   finally:
		    content_file.close()
		 text_file.close()
		 all_files_summary[filename+"/"+datetime] = str(best_value)+", "+str(pos_value)+", "+str(neg_value)+", "+str(best_K)+", "+str(iteration)+", "+str(best_time)+", "+str(global_time)+", "+best_param+", "+str(total_iter)
		 # calcula a media dos valores de I(P) da fase de construcao para este set de arquivos de resultado
		 local_avg_ip_const = local_avg_ip_const / local_avg_count
		 # armazena os valores de todas as execucoes de um mesmo grafo para calculo da media
		 if filename == previous_filename:
		    avg_ip_const = avg_ip_const + local_avg_ip_const
		    avg_value = avg_value + best_value
		    avg_k = avg_k + best_K
		    avg_time = avg_time + global_time
		    avg_iter = avg_iter + total_iter
		    avg_count = avg_count + 1
		 else:
		    if avg_count > 0:
		        print "storing " + previous_filename
		        avg_file_summary[previous_filename] = str(avg_ip_const / avg_count) + ", " + str(avg_value / avg_count)+", "+str(avg_k / avg_count)+", "+str(avg_time / avg_count)+", "+str(avg_iter / avg_count)+", "+str(avg_count)
		        print "average execution times for file " + previous_filename
			tdir = "./times"
                        if not os.path.exists(tdir):
				os.makedirs(tdir)
			times_file = open(tdir + "/" + previous_filename + "-executionTimes.txt", "w")
			for key, value in sorted(timeInterval.items()):
      				times_file.write(str(key) + "," + str(value / timeCount[key]) + "\n")
			times_file.close()
			timeInterval = dict()
			timeCount = dict()
		    avg_ip_const = local_avg_ip_const 
		    avg_value = best_value
		    avg_k = best_K
		    avg_time = global_time
		    avg_iter = total_iter
		    avg_count = 1      

		 # captura o melhor resultado dadas todas as execucoes de um mesmo grafo
		 if best_file_summary.has_key(filename):
		    element = best_file_summary[filename]
		    value = float(element[0:element.find(',')-1])
		    if(best_value < value):
		       best_file_summary[filename] = str(all_files_summary[filename+"/"+datetime])
		 else:
		    best_file_summary[filename] = str(all_files_summary[filename+"/"+datetime])

		 previous_filename = filename
	 # varre os arquivos da pasta em busca dos intervalos de tempo	
	 for filename in files:
		filename = os.path.join(root, filename)
		if "timeIntervals" in filename:
			time_list = dict()
			input_file_t = open(filename, "r")
		 
			content_t = input_file_t.read()
			reader_t = csv.reader(StringIO.StringIO(content_t), delimiter=',')

			last_time = 0.0
			for row in reader_t:
				if len(row) > 3:
					time_list[float(row[0])] = row
					last_time = float(row[0])
			input_file_t.close
			
			# print "Time intervals for file " + filename
			# print "Last time is " + str(last_time)
			time = 0.0
			time_list_keys = time_list.keys()
			time_list_keys.sort()
			# print "Computed times are " + str(time_list_keys)
			# print "Time slot(s), I(P), Details (Time(s), I(P), I(P)+, I(P)-, k, Iter)"
			graph_name = filename[:filename.rfind("/")]
			graph_name = graph_name[:graph_name.rfind("/")]
			graph_name = graph_name[graph_name.rfind("/")+1:]
			print "graph name is " + graph_name  
			while time <= last_time:
				index = min(time_list_keys, key=lambda x: q(x < time, time - x, x))
				timeInterval[time] = float(timeInterval.get(time, 0.0)) + float(time_list[index][1])
				timeCount[time] = timeCount.get(time, 0.0) + 1
				# str(time) + ", " + str(time_list[index][1]) + ", " + str(time_list[index])
				time += 10.0
			timeInterval[time] = float(timeInterval.get(time, 0.0)) + float(time_list[last_time][1])
			timeCount[time] = timeCount.get(time, 0.0) + 1
			# str(last_time) + ", " + str(time_list[last_time][1]) + ", " + str(time_list[last_time])

   if avg_count > 0:
      print "storing " + previous_filename
      avg_file_summary[previous_filename] = str(avg_ip_const / avg_count)+", "+str(avg_value / avg_count)+", "+str(avg_k / avg_count)+", "+str(avg_time / avg_count)+", "+str(avg_iter / avg_count)+", "+str(avg_count)
      print "average execution times for file " + previous_filename
      tdir = "./times"
      if not os.path.exists(tdir):
	      os.makedirs(tdir)
      times_file = open(tdir + "/" + previous_filename + "-executionTimes.txt", "w")
      for key, value in sorted(timeInterval.items()):
    	times_file.write(str(key) + "," + str(value / timeCount[key]) + "\n")
      times_file.close()
      timeInterval = dict()
      timeCount = dict()
   
   result_file = open(folder + "/summary.txt", "w")
   print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter"
   result_file.write("Filename, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter\n")
   for key in sorted(all_files_summary.iterkeys()):
      print "%s, %s" % (key, all_files_summary[key])
      result_file.write("%s, %s\n" % (key, all_files_summary[key]))
   result_file.close()
   print "------ Best results:"
   for key in sorted(best_file_summary.iterkeys()):
      print "%s, %s" % (key, best_file_summary[key])
   print "------ Average results:"
   print "Instance, Avg I(P) const, Avg I(P), Avg K, Avg Time(s), Avg Iter, Num executions"
   for key in sorted(avg_file_summary.iterkeys()):
      print "%s, %s" % (key, avg_file_summary[key])


if __name__ == "__main__":
   main(sys.argv[1:])
