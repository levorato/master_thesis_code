# script para processamento do resultado das instancias House of Cunha (Congresso Nacional)
# associa os numeros dos vertices aos nomes dos deputados / partidos / estados

# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

# .dcode files location is the same folder where the instance files are located

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import math
import pandas as pd

import HTML


# gambit para python 2.4
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
    print [tryint(c) for c in re.split('([0-9]+)', s)][1]
    return int([tryint(c) for c in re.split('([0-9]+)', s)][1])


def natsorted(l):
    """ Sort the given list in the way that humans expect.
    """
    temp = [(alphanum_key(x), x) for x in l]
    temp.sort()


#def process_result_for_year(year):



def main(argv):

    csv.field_size_limit(1000000000)

    folder = ''
    filter = ''
    instances_path = ''
    threshold = float(90.0)  # percentual threshold for the quantity of relashionships
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["folder=", "filter=", "instancespath="])
    except getopt.GetoptError:
        print 'process_congress_graphs.py --folder <folder> --filter <filter> --instancespath <instances_path>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'process_congress_graphs.py --folder <folder> --filter <filter> --instancespath <instances_path>'
            sys.exit()
        elif opt in ("-i", "--folder"):
            folder = arg
        if opt in ("-f", "--filter"):
            filter = arg
        if opt in ("-p", "--instancespath"):
            instances_path = arg

    if (folder == ''):
        print 'Please specify the results dir'
        sys.exit()
    if (instances_path == ''):
        print 'Please specify the graph instances path'
        sys.exit()
    print 'Input dir is ', folder
    print 'File filter is ', filter
    print 'Graph instances dir is ', instances_path
    print 'Mediation threshold is ', threshold

 #   for year in [2010, 2011, 2012, 2013, 2014, 2015, 2016]:
 #       process_result_for_year(year)

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
    all_cc_results_df = pd.DataFrame({'A': []})
    all_srcc_results_df = pd.DataFrame({'A': []})

    previous_filename = ''
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()
        print "Processing folder " + ''.join(root)
        if (len(files)):
            file_list = []
            file_list.extend(glob.glob(root + "/*-result.txt"))
            file_count = len(file_list) - 1

            while file_count >= 0:
                print "Processing file " + file_list[file_count] + "\n"
                content_file = open(file_list[file_count], 'r')
                prefix = file_list[file_count]
                result_file_path = prefix[0:prefix.rfind('.txt')]
                result_file = open(result_file_path + '-addinfo.txt', "w")
                result_file.write("Additional results info\n")
                result_file_name = result_file_path[result_file_path.rfind(os.path.sep) + 1:]

                # lookup deputy name, party and state from csv file 'congress-<year>-dcode.csv'
                country_full_name = dict()
                graphfile = prefix[0:prefix.rfind(os.path.sep)]
                graphfile = graphfile[0:graphfile.rfind(os.path.sep)]
                # graphfile is the graph file currently being processed
                graphfile = graphfile[graphfile.rfind(os.path.sep) + 1:]
                file_prefix = graphfile[0:graphfile.rfind('-v')]
                graph_version = graphfile[graphfile.rfind('-')+1:graphfile.rfind('.')]
                dcode_filename = os.path.join(instances_path, file_prefix + '-dcode.csv')
                print file_prefix
                current_year = file_prefix[file_prefix.find('-')+1:]
                print current_year
                dcode_df = pd.read_csv(dcode_filename, encoding="utf-8-sig", sep=';', header=0,
                                       names=['id_dep', 'id_vertex', 'name', 'party', 'state'])

                if result_file_name <> previous_filename:
                    print 'Reading graph file: {0}'.format(str(instances_path + os.path.sep + graphfile))
                    with open(instances_path + os.path.sep + graphfile, 'rb') as csvfile:
                        dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=" \t")
                        csvfile.seek(0)
                        r = csv.reader(csvfile, dialect)
                        line1 = r.next()
                        values = [col for col in line1]
                        header_str = values[0]
                        if header_str.find('\t') > 0:
                            header_str = header_str[0:header_str.find('\t')]
                        if header_str.find(' ') > 0:
                            N = int(header_str[:header_str.find(' ')])
                        else:
                            N = int(header_str)
                            # print 'N is {0}'.format(str(N))
                        matrix = [[0 for x in xrange(N)] for x in xrange(N)]
                        pos_edge_sum = 0.0
                        neg_edge_sum = 0.0
                        for line in r:
                            i = int(line[0])
                            j = int(line[1])
                            w = float(line[2])
                            # congress graphs are undirected
                            matrix[i][j] = w
                            matrix[j][i] = w
                            # print '({0}, {1}) = {2}\n'.format(str(i), str(j), str(w))
                            if (w < 0):
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
                        vertex_contrib = [[float(0.0), float(0.0)] for x in xrange(N)]
                        pos_contrib_sum = float(0.0)
                        neg_contrib_sum = float(0.0)
                        started_cluster_analysis = False
                        processed_imb_analysis = False
                        started_imb_analysis = False
                        partitionNumber = 0
                        cluster = [int(0) for x in xrange(N)]
                        vertex_cluster_tuples = []
                        vertices_in_cluster = []
                        for row in reader:
                            linestring = ''.join(row)
                            if (linestring == ''):
                                continue
                            if (linestring.startswith(' ')):
                                processed = True
                                vertex_list = linestring[linestring.find('[') + 1:linestring.rfind(']') - 1].strip().split()
                                vertices_in_cluster.append(vertex_list)
                                #print 'List of vertices in cluster {0}: {1}'.format(partitionNumber, vertex_list)
                                for vertex in vertex_list:
                                    cluster[int(vertex)] = partitionNumber
                                    vertex_cluster_tuples.append((int(vertex), partitionNumber))
                                partitionNumber += 1
                            else:
                                if (not processed):
                                    result_file.write(str(row[0]) + '\n')
                                else:
                                    # processes imbalance analysis data (out edges contribution)
                                    if (linestring.startswith('Cluster')):
                                        if (not started_cluster_analysis):
                                            started_cluster_analysis = True
                                        # initialize matrix with cluster contribution
                                        numberOfClusters = partitionNumber
                                        clusterImbMatrix = [[float(0.0) for x in xrange(numberOfClusters)] for x in
                                                            xrange(numberOfClusters)]
                                        line = 0
                                    else:
                                        if (not processed_imb_analysis):
                                            if (linestring.startswith('Vertex')):
                                                continue
                                            if (linestring.startswith('Imbalance')):
                                                if (not started_imb_analysis):
                                                    started_imb_analysis = True
                                                else:
                                                    processed_imb_analysis = True
                                                continue
                                            # print 'line: {0}'.format(str(linestring))
                                            vertex = int(linestring[:linestring.find(',')])
                                            pos_contrib = float(
                                                linestring[linestring.find(',') + 1:linestring.rfind(',')])
                                            neg_contrib = float(linestring[linestring.rfind(',') + 1:])
                                            vertex_contrib[vertex][0] = pos_contrib
                                            vertex_contrib[vertex][1] = neg_contrib
                                            pos_contrib_sum += pos_contrib
                                            neg_contrib_sum += neg_contrib
                                            # print '{0} {1}'.format(str(neg_contrib_sum), str(neg_contrib))
                                        if (started_cluster_analysis):
                                            # reads matrix with cluster contribution to imbalance
                                            reader3 = csv.reader(StringIO.StringIO(linestring), delimiter=',')
                                            column = 0
                                            for elem in reader3.next():
                                                if (str(elem).strip() != ''):
                                                    clusterImbMatrix[line][column] = float(str(elem))
                                                    column += 1
                                            line += 1

                        cluster_df = pd.DataFrame(vertex_cluster_tuples, columns=['id_vertex', 'id_cluster'])
                        result_df = pd.merge(dcode_df, cluster_df, on='id_vertex')
                        result_df['year'] = current_year
                        result_df['version'] = graph_version
                        result_df = result_df.set_index(['id_vertex', 'year', 'version'])
                        print len(result_df.index)


                        # from this point on, makes calculations based on the graph (imbalance)
                        if (result_file_name.startswith('rcc')):  # process mediation info from rcc result
                            # for each partition, tries to find a partition where most of the external edges are positive
                            # print "Cluster,%IntPosEdges,%IntNegEdges,%ExtPosEdges,%ExtNegEdges"
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
                                    PIntPosEdges = float(numberOfIntPosEdges) * 100 / totalNumberOfIntEdges
                                    PIntNegEdges = float(numberOfIntNegEdges) * 100 / totalNumberOfIntEdges
                                else:
                                    PIntPosEdges = PIntNegEdges = 0
                                if totalNumberOfExtEdges > 0:
                                    PExtPosEdges = float(numberOfExtPosEdges) * 100 / totalNumberOfExtEdges
                                    PExtNegEdges = float(numberOfExtNegEdges) * 100 / totalNumberOfExtEdges
                                else:
                                    PExtPosEdges = PExtNegEdges = 0
                                    # print str(c) + ",%.2f,%.2f,%.2f,%.2f" % (PIntPosEdges, PIntNegEdges, PExtPosEdges, PExtNegEdges)

                                # internal pos + external pos : "plus mediators" fig 2, according to Doreian et. al
                                # internal neg + external pos : "plus mutually hostile mediators" fig 3
                                # maybe internal neg + external neg would be "internal subgroup hostility" ???
                                if (PIntPosEdges > threshold and PExtPosEdges > threshold):
                                    PlusMediators[graphfile] = (
                                    "Cluster " + str(c + 1) + str(" (%%IntPosEdges = %.2f" % (PIntPosEdges)) + str(
                                        " and %%ExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
                                if (PIntNegEdges > threshold and PExtPosEdges > threshold):
                                    PlusMutuallyHostileMediators[graphfile] = (
                                    "Cluster " + str(c + 1) + str(" (%%IntNegEdges = %.2f" % (PIntNegEdges)) + str(
                                        " and %%ExtPosEdges = %.2f" % (PExtPosEdges)) + ")")
                                if (PIntNegEdges > threshold and PExtNegEdges > threshold):
                                    InternalSubgroupHostility[graphfile] = (
                                    "Cluster " + str(c + 1) + str(" (%%IntNegEdges = %.2f" % (PIntNegEdges)) + str(
                                        " and %%ExtNegEdges = %.2f" % (PExtNegEdges)) + ")")

                        # export cluster imb matrix to html
                        matrixline = ['Cluster']
                        for line in xrange(1, numberOfClusters + 1):
                            matrixline.append('%d' % line)
                        matrixline.append('Sum')
                        t = HTML.Table(header_row=matrixline)
                        for line in xrange(0, numberOfClusters):
                            matrixline = ['<b>' + str(line + 1) + '</b>']
                            sum = float(0.0)
                            for column in xrange(0, numberOfClusters):
                                if (clusterImbMatrix[line][column] > 0):
                                    matrixline.append('%.4f' % clusterImbMatrix[line][column])
                                else:
                                    matrixline.append(
                                        '<font color=\"red\">%.4f</font>' % clusterImbMatrix[line][column])
                                sum += abs(clusterImbMatrix[line][column])
                            matrixline.append('%.4f' % sum)
                            t.rows.append(matrixline)
                        if (result_file_name.startswith('cc')):
                            cc_cluster_imb_matrix[graphfile] = str(t)
                        else:
                            rcc_cluster_imb_matrix[graphfile] = str(t)
                        # end process *-result.txt files


                        # process *-Node0-*-iterations.csv file
                        csv_file_list = []
                        if (result_file_name.startswith('cc')):
                            prfx = 'CC'
                        else:
                            prfx = 'RCC'
                            pos_imbalance = pos_contrib_sum
                            neg_imbalance = neg_contrib_sum
                        csv_file_list.extend(glob.glob(root + os.path.sep + prfx + '-Node0-*-iterations.csv'))
                        if not csv_file_list:  # if there is no file, I(P) is probably zero
                            with open(root + os.path.sep + result_file_name + '.txt', 'r') as r_csv_file:
                                r2 = csv.reader(r_csv_file)
                                line = r2.next()
                                while (line[0].find('I(P) = ') < 0):
                                    line = r2.next()
                                # extracts total, positive and negative imbalance of the best result
                                imbalance = float(line[0][line[0].find('=') + 1:])
                                if (result_file_name.startswith('cc')):
                                    pos_imbalance = imbalance
                                    neg_imbalance = 0
                                    print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance),
                                                                                str(neg_imbalance))
                                else:  # for RCC problem the numbers are internal and external imbalance values
                                    int_imbalance = imbalance
                                    ext_imbalance = 0
                                    print 'imbalances {0} (int){1} (ext){2}'.format(str(imbalance), str(int_imbalance),
                                                                                    str(ext_imbalance))
                                    print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance),
                                                                                str(neg_imbalance))
                        else:
                            with open(csv_file_list[0], 'r') as r_csv_file:
                                r2 = csv.reader(r_csv_file)
                                line = r2.next()
                                while (not line[0].startswith('Best')):
                                    line = r2.next()
                                # extracts total, positive and negative imbalance of the best result
                                imbalance = float(line[1])
                                if (result_file_name.startswith('cc')):
                                    pos_imbalance = float(line[2])
                                    neg_imbalance = float(line[3])
                                    print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance),
                                                                                str(neg_imbalance))
                                else:  # for RCC problem the numbers are internal and external imbalance values
                                    int_imbalance = float(line[2])
                                    ext_imbalance = float(line[3])
                                    print 'imbalances {0} (int){1} (ext){2}'.format(str(imbalance), str(int_imbalance),
                                                                                    str(ext_imbalance))
                                    print 'imbalances {0} (+){1} (-){2}'.format(str(imbalance), str(pos_imbalance),
                                                                                str(neg_imbalance))

                        # for line in clustering_full_names:
                        #     result_file.write(line)
                        # result_file.write('\r\n')
                        # for line in clustering_names:
                        #     result_file.write(line)
                        # result_file.write('\r\n')
                        # for line in clustering_numbers:
                        #     result_file.write(line)

                        # Codigo de criacao das tabelas de deputados por ano
                        deputies_in_cluster = []
                        result_df['party'] = result_df['party'].map(lambda x: x.strip())
                        result_df['state'] = result_df['state'].map(lambda x: x.strip())
                        result_df['name'] = result_df['name'].map(lambda x: x.strip())
                        result_df['nome_e_partido'] = result_df['name'] + ' (' + result_df['party']  + '-' + result_df['state'] + ')'
                        deputies_in_cluster_series = result_df[['id_cluster', 'nome_e_partido']].groupby(['id_cluster'])['nome_e_partido'].apply(lambda x: u', '.join(x).encode('utf-8').strip())
                        deputies_in_cluster_df = deputies_in_cluster_series.to_frame()
                        deputies_in_cluster_df['#'] = result_df[['id_cluster', 'nome_e_partido']].groupby(['id_cluster'])['nome_e_partido'].apply(lambda x: len(x))
                        #print deputies_in_cluster_df
                        # agrega todos os resultados (por ano e por versao do grafo) em um dataframe para o CC e outro para o SRCC
                        if (result_file_name.startswith('cc')):
                            if all_cc_results_df.empty:
                                all_cc_results_df = result_df
                            else:
                                all_cc_results_df = all_cc_results_df.append(result_df)
                            print len(all_cc_results_df.index)
                        else:
                            if all_srcc_results_df.empty:
                                all_srcc_results_df = result_df
                            else:
                                all_srcc_results_df = all_srcc_results_df.append(result_df)
                            print len(all_srcc_results_df.index)
                            # print all_cc_results_df

                        if (result_file_name.startswith('cc')):
                            pd.set_option('display.max_colwidth', -1)
                            cc_result_summary[graphfile] = deputies_in_cluster_df.to_html(columns = [1, 0],
                                                                                          formatters = {'nome_e_partido': (lambda x: u', '.join(x).encode('utf-8').strip())})
                            cc_imbalance_summary[graphfile] = [str(imbalance), str(pos_edge_sum + neg_edge_sum), str(
                                round(100 * imbalance / (pos_edge_sum + neg_edge_sum), 2)), str(pos_imbalance),
                                                               str(pos_edge_sum),
                                                               str(round(100 * pos_imbalance / pos_edge_sum, 2)),
                                                               str(neg_imbalance), str(neg_edge_sum),
                                                               str(round(100 * neg_imbalance / neg_edge_sum, 2))]
                        else:  # rcc
                            rcc_result_summary[graphfile] = deputies_in_cluster_df.to_html()
                            rcc_imbalance_summary[graphfile] = [str(imbalance), str(pos_edge_sum + neg_edge_sum), str(
                                round(100 * imbalance / (pos_edge_sum + neg_edge_sum), 2)), str(pos_imbalance),
                                                                str(pos_edge_sum),
                                                                str(round(100 * pos_imbalance / pos_edge_sum, 2)),
                                                                str(neg_imbalance), str(neg_edge_sum),
                                                                str(round(100 * neg_imbalance / neg_edge_sum, 2)),
                                                                str(int_imbalance), str(
                                    round(100 * int_imbalance / (pos_edge_sum + neg_edge_sum), 2)), str(ext_imbalance),
                                                                str(round(
                                                                    100 * ext_imbalance / (pos_edge_sum + neg_edge_sum),
                                                                    2))]
                    finally:
                        content_file.close()
                        result_file.close()
                    file_count = file_count - 1
                    previous_filename = graphfile
                # end process results  # TODO implementar tabela com os valores de contribuicao pos/neg entre clusters

    print "\nSuccessfully processed all graph files.\n"

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(os.path.join(folder, 'congress-all_results.xlsx'), engine='xlsxwriter')

    all_cc_results_df['party'] = all_cc_results_df['party'].map(lambda x: x.strip())
    all_cc_results_df['state'] = all_cc_results_df['state'].map(lambda x: x.strip())
    all_cc_results_df['name'] = all_cc_results_df['name'].map(lambda x: x.strip())

    # Convert the dataframe to an XlsxWriter Excel object.
    all_cc_results_df.to_excel(writer, sheet_name='CC')
    all_srcc_results_df.to_excel(writer, sheet_name='SRCC')

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
    print "\nExported all results to Excel file.\n"

    # exports summary to HTML file
    # Exporting CC and RCC results
    mhtml = '<!DOCTYPE html>\n<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://genshi.edgewall.org/">\n<head><title>House of Cunha CC and RCC Summary</title>'
    mhtml += '</head><body class="index">'
    mhtml += '<h1>Brazilian House of Deputies Voting Data - CC and RCC Imbalance Analysis</h1>'
    startyear = 2010
    partial_html_cc = []
    partial_html_rcc = []

    for key, value in sorted(cc_result_summary.iteritems()):
        html = '<h2>CC - ' + str(key) + '</h2><p>'

        imbalance_array = cc_imbalance_summary[key]
        t = HTML.Table(header_row=['I(P)', 'Edge sum', 'I(P)%', 'I(P)+/I(P)ext', 'Pos edge sum', 'I(P)+%', 'I(P)-/I(P)int',
                                   'Neg edge sum', 'I(P)-%'])
        t.rows.append(imbalance_array)
        html += str(t) + '</p><p>'

        # t = HTML.Table(header_row=['Partition', '#', 'Deputies'])
        # count = 0
        # partition_tuples = []
        # for item in value:
        #     partition_tuples.append([str(count + 1), len(item), str(item)])
        #     count += 1
        # # sort partitions by partition size
        # for item in sorted(partition_tuples, key=lambda partition: partition[1]):
        #     t.rows.append(item)

        html += u''.join(value).encode('utf-8').strip()
        html += '</p><br/>Cluster to cluster imbalance contribution matrix'
        # append imbalance matrix
        html += cc_cluster_imb_matrix[key] + '<br/>'
        startyear += 1
        partial_html_cc.append(html)
    html = '<h2>CC - Frequency of deputies per cluster per year</h2><p>'
    t = HTML.Table(header_row=['Year', 'Partition sizes'])
    year = 2010
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

    startyear = 2010

    for key, value in sorted(rcc_result_summary.iteritems()):
        html = '<h2>SRCC - ' + str(key) + '</h2><p>'

        imbalance_array = rcc_imbalance_summary[key]
        t = HTML.Table(
            header_row=['RI(P)', 'Edge sum', 'RI(P)%', 'RI(P)+', 'Pos edge sum', 'RI(P)+%', 'RI(P)-', 'Neg edge sum',
                        'RI(P)-%', 'RI(P)int', 'RI(P)int%', 'RI(P)ext', 'RI(P)ext%'])
        t.rows.append(imbalance_array)
        html += str(t) + '</p><p>'

        # t = HTML.Table(header_row=['Partition', '#', 'Deputies'])
        # count = 0
        # partition_tuples = []
        # for item in value:
        #     partition_tuples.append([str(count + 1), len(item), str(item)])
        #     count += 1
        #
        # for item in sorted(partition_tuples, key=lambda partition: partition[1]):
        #     t.rows.append(item)

        #html += str(t)
        html += u''.join(value).encode('utf-8').strip()
        html += '</p><br/><h4>Mediation analysis</h4>'

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
    year = 2010
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

    with open(os.path.join(folder, 'congress-cc-rcc-summary.html'), 'w') as rccfile:
        rccfile.write(mhtml)
        for i in range(0, len(partial_html_cc)):
            rccfile.write(partial_html_cc[i])
            rccfile.write(partial_html_rcc[i])
        rccfile.write('</body>\n</html>')
        print 'Saved CC and RCC Congress summary results.'

if __name__ == "__main__":
    main(sys.argv[1:])
