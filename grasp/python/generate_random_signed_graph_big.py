# *****************************************************
# Instance generator for random symmetric signed graphs
# INPUT:
#       n number of vertices
#       d graph density ( 2m/n(n-1) )
#       dn negative percentual (|E-|/m)
#       
#       where m=|E+|+|E-|
# OUTPUT
#       undirected graph with NO parallel edges
# OUTPUT FORMAT
#       <n> <|E+ \cup E-|>
#       <i> <j> <s>
#       .
#       .
#       .
# ****************************************************/

import sys
import csv
import StringIO
import os
import math
import argparse
import collections
import random
import time
from scipy.sparse import dok_matrix

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __str__(self):
        return "range(" + str(self.start) + ", " + str(self.end) + ")"


edge = collections.namedtuple('Edge', ['x', 'y'])

EPS = 0.0001


# Returns a random edge from the graph
def pick_random_edge(n):
    a = random.randint(0, n - 1)
    b = random.randint(0, n - 1)
    while b == a:
        b = random.randint(0, n - 1)
    return edge(a, b)


def generate_random_graph(n, d, dneg):
    M = dok_matrix((n, n))
    # [[0 for x in xrange(n)] for x in xrange(n)]

    m = 0  # number of edges
    mneg = 0  # number of neg edges
    mpos = 0  # number of pos edges

    m_min = int(math.ceil(d * (n * (n - 1) / 2.0)))
    mneg_min = int(math.ceil(dneg * (float(d) * (n * (n - 1) / 2.0))))
    edge_list = []

    previous_total = -1
    percentage = 0
    threshold = int(math.floor(mneg_min / 10.0))
    print "Generating negative edges..."
    while mneg < mneg_min:
        e = pick_random_edge(n)
        i = e.x
        j = e.y

        if M[i,j] == 0:
            mneg += 1
            m += 1
            M[i,j] = -1
            M[j,i] = -1
            edge_list.append('({0},{1})-1\n'.format(str(i + 1), str(j + 1)))
            percentage = int(math.ceil(100 * (float(mneg) / mneg_min)))
            if mneg % threshold < EPS and percentage != previous_total:
                print str(percentage) + " % ",
                previous_total = percentage
    print "\nGenerated {0} negative edges.".format(str(mneg))

    previous_total = -1
    percentage = 0
    threshold = int(math.floor(m_min - mneg / 10.0))
    m = mneg
    print "Generating positive edges..."
    while m < m_min:
        e = pick_random_edge(n)
        i = e.x
        j = e.y

        if M[i,j] == 0:
            mpos += 1
            m += 1
            M[i,j] = 1
            M[j,i] = 1
            edge_list.append('({0},{1})1\n'.format(str(i + 1), str(j + 1)))
            percentage = int(math.ceil(100 * (float(m - mneg) / (m_min - mneg))))
            if m % threshold < EPS and percentage != previous_total:
                print str(percentage) + " % ",
                previous_total = percentage
    print "\nGenerated the remaining {0} (positive) edges.".format(mpos)

    # Writes output file in XPRESS Mosel format (.mos)
    # --------------------------------------------------
    # Stores graph file in output folder
    # Vertex numbers start at 1
    filename_prefix = "file_" + str(n) + "_" + str(d) + "_" + str(dneg)
    filename = filename_prefix + ".g"
    directory = "output-random-sg"
    if not os.path.exists(directory):
        os.makedirs(directory)
    print "Saving output file..."
    with open(directory + "/" + filename, "w") as g_file:
        # file header
        g_file.write('people: {0}\r\n\r\n'.format(str(n)))
        g_file.write('VarErr: 0.5\r\n\r\n')
        vertex_names = []
        g_file.write('Names: [')
        for x in xrange(1, n + 1):
            vertex_names.append(str(x))
            vertex_names.append(" ")
            if x % 20 == 0:
                vertex_names.append('\n')
            if len(vertex_names) >= 4096:
                g_file.write(''.join(vertex_names))
                vertex_names = []
        g_file.write(''.join(vertex_names))
        g_file.write(']\r\n')
        # graph contents
        previous_total = -1
        percentage = 0
        count = 0
        list_size = len(edge_list)
        threshold = int(math.floor(list_size / 10.0))
        g_file.write('Mrel: [ \n')
        for item in edge_list:
            g_file.write(item)
            count += 1
            percentage = int(math.ceil(100 * (float(count) / list_size)))
            if count % threshold < EPS and percentage != previous_total:
                print str(percentage) + " % ",
                previous_total = percentage
        g_file.write(']\n')

    print "\nNumber of negative edges: {0} ({1}%), number of positive edges: {2} ({3}%), total edges: {4} ({5}%).".format(
        str(mneg), str(100 * mneg / m), str(mpos), str(100 * mpos / m), str(m), str(100 * m / (n * (n - 1) / 2)))
    print "Output file successfully generated."


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, required=True, help="number of nodes")
    parser.add_argument('-d', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="density of edges of the graph")
    parser.add_argument('-nd', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="negative density of edges of the graph")
    args = parser.parse_args()

    # measure elapsed time during graph generation
    start = time.time()

    print "Will generate a completely random graph of size {0}, density = {1} and negative density = {2}".format(args.n, args.d, args.nd)
    # Call random graph generation procedure
    generate_random_graph(args.n, args.d, args.nd)

    end = time.time()
    elapsed = end - start
    print "Graph generation took {0:.2f} seconds.".format(elapsed)


if __name__ == "__main__":
    main(sys.argv[1:])
