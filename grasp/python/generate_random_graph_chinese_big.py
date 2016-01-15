# ========================================================================
# Random graph generation script, according to the article
# 'Community Mining from Signed Social Networks'
# (URL: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4302742)
# Author: Mario Costa Levorato Junior
# ========================================================================

import sys
import csv
import StringIO
import glob
import os
import math
from array import *
import random
import argparse
import collections
from random import randint
from time import sleep
import time
from scipy.sparse import dok_matrix

# Global variables
mycluster = []
cluster_node_list = []
vertex_list = []
matrix = dok_matrix  # original directed graph from file


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


# Returns a random node of in-degree + out-degree < k
def pick_random_vertex(v_list, indegree, outdegree, k):
    if len(v_list) == 0:
        # print 'Empty vertex list!!!\n'
        return -1
    a = random.choice(v_list)
    v_list.remove(a)
    while (indegree[a] + outdegree[a]) >= k:
        if len(v_list) == 0:
            # print 'Empty vertex list!!!\n'
            return -1
        a = random.choice(v_list)
        v_list.remove(a)
    return a


# Returns a new random edge that connects 2 nodes of the same cluster
# Each node must have in-degree + out-degree < k ???
def pick_random_internal_edge(N, indegree, outdegree, k):
    # Global variables
    global mycluster
    global cluster_node_list
    global matrix

    vertex_list = range(N)
    random.shuffle(vertex_list)

    found = False
    while not found:
        # Selects a random node of in-degree + out-degree < k
        a = pick_random_vertex(vertex_list, indegree, outdegree, k)
        if a < 0:
            return None
        # print "1- Selected vertex %s" % str(a)
        # Selects a neighbour of a (= same cluster) with degree < k
        # print 'mycluster of {0} is {1}'.format(a, mycluster[a])
        neighbour_list = cluster_node_list[mycluster[a]]
        random.shuffle(neighbour_list)
        for b in neighbour_list:
            if b != a and matrix[a, b] == 0 and indegree[b] + outdegree[b] < k:
                found = True
                break

    return edge(a, b)


# Returns a new random edge that connects 2 nodes of different clusters
# Each node must have degree < k
def pick_random_external_edge(N, indegree, outdegree, k, c):
    # Global variables
    global mycluster
    global matrix
    global vertex_list
    global cluster_node_list

    cluster_list = range(c)
    random.shuffle(cluster_list)

    # removes the vertices that have exceeded their maximum degrees
    for a in vertex_list:
        if indegree[a] + outdegree[a] >= k:
            vertex_list.remove(a)
            # print 'removed node {0}\n'.format(str(a))
    random.shuffle(vertex_list)

    found = False
    i = 0
    while i < len(vertex_list):
        # Selects a random node of degree < k which belongs to vertex_list
        if len(vertex_list) == 0:
            return None
        a = vertex_list[i]
        i += 1
        # print 'selected node {0}\n'.format(str(a))
        # Selects a vertex b, non-neighbour of a (different cluster) with degree < k
        # vertex b belongs to cluster c1 in cluster_list
        for c1 in cluster_list:
            if c1 != mycluster[a]:
                random.shuffle(cluster_node_list[c1])
                for b in cluster_node_list[c1]:
                    if matrix[a, b] == 0 and (indegree[b] + outdegree[b] < k):
                        found = True
                        break
            if found:
                break
        if found:
            break
    # end while
    if not found:
        print 'b not found!'
        return None
    else:
        return edge(a, b)


def CCObjectiveFunction(N):
    # Global variables
    global mycluster
    global matrix
    positiveSum = float(0)
    negativeSum = float(0)

    # For each vertex i
    for i in xrange(N):
        # Find out to which cluster vertex i belongs
        k = mycluster[i]
        # For each out edge of i
        for j in xrange(N):
            if matrix[i, j] != 0:
                weight = float(matrix[i, j])
                sameCluster = (k == mycluster[j])
                if weight < 0 and sameCluster:  # negative edge
                    # i and j are in the same cluster
                    negativeSum += weight * (-1)
                elif weight > 0 and not sameCluster:  # positive edge
                    # i and j are NOT in the same cluster
                    positiveSum += weight

    print 'CC calculated value. Obj = {0}\n'.format(str(positiveSum + negativeSum))
    return positiveSum + negativeSum


def SG(c, n, k, p_in, p_minus, p_plus):
    # Global variables
    global mycluster
    global cluster_node_list
    global matrix
    global vertex_list
    # Variables used in graph generation
    N = n * c
    print 'Generating random graph with {0} vertices...'.format(str(N))
    success = False

    while not success:

        # Array that controls to which cluster a vertex belongs (size = N)
        mycluster = [0 for x in xrange(N)]
        # List of nodes inside a cluster
        cluster_node_list = [[] for x in xrange(c)]
        # Graph adjacency matrix (with dimensions N x N)
        # Creates a list containing N lists initialized to 0 (index begins with 0)
        matrix = dok_matrix((N, N))
        # Arrays that control the degree of each vertex
        indegree = [0 for x in xrange(N)]
        outdegree = [0 for x in xrange(N)]

        # builds a list with all the N vertices of the graph
        vertex_list = range(N)

        # 1- randomly chooses n vertices to insert in each of the c clusters
        print "Step 1: randomly chooses n vertices to insert in each of the c clusters..."
        for cluster_number in xrange(c):
            for i in xrange(n):
                # randomly removes a vertex from the vertex_list
                v = random.choice(vertex_list)
                vertex_list.remove(v)
                mycluster[v] = cluster_number
                cluster_node_list[cluster_number].append(v)
        assert len(vertex_list) == 0, "All vertices must belong to a cluster"

        # assertion tests / sanity check
        cluster_count = [0 for x in xrange(c)]
        for v in xrange(N):
            cluster_count[mycluster[v]] += 1
        for x in xrange(c):
            assert cluster_count[x] == n, "There must be n vertices in each cluster"
            assert len(cluster_node_list[x]) == n, "There must be n vertices in each cluster"
        assert len(vertex_list) == 0, "Vertex list must be empty after cluster-vertex distribution"

        # 2- uses probability p_in to generate k edges (initially positive) connecting each of the vertices
        # p_in % of the edges connect nodes of the same cluster
        print "Step 2: uses probability p_in to generate k edges (initially positive) connecting each of the vertices..."
        total_edges = int(math.ceil((c * n * k) / 2))
        internal_edge_num = int(math.ceil(total_edges * p_in))
        external_edge_num = total_edges - internal_edge_num
        print 'Number of total edges of the graph is {0}'.format(str(total_edges))
        print 'Number of expected internal edges of the graph is {0}'.format(str(internal_edge_num))
        print 'Number of expected external edges of the graph is {0}'.format(str(external_edge_num))
        assert (internal_edge_num + external_edge_num == total_edges), "Sum of internal and external edges do not match"
        vertex_list = range(N)
        random.shuffle(vertex_list)

        # Guarantees that each node has exactly k edges (TODO: in-degree + out-degree = k ???)
        # 2.1- Creates the internal edges (within same cluster), initially positive
        print "Generating internal edges..."
        # count = 0
        # for i in xrange(internal_edge_num):
        # Randomly picks a pair of nodes that belong to the same cluster
        # e = pick_random_internal_edge(N, indegree, outdegree, k)
        # if e == None:
        #    print "Warn: no more internal edge found."
        #    break
        # edges are directed
        # matrix[e.x][e.y] = 1
        # outdegree[e.x] += 1
        # indegree[e.y] += 1
        # count += 1
        # for each node, generate an internal cluster edge
        int_count = 0
        # each cluster will get (internal_edge_num) / c internal edges
        edge_count = int(math.floor((internal_edge_num) / c))
        rest = internal_edge_num - (c * edge_count)
        while int_count < c * edge_count:
            inserted = False
            print 'there are {0} internal edges so far'.format(int_count)
            for cluster in xrange(c):
                print 'must generate {0} internal edges for cluster {1}'.format(str(edge_count), str(cluster))
                i = 0
                previous_total = -1
                percentage = 0
                while i < edge_count:
                    random.shuffle(cluster_node_list[cluster])
                    for node in cluster_node_list[cluster]:
                        if i >= edge_count:
                            break
                        if indegree[node] + outdegree[node] < k:
                            # selects another node in the same cluster
                            # print 'selected node {0} with outdegree = {1}\n'.format(str(node), str(outdegree[node]))
                            neighbour_list = cluster_node_list[mycluster[node]]
                            random.shuffle(neighbour_list)
                            for b in neighbour_list:
                                if indegree[node] + outdegree[node] >= k:
                                    break
                                if b != node and matrix[node, b] == 0 and (indegree[b] + outdegree[b] < k):
                                    # Uses a random uniform distribution to generate edges' weights in the [0, 1] interval
                                    matrix[node, b] = round(random.uniform(EPS, 1), 4)
                                    outdegree[node] += 1
                                    indegree[b] += 1
                                    if (indegree[b] + outdegree[b] > k) or (indegree[node] + outdegree[node] > k):
                                        print 'violation!\n'
                                    int_count += 1
                                    i += 1
                                    inserted = True
                                    break
                                    # sys.stdout.write('.')
                        # display percentage of work completed
                        threshold = int(math.floor(edge_count / 10.0))
                        percentage = int(100 * (float(i) / edge_count))
                        if i % threshold < EPS and percentage != previous_total:
                            print str(percentage) + " % ",
                            previous_total = percentage
                if not inserted:
                    print 'No more internal edge found! Count = {0}'.format(str(int_count))
                    break
                print "100 % completed."
        print 'internal edge count is {0}'.format(str(int_count))

        # tries to add more internal edges by locating the vertices whose degree is smaller than k
        for node in xrange(N):
            while int_count < internal_edge_num:
                # gets a vertex with indegree + outdegree < k
                executed = False
                for v in cluster_node_list[mycluster[node]]:
                    if int_count >= internal_edge_num:
                        break
                    if matrix[node, v] == 0 and (indegree[v] + outdegree[v] < k):
                        # creates an internal edge (node,v)
                        matrix[node, v] = round(random.uniform(EPS, 1), 4)
                        outdegree[node] += 1
                        indegree[v] += 1
                        int_count += 1
                        print 'Added one more.'
                        executed = True
                if not executed:
                    print 'No more internal edge found! Count = {0}'.format(str(int_count))
                    break
        print 'internal edge count is {0}'.format(str(int_count))
        print 'Generated {0} random internal edges.'.format(internal_edge_num)

        # 2.2- Creates the external edges (between different clusters), initially negative
        print "Generating external edges..."
        ext_count = 0
        previous_total = -1
        percentage = 0
        for i in xrange(external_edge_num):
            # Randomly picks a pair of nodes that do not belong to the same cluster
            e = pick_random_external_edge(N, indegree, outdegree, k, c)
            if e is None:
                print "Warn: no more external edge found."
                break
            # edges are directed
            # TODO: generate edges' weight randomly in the [0, 1] interval
            matrix[e.x, e.y] = round(random.uniform(-1, -EPS), 4)
            outdegree[e.x] += 1
            indegree[e.y] += 1
            ext_count += 1

            threshold = int(math.floor(external_edge_num / 10.0))
            percentage = int(100 * (float(ext_count) / external_edge_num))
            if ext_count % threshold < EPS and percentage != previous_total:
                print str(percentage) + " % ",
                previous_total = percentage
            # sys.stdout.write('.')

        print "100 % completed."
        print '\nGenerated {0} random external edges.'.format(ext_count)
        external_edge_num = ext_count

        print 'Now there are {0} internal edges and {1} external edges'.format(str(int_count), str(ext_count))
        # internal_edge_num = count
        # external_edge_num = total_edges - internal_edge_num

        # DISABLED - restarts graph generation until the correct number of edges is obtained
        success = True
        # for v in range(N):
        #     # assert degree[v] == k, "The degree of each vertex must be equal to k"
        #     if indegree[v] + outdegree[v] != k:
        #         print 'Degrees of {0}: {1} = {2} = {3}\n'.format(str(v), str(indegree[v]), str(outdegree[v]),
        #                                                          str(indegree[v] + outdegree[v]))
        #         success = False
        #         print "\nRetrying...\n"
        #         break
        #     if not success:
        #         # sleep a random number of seconds (between 0s and 3s)
        #         wait_time = randint(0, 2)
        #         print "Waiting {0} second(s)...".format(str(wait_time))
        #         sleep(wait_time)
    # end while not success

    print 'Number of internal edges of the graph is {0}'.format(str(internal_edge_num))
    print 'Number of external edges of the graph is {0}'.format(str(external_edge_num))

    # this assertion was disabled
    # total_internal_edges = 0
    # for c1 in xrange(c):
    #     for v1 in cluster_node_list[c1]:
    #         for v2 in cluster_node_list[c1]:
    #             if matrix[v1, v2] != 0:
    #                 total_internal_edges += 1
    # assert total_internal_edges == internal_edge_num, "There must be (c x n x k x p_in) edges connecting pairs of nodes of the same cluster"

    # 3- uses probability p- to generate the negative sign of the internal clusters' edges (connecting vertices inside a given cluster)
    print "Step 3: uses probability p- to generate the negative sign of the internal clusters' edges..."
    neg_links_num = int(math.ceil(internal_edge_num * p_minus))
    print "Converting %s random internal cluster edges to negative sign." % str(neg_links_num)
    # Traverses all internal clusters' edges
    count = neg_links_num
    previous_total = -1
    percentage = 0
    # equally divides the number of internal negative edges between all clusters / communities
    c1 = 0
    negative_edges_per_cluster = int(math.ceil(neg_links_num / c))
    print "There will be {0} random negative internal cluster edges per cluster.".format(negative_edges_per_cluster)
    while count > 0:
    # for c1 in xrange(c):
        negative_edges_in_cluster = 0
        random.shuffle(cluster_node_list[c1])
        # print "cluster {0} => ".format(c1),
        for v1 in cluster_node_list[c1]:
            for v2 in cluster_node_list[c1]:
                if count == 0:
                    break
                total_done = neg_links_num - count
                threshold = int(math.floor(neg_links_num / 10.0))
                percentage = int(100 * (float(total_done) / neg_links_num))
                #print str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
                if total_done % threshold < EPS and percentage != previous_total:
                    print str(percentage) + " % ",
                    previous_total = percentage
                if matrix[v1, v2] > 0:
                    matrix[v1, v2] *= -1
                    count -= 1
                    negative_edges_in_cluster += 1
                if negative_edges_in_cluster >= negative_edges_per_cluster:
                    break
            if negative_edges_in_cluster >= negative_edges_per_cluster or count == 0:
                break
        if count == 0:
            break
        # circular loop increment
        c1 = (c1 + 1) % c
    # end while count > 0
    print "100 % completed."
    assert count == 0, "3. There must be (c x n x k x pin x p-) negative links within communities"

    # 4- uses probability p+ to generate the positive sign of the external clusters' edges (connecting vertices between different clusters)
    print "Step 4: uses probability p+ to generate the positive sign of the external clusters' edges..."
    pos_links_num = int(math.ceil(external_edge_num * p_plus))
    print "Converting %s random external cluster edges to positive sign." % str(pos_links_num)
    # Traverses all external clusters' edges
    previous_total = -1
    percentage = 0
    count = pos_links_num
    conj1 = range(c)
    conj2 = range(c)
    random.shuffle(conj1)
    random.shuffle(conj2)
    # shuffles the vertices in each cluster
    for c1 in conj1:
        random.shuffle(cluster_node_list[c1])

    for c1 in conj1:
        for c2 in conj2:
            if c1 != c2 and count > 0:
                for v1 in cluster_node_list[c1]:
                    for v2 in cluster_node_list[c2]:
                        if v1 != v2 and matrix[v1, v2] < 0:
                            total_done = pos_links_num - count
                            threshold = int(math.floor(pos_links_num / 10.0))
                            percentage = int(100 * (float(total_done) / pos_links_num))
                            if total_done % threshold < EPS and percentage != previous_total:
                                print str(percentage) + " % ",
                                previous_total = percentage
                            if count > 0:
                                matrix[v1, v2] *= -1
                                count -= 1
                        if count == 0:
                            break
                    if count == 0:
                        break
                if count == 0:
                    break
        if count == 0:
            break
    print "100 % completed."
    assert count == 0, "4. There must be (c x n x k x (1 - p_in) x p+) positive links outside communities"

    # Writes output file in XPRESS Mosel format (.mos)
    # --------------------------------------------------
    # file_content[div_a] += str(vertex_a-size*div_a)+"\t"+str(vertex_b-size*div_a)+"\t"+str(value)+"\n"
    # Stores graph file in output folder
    # Vertex numbers start at 1
    filename_prefix = "c" + str(c) + "n" + str(n) + "k" + str(k) + "pin" + str(p_in) + "p-" + str(p_minus) + "p+" + str(
        p_plus)
    filename = filename_prefix + ".g"
    directory = "output"
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(directory + "/" + filename, "w") as g_file:
        # file header
        g_file.write('people: {0}\r\n\r\n'.format(str(N)))
        g_file.write('VarErr: 0.5\r\n\r\n')
        vertex_names = ""
        for x in xrange(1, N + 1):
            vertex_names += str(x) + " "
        g_file.write('Names: [{0}]\r\n\r\n'.format(vertex_names))
        # graph contents
        edge_list = ""
        for i in xrange(N):
            for j in xrange(N):
                if matrix[i, j] != 0:
                    edge_list += '({0},{1}){2} \r\n'.format(str(i + 1), str(j + 1), str(matrix[i,j]))
        g_file.write('Mrel: [ {0}]'.format(edge_list))

    # Writes output file with additional information about graph generation
    imbalance = CCObjectiveFunction(N)
    with open(directory + "/" + filename_prefix + "-info.txt", "w") as t_file:
        t_file.write('n: {0} vertices per cluster\r\n'.format(str(n)))
        t_file.write('c (clusters): {0}\r\n'.format(str(c)))
        t_file.write('N: {0} total vertices\r\n'.format(str(N)))
        t_file.write('I(P) = {0}\r\n'.format(str(imbalance)))
        for cl in xrange(c):
            t_file.write('Cluster {0}: '.format(str(cl)))
            count = 0
            for node in cluster_node_list[cl]:
                t_file.write('{0} '.format(str(node)))
                count += 1
            assert count == n
            t_file.write('\r\n')

    print "\nOutput files successfully generated: {0}".format(filename)


def main(argv):
    c = n = k = 0
    p_in = p_minus = p_plus = -1.0
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=int, required=True, help="<number of clusters>")
    parser.add_argument('-n', type=int, required=True, help="<number of vertices in each cluster>")
    parser.add_argument('-k', type=int, required=True, help="<degree of each vertex>")
    parser.add_argument('-pin', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="<probability of each node connecting other nodes in the same community>")
    parser.add_argument('-p_minus', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="<probability of negative links appearing within communities>")
    parser.add_argument('-p_plus', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="<probability of of positive links appearing between communities>")
    args = parser.parse_args()

    # measure elapsed time during graph generation
    start = time.time()

    # Call random graph generation procedure
    SG(args.c, args.n, args.k, args.pin, args.p_minus, args.p_plus)

    end = time.time()
    elapsed = end - start
    print "Graph generation took {0:.2f} seconds.".format(elapsed)

if __name__ == "__main__":
    main(sys.argv[1:])
