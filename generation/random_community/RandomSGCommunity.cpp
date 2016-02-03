/*  Random graph generation script, according to the article
 *  'Community Mining from Signed Social Networks'
 *  (URL: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4302742 )
 *   Author: Mario Costa Levorato Junior
 */
#include <fstream>      // fstream
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include <iterator>     // ostream_operator
#include <random>
#include <iterator>

#include <boost/tokenizer.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "RandomSGCommunity.h"

// TODO CHANGE TO APPLICATION CLI PARAMETER
// Global variables
#define EPS 0.0001

const int generation::InputMessage::TERMINATE_MSG_TAG = 0;
const int generation::InputMessage::DONE_MSG_TAG = 1;

using namespace boost::numeric::ublas;
namespace fs = boost::filesystem;
namespace mpi = boost::mpi;

// http://stackoverflow.com/questions/6942273/get-random-element-from-container-c-stl
template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

namespace generation {

using namespace boost::algorithm;
using namespace boost;
using namespace std;

RandomSGCommunity::RandomSGCommunity() {
	
}

RandomSGCommunity::~RandomSGCommunity() {
	// TODO Auto-generated destructor stub
}

// Returns a random node of in-degree + out-degree < k
long RandomSGCommunity::pick_random_vertex(vector<long>& indegree, vector<long>& outdegree, const long& k)
{
    if(vertex_list.size() == 0) {
        // BOOST_LOG_TRIVIAL(info) << "Empty vertex list!!!\n";
        return -1;
	}
	long a = *select_randomly(vertex_list.begin(), vertex_list.end());
    vertex_list.remove(a);
    while(indegree[a] + outdegree[a] >= k) {
        if(vertex_list.size() == 0) {
            // BOOST_LOG_TRIVIAL(info) << "Empty vertex list!!!\n";
            return -1;
		}
		a = *select_randomly(vertex_list.begin(), vertex_list.end());
        vertex_list.remove(a);
	}
    return a;
}

// Returns a new random edge that connects 2 nodes of the same cluster
// Each node must have in-degree + out-degree < k ???
Edge RandomSGCommunity::pick_random_internal_edge(const long& N, const long& indegree, const long& outdegree, const long& k)
{
    std::random_shuffle(vertex_list.begin(), vertex_list.end());
	
    bool found = false;
    while(not found) {
        // Selects a random node of in-degree + out-degree < k
        long a = pick_random_vertex(vertex_list, indegree, outdegree, k);
        if(a < 0)  return None;
        // BOOST_LOG_TRIVIAL(info) << "1- Selected vertex " << a;
        // Selects a neighbour of a (= same cluster) with degree < k
        // BOOST_LOG_TRIVIAL(info) << "mycluster of " << a << " is " << mycluster[a];
        vector<long> neighbour_list = cluster_node_list[mycluster[a]];
		std::random_shuffle(neighbour_list.begin(), neighbour_list.end());
        for(vector<long>::iterator v_it = neighbour_list.begin(); v_it != neighbour_list.end(); v_it++)
		{
			long b = *v_it;
            if(b != a and matrix(a, b) == 0 and indegree[b] + outdegree[b] < k) {
                found = true;
                break;
			}
		}
	}
    return Edge(a, b);
}

// Returns a new random edge that connects 2 nodes of different clusters
// Each node must have degree < k
Edge RandomSGCommunity::pick_random_external_edge(const long& N, const long& indegree, const long& outdegree, 
		const long& k, const long& c)
{
    std::vector<long> cluster_list;
	boost::push_back(cluster_list, boost::irange(0, c));
    std::random_shuffle(cluster_list.begin(), cluster_list.end());

    // removes the vertices that have exceeded their maximum degrees
	for(vector<long>::iterator v_it = vertex_list.begin(); v_it != vertex_list.end(); v_it++)
	{
		long a = *v_it;
		if(indegree[a] + outdegree[a] >= k) {
            vertex_list.remove(a);
            // BOOST_LOG_TRIVIAL(info) << "removed node " << a;
		}
	}
    std::random_shuffle(vertex_list.begin(), vertex_list.end());

    bool found = false;
    long i = 0;
    while(i < vertex_list.size()) {
        // Selects a random node of degree < k which belongs to vertex_list
        if(vertex_list.size() == 0)  return None;
        long a = vertex_list[i];
        i++;
        // BOOST_LOG_TRIVIAL(info) << 'selected node {0}\n'.format(str(a))
        // Selects a vertex b, non-neighbour of a (different cluster) with degree < k
        // vertex b belongs to cluster c1 in cluster_list
		for(long c1 = 0; c1 < cluster_node_list.size(); c1++)
		{
			if(c1 != mycluster[a]) {
				vector<long> nodeList = cluster_node_list[c1];
				std::random_shuffle(nodeList.begin(), nodeList.end());
				for(vector<long>::iterator v_it = nodeList.begin(); v_it != nodeList.end(); v_it++)
				{
					long b = *v_it;
					if(matrix(a, b) == 0 and (indegree[b] + outdegree[b] < k) {
                        found = true;
                        break;
					}
				}
			}
            if(found)  break;
		}
        if(found)  break;
    } // end while
    if(not found) {
        BOOST_LOG_TRIVIAL(info) << "b not found!";
        return None;
	} else {
        return Edge(a, b);
	}
}

double RandomSGCommunity::CCObjectiveFunction(const long& N)
{
    double positiveSum = double(0.0);
    double negativeSum = double(0.0);

    // For each vertex i
    for(long i = 0; i < N; i++) {
        // Find out to which cluster vertex i belongs
        long k = mycluster[i];
        // For each out edge of i
        for(long j = 0; j < N; j++) {
            if(matrix(i, j) != 0) {
                double weight = matrix(i, j);
                bool sameCluster = (k == mycluster[j]);
                if(weight < 0 and sameCluster) {  // negative edge
                    // i and j are in the same cluster
                    negativeSum += weight * (-1);
				} else if (weight > 0 and not sameCluster) {  // positive edge
                    // i and j are NOT in the same cluster
                    positiveSum += weight;
				}
			}
		}
	}

    BOOST_LOG_TRIVIAL(info) << 'CC calculated value. Obj = {0}\n'.format(str(positiveSum + negativeSum))
    return positiveSum + negativeSum;
}

void RandomSGCommunity::SG(c, n, k, p_in, p_minus, p_plus) 
{
    // Variables used in graph generation
    long N = n * c;
    BOOST_LOG_TRIVIAL(info) << "Generating random graph with " << N << " vertices...";
    bool success = false;

    // measure elapsed time during graph generation
    start = time.time()

    while(not success)
	{
        // Array that controls to which cluster a vertex belongs (size = N)
		mycluster = vector<long>(N, 0);
        // List of nodes inside a cluster
        //cluster_node_list = [[] for x in xrange(c)]
        // Graph adjacency matrix (with dimensions N x N)
        // Creates a list containing N lists initialized to 0 (index begins with 0)
        matrix = std::vector< std::vector< std::pair<long, char> > >(N);
        // Arrays that control the degree of each vertex
        indegree = vector<long>(N, 0);
        outdegree = vector<long>(N, 0);

        // builds a list with all the N vertices of the graph
		std::list<long> vertex_list;
		boost::push_back(vertex_list, boost::irange(0, N));
		
        // 1- randomly chooses n vertices to insert in each of the c clusters
        BOOST_LOG_TRIVIAL(info) << "Step 1: randomly chooses n vertices to insert in each of the c clusters...";
        for(long cluster_number = 0; cluster_number < c; cluster_number++) {
			for(long i = 0; i < n; i++) {
                // randomly removes a vertex from the vertex_list
                long v = *select_randomly(vertex_list.begin(), vertex_list.end());
                vertex_list.remove(v);
                mycluster[v] = cluster_number;
                cluster_node_list[cluster_number].push_back(v);
			}
		}
        assert(vertex_list.size() == 0 && "All vertices must belong to a cluster");

        // assertion tests / sanity check
        cluster_count = vector<long>(c, 0);
        for(long v = 0; v < N; v++) {
            cluster_count[mycluster[v]] += 1;
		}
        for(long x = 0; x < c; x++) {
            assert(cluster_count[x] == n && "There must be n vertices in each cluster");
            assert(len(cluster_node_list[x]) == n && "There must be n vertices in each cluster");
		}
        assert(vertex_list.size() == 0 && "Vertex list must be empty after cluster-vertex distribution");

        // 2- uses probability p_in to generate k edges (initially positive) connecting each of the vertices
        // p_in % of the edges connect nodes of the same cluster
        BOOST_LOG_TRIVIAL(info) << "Step 2: uses probability p_in to generate k edges (initially positive) connecting each of the vertices...";
        long total_edges = long(std::ceil((c * n * k) / double(2.0)));
        long internal_edge_num = long(std::ceil(total_edges * p_in));
        long external_edge_num = total_edges - internal_edge_num;
        BOOST_LOG_TRIVIAL(info) << "Number of total edges of the graph is " << total_edges;
        BOOST_LOG_TRIVIAL(info) << "Number of expected internal edges of the graph is " << internal_edge_num;
        BOOST_LOG_TRIVIAL(info) << "Number of expected external edges of the graph is " << external_edge_num;
        assert((internal_edge_num + external_edge_num == total_edges) && "Sum of internal and external edges do not match");
        boost::push_back(vertex_list, boost::irange(0, N));
        std::random_shuffle(vertex_list.begin(), vertex_list.end());

        // Guarantees that each node has exactly k edges (TODO: in-degree + out-degree = k ???)
        // 2.1- Creates the internal edges (within same cluster), initially positive
        BOOST_LOG_TRIVIAL(info) << "Generating internal edges...";
        // count = 0
        // for i in xrange(internal_edge_num):
        // Randomly picks a pair of nodes that belong to the same cluster
        // e = pick_random_internal_edge(N, indegree, outdegree, k)
        // if e == None:
        //    BOOST_LOG_TRIVIAL(info) << "Warn: no more internal edge found."
        //    break
        // edges are directed
        // matrix[e.x][e.y] = 1
        // outdegree[e.x] += 1
        // indegree[e.y] += 1
        // count += 1
        // for each node, generate an internal cluster edge
        long int_count = 0;
        // each cluster will get (internal_edge_num) / c internal edges
        long edge_count = int(math.floor((internal_edge_num) / c));
        long rest = internal_edge_num - (c * edge_count);
        while(int_count < c * edge_count) {
            inserted = False
            BOOST_LOG_TRIVIAL(info) << "there are " << int_count << " internal edges so far";
            for(long cluster = 0; cluster < c; cluster++) {
                BOOST_LOG_TRIVIAL(info) << "must generate " << edge_count << " internal edges for cluster " << cluster;
                long i = 0;
                int previous_total = -1;
                int percentage = 0;
                while(i < edge_count) {
					std::vector<long> nodeList = cluster_node_list[cluster];
                    std::shuffle(nodeList.begin(), nodeList.end());
					for(std::vector<long>::iterator v_it = nodeList.begin(); v_it != nodeList.end(); v_it++) {
						node = *v_it;
                        if(i >= edge_count)  break;
                        if(indegree[node] + outdegree[node] < k) {
                            // selects another node in the same cluster
                            // BOOST_LOG_TRIVIAL(info) << 'selected node {0} with outdegree = {1}\n'.format(str(node), str(outdegree[node]))
                            std::vector<long> neighbour_list = cluster_node_list[mycluster[node]];
							std::shuffle(neighbour_list.begin(), neighbour_list.end());
							for(std::vector<long>::iterator vb_it = neighbour_list.begin(); vb_it != neighbour_list.end(); vb_it++) {
								long b = *vb_it;
                                if(indegree[node] + outdegree[node] >= k)  break;
                                if(b != node and matrix(node, b) == 0 and (indegree[b] + outdegree[b] < k)) {
                                    // Uses a random uniform distribution to generate edges' weights in the [0, 1] interval - DISABLED
                                    matrix(node, b) = 1; // round(random.uniform(EPS, 1), 4)
                                    outdegree[node]++;
                                    indegree[b]++;
                                    if((indegree[b] + outdegree[b] > k) or (indegree[node] + outdegree[node] > k)) {
                                        BOOST_LOG_TRIVIAL(error) << "violation!";
									}
                                    int_count++;
                                    i++;
                                    inserted = true;
                                    break;
								}
							}
                        }
                        // display percentage of work completed
                        threshold = int(math.floor(edge_count / 10.0));
                        percentage = int(std::ceil(100 * (double(i) / edge_count)));
                        if(i % threshold < EPS and percentage != previous_total) {
                            BOOST_LOG_TRIVIAL(info) << percentage << " % ";
                            previous_total = percentage;
						}
					}
				}
                if(not inserted) {
                    BOOST_LOG_TRIVIAL(info) << "No more internal edge found! Count = " << int_count;
                    break;
				}
                BOOST_LOG_TRIVIAL(info) << " completed.";
			}
		}
        BOOST_LOG_TRIVIAL(info) << "internal edge count is " << int_count;

        // tries to add more internal edges by locating the vertices whose degree is smaller than k
        for(long node = 0; node < N; node++) {
            while(int_count < internal_edge_num) {
                // gets a vertex with indegree + outdegree < k
                executed = false;
				std::vector<long> nodeList = cluster_node_list[mycluster[node]];
				for(std::vector<long>::iterator v_it = nodeList.begin(); v_it != nodeList.end(); v_it++) {
					long v = *v_it;
					if(int_count >= internal_edge_num)  break;
                    if(matrix(node, v) == 0 and (indegree[v] + outdegree[v] < k)) {
                        // creates an internal edge (node,v)
                        matrix(node, v) = 1; // round(random.uniform(EPS, 1), 4)
                        outdegree[node] += 1;
                        indegree[v] += 1;
                        int_count += 1;
                        BOOST_LOG_TRIVIAL(info) << "Added one more.";
                        executed = true;
					}
				}
                if(not executed) {
                    BOOST_LOG_TRIVIAL(info) << "No more internal edge found! Count = " << int_count;
                    break;
				}
			}
		}
        BOOST_LOG_TRIVIAL(info) << "internal edge count is " << int_count;
        BOOST_LOG_TRIVIAL(info) << "Generated " << internal_edge_num << " random internal edges.";

        // 2.2- Creates the external edges (between different clusters), initially negative
        BOOST_LOG_TRIVIAL(info) << "Generating external edges...";
        long ext_count = 0;
        int previous_total = -1;
        int percentage = 0;
        for(long i = 0; i < external_edge_num; i++) {
            // Randomly picks a pair of nodes that do not belong to the same cluster
            Edge e = pick_random_external_edge(N, indegree, outdegree, k, c);
            if e is None: {
                BOOST_LOG_TRIVIAL(warn) << "Warn: no more external edge found.";
                break;
			}
            // edges are directed
            // TODO: generate edges' weight randomly in the [0, 1] interval
            matrix[e.x, e.y] = -1; // round(random.uniform(-1, -EPS), 4)
            outdegree[e.x] += 1;
            indegree[e.y] += 1;
            ext_count += 1;

            threshold = int(std::floor(external_edge_num / double(10.0)));
            percentage = int(std::ceil(100 * (double(ext_count) / external_edge_num)));
            if(ext_count % threshold < EPS and percentage != previous_total) {
                BOOST_LOG_TRIVIAL(info) << percentage << " % ";
                previous_total = percentage;
			}
		}

        BOOST_LOG_TRIVIAL(info) << " completed."
        BOOST_LOG_TRIVIAL(info) << "Generated " << ext_count << " random external edges.";
        external_edge_num = ext_count;

        BOOST_LOG_TRIVIAL(info) << "Now there are " << int_count << " internal edges and " << ext_count << " external edges";
        // internal_edge_num = count
        // external_edge_num = total_edges - internal_edge_num

        // DISABLED - restarts graph generation until the correct number of edges is obtained
        success = true;
	} // end while not success

    BOOST_LOG_TRIVIAL(info) << "Number of internal edges of the graph is " << internal_edge_num;
    BOOST_LOG_TRIVIAL(info) << "Number of external edges of the graph is " << external_edge_num;

    // this assertion was disabled
    // total_internal_edges = 0
    // for c1 in xrange(c):
    //     for v1 in cluster_node_list[c1]:
    //         for v2 in cluster_node_list[c1]:
    //             if matrix[v1, v2] != 0:
    //                 total_internal_edges += 1
    // assert total_internal_edges == internal_edge_num, "There must be (c x n x k x p_in) edges connecting pairs of nodes of the same cluster"

    // 3- uses probability p- to generate the negative sign of the internal clusters' edges (connecting vertices inside a given cluster)
    BOOST_LOG_TRIVIAL(info) << "Step 3: uses probability p- to generate the negative sign of the internal clusters' edges...";
    long neg_links_num = long(std::ceil(internal_edge_num * p_minus));
    BOOST_LOG_TRIVIAL(info) << "Converting " << neg_links_num << " random internal cluster edges to negative sign.";
    // Traverses all internal clusters' edges
    long count = neg_links_num;
    int previous_total = -1;
    int percentage = 0;
    // equally divides the number of internal negative edges between all clusters / communities
    long c1 = 0;
    long negative_edges_per_cluster = long(std::ceil(neg_links_num / double(c)));
    BOOST_LOG_TRIVIAL(info) << "There will be " << negative_edges_per_cluster << " random negative internal cluster edges per cluster.";
    while(count > 0) {
    // for c1 in xrange(c):
        long negative_edges_in_cluster = 0;
		// BOOST_LOG_TRIVIAL(info) << "cluster {0} => ".format(c1),
        std::vector<long> nodeList1 = cluster_node_list[c1];
		std::shuffle(nodeList1.begin(), nodeList1.end());
		std::vector<long> nodeList2 = cluster_node_list[c1];
		std::shuffle(nodeList2.begin(), nodeList2.end());
		for(std::vector<long>::iterator v1_it = nodeList1.begin(); v1_it != nodeList1.end(); v1_it++) {
			long v1 = *v1_it;
            for(std::vector<long>::iterator v2_it = nodeList2.begin(); v2_it != nodeList2.end(); v2_it++) {
				long v2 = *v2_it;
                if(count == 0)  break;
                total_done = neg_links_num - count;
                threshold = int(math.floor(neg_links_num / double(10.0)));
                percentage = int(std::ceil(100 * (double(total_done) / neg_links_num)));
                // BOOST_LOG_TRIVIAL(info) << str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
                if(total_done % threshold < EPS and percentage != previous_total) {
                    BOOST_LOG_TRIVIAL(info) << percentage << " % ";
                    previous_total = percentage;
				}
                if(matrix(v1, v2) > 0) {
                    matrix(v1, v2) *= -1;
                    count -= 1;
                    negative_edges_in_cluster += 1;
				}
                if(negative_edges_in_cluster >= negative_edges_per_cluster)  break;
			}
            if(negative_edges_in_cluster >= negative_edges_per_cluster or count == 0)  break;
		}
        if(count == 0)  break;
        // circular loop increment
        c1 = (c1 + 1) % c;
	} // end while count > 0
    BOOST_LOG_TRIVIAL(info) << " completed.";
    assert(count == 0 && "3. There must be (c x n x k x pin x p-) negative links within communities");

    // 4- uses probability p+ to generate the positive sign of the external clusters' edges (connecting vertices between different clusters)
    BOOST_LOG_TRIVIAL(info) << "Step 4: uses probability p+ to generate the positive sign of the external clusters' edges...";
    long pos_links_num = long(std::ceil(external_edge_num * p_plus));
    BOOST_LOG_TRIVIAL(info) << "Converting " << pos_links_num << " random external cluster edges to positive sign.";
    // Traverses all external clusters' edges
    previous_total = -1;
    percentage = 0;
    count = pos_links_num;
	std::vector<long> conj1, conj2;
	boost::push_back(conj1, boost::irange(0, c));
    std::random_shuffle(conj1.begin(), conj1.end());
    boost::push_back(conj2, boost::irange(0, c));
    std::random_shuffle(conj2.begin(), conj2.end());
    // shuffles the vertices in each cluster
    for(long c1 = 0; c1 < c; c1++) {
        std::shuffle(cluster_node_list[c1].begin(), cluster_node_list[c1].end());
	}

	for(std::vector<long>::iterator c1_it = conj1.begin(); c1_it != conj1.end(); c1_it++) {
		long c1 = *c1_it;
		for(std::vector<long>::iterator c2_it = conj2.begin(); c2_it != conj2.end(); c2_it++) {
			long c2 = *c2_it;
			if(c1 != c2 and count > 0) {
				std::vector<long> nodeList1 = cluster_node_list[c1];
				std::vector<long> nodeList2 = cluster_node_list[c1];
				for(std::vector<long>::iterator v1_it = nodeList1.begin(); v1_it != nodeList1.end(); v1_it++) {
					long v1 = *v1_it;
					for(std::vector<long>::iterator v2_it = nodeList2.begin(); v2_it != nodeList2.end(); v2_it++) {
						long v2 = *v2_it;
                        if(v1 != v2 and matrix(v1, v2) < 0) {
                            total_done = pos_links_num - count;
                            threshold = int(math.floor(pos_links_num / double(10.0)));
                            percentage = int(std::ceil(100 * (double(total_done) / pos_links_num)));
                            if(total_done % threshold < EPS and percentage != previous_total) {
                                BOOST_LOG_TRIVIAL(info) << percentage << " % ";
                                previous_total = percentage;
							}
                            if(count > 0) {
                                matrix(v1, v2) *= -1;
                                count -= 1;
							}
						}
                        if(count == 0)  break;
					}
                    if(count == 0)  break;
				}
                if(count == 0)  break;
			}
		}
        if(count == 0)  break;
	}
    BOOST_LOG_TRIVIAL(info) << " completed.";
    assert(count == 0 && "4. There must be (c x n x k x (1 - p_in) x p+) positive links outside communities");

    // Writes output file in XPRESS Mosel format (.mos)
    // --------------------------------------------------
    // file_content[div_a] += str(vertex_a-size*div_a)+"\t"+str(vertex_b-size*div_a)+"\t"+str(value)+"\n"
    // Stores graph file in output folder
    // Vertex numbers start at 1
    filename_prefix = "c" + str(c) + "n" + str(n) + "k" + str(k) + "pin" + str(p_in) + "p-" + str(p_minus) + "p+" + str(
        p_plus)
    filename = filename_prefix + ".g"
    directory = "output"
    if not os.path.exists(directory):
        os.makedirs(directory)
    BOOST_LOG_TRIVIAL(info) << "Saving output graph file...";

    // write the matrix to matrix market file format
    scipy.io.mmwrite(directory + "/" + filename, matrix)

    BOOST_LOG_TRIVIAL(info) << "Graph file generated: " << filename;
    end = time.time()
    elapsed = end - start
    BOOST_LOG_TRIVIAL(info) << "Graph generation took " << elapsed << " seconds.";

    // Writes output file with additional information about graph generation
    double imbalance = CCObjectiveFunction(N);
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
            assert(count == n);
            t_file.write('\r\n')

    BOOST_LOG_TRIVIAL(info) << "Info file generated: " << filename_prefix << "-info.txt";

    BOOST_LOG_TRIVIAL(info) << "\nOutput files successfully generated: " << filename;
}

bool RandomSGCommunity::generateRandomSG(const long& c, const long& n, const long& k,
		const double& p_in, const double& p_minus, const double& p_plus,
		const unsigned int &myRank, const unsigned int &numProcessors) {
	BOOST_LOG_TRIVIAL(info) << "Processing folder " << folder;
	
	// read the files in the specified folder
	fs::path inputDir(folder);
	fs::directory_iterator end_iter;
	if (!fs::exists(inputDir) || !fs::is_directory(inputDir)) {
		BOOST_LOG_TRIVIAL(fatal) << "Input file directory not found. Exiting." << endl;
		return false;
	}
	std::vector<fs::path> fileList;
	for( fs::directory_iterator dir_iter(inputDir) ; dir_iter != end_iter ; ++dir_iter) {
		string ext = dir_iter->path().extension().string();
		string filename = dir_iter->path().filename().string();
		boost::algorithm::to_lower(filename);
		if ((fs::is_regular_file(dir_iter->status())) &&
				filename == filter) {
			fs::path filePath = *dir_iter;
			fileList.push_back(filePath);
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Found " << fileList.size() << " file(s).";
	// output folder
	stringstream path_ss;
	path_ss << folder << boost::filesystem::path::preferred_separator << "unweightedSG";
	boost::filesystem::path outputPath(path_ss.str());
	boost::system::error_code returnedError;
	boost::filesystem::create_directories( outputPath, returnedError );
	
	// read the movielens dataset voting file from CSV
	for(unsigned int i = 0; i < fileList.size(); i++) {
		boost::timer::cpu_timer timer;
		boost::timer::cpu_times start_time, end_time;
		double timeSpent = 0.0;
		// measure conversion execution time
		timer.start();
		start_time = timer.elapsed();

		fs::path filePath = fileList.at(i);
		string filename = filePath.parent_path().filename().string();
		stringstream file_ss;
		file_ss << path_ss.str() << boost::filesystem::path::preferred_separator << filename << ".g";

		long max_user_id = 0, max_movie_id = 0;
		readMovieLensCSVFile(filePath.string(), max_user_id, max_movie_id);
		// process the dataset file and generate signed graph
		if(generateSGFromMovieRatings(max_user_id, max_movie_id, file_ss.str(), myRank, numProcessors)) {
			if(myRank == 0) {
				BOOST_LOG_TRIVIAL(info) << "\nCreated output signed graph file " << file_ss.str();
			}
			// Stops the timer and stores the elapsed time
			timer.stop();
			end_time = timer.elapsed();
			timeSpent = (end_time.wall - start_time.wall) / double(1000000000);
			BOOST_LOG_TRIVIAL(info) << "Graph generation took " << timeSpent << " seconds.";
		} else {
			BOOST_LOG_TRIVIAL(fatal) << "Error generating signed graph file " << file_ss.str() << endl;
		}
	}
	BOOST_LOG_TRIVIAL(info) << "Done.";
	return true;
}

bool RandomSGCommunity::generateSGFromMovieRatings(const long& max_user_id, const long& max_movie_id,
		const string& outputFileName, const unsigned int &myRank, const unsigned int &numProcessors) {
	// the signed graph representing the voting in movie lens dataset
	int previous_total = -1;
	int percentage = 0;

	// compressed_matrix<long> common_rating_count(max_user_id, max_user_id);
	// compressed_matrix<long> common_similar_rating_count(max_user_id, max_user_id);
	BOOST_LOG_TRIVIAL(info) << "Processing graph with " << max_user_id << " users and " << max_movie_id << " movies...";

	// Parallel processing with MPI *************************************
	// Splits the processing in (n / numProcessors) chunks,
	// to be consumed by numProcessors processes
	long MPIChunkSize = long(std::floor(double(max_user_id) / numProcessors));
	long MPIRemainingVertices = max_user_id % numProcessors;
	long initialProcessorUserIndex = myRank * MPIChunkSize;
	long finalProcessorUserIndex = (myRank + 1) * MPIChunkSize - 1;
	if(MPIRemainingVertices > 0 and myRank == numProcessors - 1) {
		finalProcessorUserIndex = max_user_id - 1;
		MPIChunkSize = finalProcessorUserIndex - initialProcessorUserIndex + 1;
	}
	BOOST_LOG_TRIVIAL(info) << "This is process " << myRank << ". Will process user range [" << initialProcessorUserIndex << ", " << finalProcessorUserIndex << "]";

	// Begin dividing the MPI Chunk into smaller chunks for local processing
	long chunkSize = long(std::floor(double(MPIChunkSize) / NUMBER_OF_CHUNKS));
	long remainingVertices = MPIChunkSize % NUMBER_OF_CHUNKS;
	std::ostringstream out;
	long edgeCount = 0;
	for(int i = 0; i < NUMBER_OF_CHUNKS; i++) {
		// will process only the range (initialUserIndex <= user_a <= finalUserIndex)
		long initialUserIndex = initialProcessorUserIndex + i * chunkSize;
		long finalUserIndex = initialProcessorUserIndex + (i + 1) * chunkSize - 1;
		if(remainingVertices > 0 and i == NUMBER_OF_CHUNKS - 1) {
			finalUserIndex = finalProcessorUserIndex;
			chunkSize = finalUserIndex - initialUserIndex + 1;
		}
		BOOST_LOG_TRIVIAL(info) << "\nProcessing user range [" << initialUserIndex << ", " << finalUserIndex << "]";
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > > common_rating_count(chunkSize, max_user_id);;
		generalized_vector_of_vector< long, row_major, boost::numeric::ublas::vector<mapped_vector<long> > > common_similar_rating_count(chunkSize, max_user_id);

		BOOST_LOG_TRIVIAL(info) << "Begin movie list traversal (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")...";
		cout << "\nBegin movie list traversal (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")..." << endl;
		for(long movie_id = 0; movie_id < max_movie_id; movie_id++) {
			std::vector< std::pair<long, char> > ratingList = movie_users[movie_id];
			// cout << movie_id << " with size " << ratingList.size() << endl;
			for(long x = 0; x < ratingList.size(); x++) {
				std::pair<long, char> user_rating_x = ratingList[x];
				long user_a = user_rating_x.first;
				char rating_a = user_rating_x.second;
				if(initialUserIndex <= user_a and user_a <= finalUserIndex) {
					for(long y = x + 1; y < ratingList.size(); y++) {
						std::pair<long, char> user_rating_y = ratingList[y];
						long user_b = user_rating_y.first;
						char rating_b = user_rating_y.second;
						if(user_a != user_b) {
							common_rating_count(user_a - initialUserIndex, user_b) += 1;
							if((rating_a <= BAD_MOVIE and rating_b <= BAD_MOVIE)
									or (rating_a >= GOOD_MOVIE and rating_b >= GOOD_MOVIE)
									or (rating_a == REGULAR_MOVIE and rating_b == REGULAR_MOVIE)) {
								// users A and B have the same opinion about the movie
								common_similar_rating_count(user_a - initialUserIndex, user_b) += 1;
							}
						}
					}
				}
			}
			// display status of processing done
			int threshold = int(std::floor(max_movie_id / 10.0));
			int percentage = int(std::ceil(100 * (double(movie_id) / max_movie_id)));
			if((movie_id % threshold < EPS) and (percentage != previous_total)) {
				cout << percentage << " % ";
				cout.flush();
				BOOST_LOG_TRIVIAL(info) << percentage << " % ";
				previous_total = percentage;
			}
		}

		BOOST_LOG_TRIVIAL(info) << "\nBegin edge generation (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")...";
		cout << "\nBegin edge generation (step " << (i+1) << " of " << NUMBER_OF_CHUNKS << ")..." << endl;
		long count = (finalUserIndex - initialUserIndex + 1) * max_user_id;
		previous_total = -1;
		percentage = 0;
		for(long user_a = initialUserIndex; user_a <= finalUserIndex; user_a++) {
			for(long user_b = 0; user_b < max_user_id; user_b++) {
				if(user_a != user_b) {
					if(common_rating_count(user_a - initialUserIndex, user_b) > 0) {
						double common_similar_rating_ratio = double(common_similar_rating_count(user_a - initialUserIndex, user_b)) /
								common_rating_count(user_a - initialUserIndex, user_b);
						if(common_similar_rating_ratio >= POS_EDGE_PERC) {
							// SG[user_a, user_b] = 1;
							out << user_a << " " << user_b << " 1\n";
							edgeCount++;
						}
						else if(common_similar_rating_ratio <= NEG_EDGE_PERC) {
							// SG[user_a, user_b] = -1;
							out << user_a << " " << user_b << " -1\n";
							edgeCount++;
						}
					}
				}
				// display status of processing done
				int total_done = (user_a - initialUserIndex) * max_user_id + user_b;
				int threshold = int(std::floor(count / 10.0));
				int percentage = int(std::ceil(100 * (double(total_done) / count)));
				// BOOST_LOG_TRIVIAL(info) << str(total_done) + " % " + str(threshold) + " = " + str(total_done % threshold)
				if((total_done % threshold < EPS) and (percentage != previous_total)) {
					cout << percentage << " % ";
					cout.flush();
					BOOST_LOG_TRIVIAL(info) << percentage << " % ";
					previous_total = percentage;
				}
			}
		}
	}

	// write the signed graph to output file
	stringstream partialFilename_ss;
	partialFilename_ss << outputFileName << ".part" << myRank;
	ofstream output_file(partialFilename_ss.str().c_str(), ios::out | ios::trunc);
	if(!output_file) {
		BOOST_LOG_TRIVIAL(fatal) << "Cannot open output result file to: " << outputFileName;
		return false;
	}
	if(myRank == 0) {
		output_file << max_user_id << "\t" << edgeCount << "\r\n";
	}
	output_file << out.str();
	// Close the file
	output_file.close();

	boost::mpi::communicator world;
	if(myRank == 0) {
		// wait for the message 'done' from all workers
		for(int i = 1; i < numProcessors; i++) {
			InputMessage imsg;
			boost::mpi::status msg = world.recv(boost::mpi::any_source, InputMessage::DONE_MSG_TAG, imsg);
		}
		// merge all the output files into a full graph output file
		ofstream output_full_file(outputFileName.c_str(), ios::out | ios::trunc);
		if(!output_full_file) {
			BOOST_LOG_TRIVIAL(fatal) << "Cannot open full output result file to: " << outputFileName;
			return false;
		}
		for(int i = 0; i < numProcessors; i++) {
			stringstream partialFilenameN_ss;
			partialFilenameN_ss << outputFileName << ".part" << i;
			ifstream in(partialFilenameN_ss.str().c_str());
			if (!in.is_open()){
				BOOST_LOG_TRIVIAL(fatal) << "Error opening output file number " << i;
				return false;
			}
			string fileContents = get_file_contents(partialFilenameN_ss.str().c_str());
			output_full_file << fileContents;
		}
		output_full_file.close();

	} else {
		// send a message to the leader process to inform that the task is done
		BOOST_LOG_TRIVIAL(info) << "Sending done message to leader...";
		InputMessage imsg;
		world.send(0, InputMessage::DONE_MSG_TAG, imsg);
	}

	return true;
}

std::string RandomSGCommunity::get_file_contents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::ostringstream contents;
    contents << in.rdbuf();
    in.close();
    return(contents.str());
  }
  throw(errno);
}

void RandomSGCommunity::find_and_replace(string& source, string const& find, string const& replace)
{
    for(string::size_type i = 0; (i = source.find(find, i)) != string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length();
    }
}

} /* namespace generation */
