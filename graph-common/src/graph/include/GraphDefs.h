#ifndef GRAPH_DEFS_H_
#define GRAPH_DEFS_H_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include "../../util/serialization/dynamic_bitset.hpp"

struct Edge {
    double weight;
    std::size_t vertex_index_t;
    Edge() : weight(0), vertex_index_t(0) { }
    Edge(double w) : weight(w), vertex_index_t(0) { }

    friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & weight;
		ar & vertex_index_t;
	}
};
struct Vertex {
    long id;
    std::size_t edge_index_t;
    std::size_t vertex_rank_t;

    Vertex() : id(0), edge_index_t(0) { }
    Vertex(long w) : id(w), edge_index_t(0), vertex_rank_t(0)  { }

    friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, unsigned int file_version)
	{
		ar & id;
		ar & edge_index_t;
		ar & vertex_rank_t;
	}
};


enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, properties);
	BOOST_INSTALL_PROPERTY(edge, properties);
}

#endif /* GRAPH_H_ */

