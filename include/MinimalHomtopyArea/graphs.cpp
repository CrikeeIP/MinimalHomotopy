// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/Spielwiese
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <iostream>
#include <vector>

//#include <FunctionalPlus/fplus.h>
#include <fplus/fplus.h>

#include "graphs.h"



void output_graph( const BoostGraph& G ) {
	// represent graph in DOT format and send to cout
	std::cout << "Erzeugter Graph:" << std::endl;
	write_graphviz( std::cout, G );
}


void connected_components( const BoostGraph& G ) {
	//calculate connected components. The component of vertex i is stored in components[i]
	std::vector<int> components( boost::num_vertices( G ) );
	int num = boost::connected_components( G, &components[0] );
	//Ooutput connected components
	std::cout << "\nTotal number of components: " << num << std::endl;
	for ( std::size_t i = 0; i != components.size(); ++i )
		std::cout << "Vertex " << i << " is in component " << components[i] << std::endl;
}


void planarity_tests( const BoostGraph& G ) {
	assert( num_vertices( G ) != 0 );
	//Test for planarity
	bool is_planar = boyer_myrvold_planarity_test( G );
	std::cout << "\nGraph is " << (is_planar ? "" : "not") << " planar!\n";

	//Test for planarity - obtain an embedding on the fly:
	typedef std::vector< boost::graph_traits<BoostGraph>::edge_descriptor > vec_t;
	std::vector<vec_t> embedding( num_vertices( G ) );
	if ( boost::boyer_myrvold_planarity_test( G, &embedding[0] ) )
		std::cout << "Input graph is planar" << std::endl;
	else
		std::cout << "Input graph is not planar" << std::endl;
	//TODO: print embedding
}



struct output_visitor : public boost::planar_face_traversal_visitor
{
	std::vector<std::vector<std::size_t>> faces;
	void begin_face() { /*std::cout << "New face: ";*/ faces.push_back( {} ); }
	template <typename Vertex> void next_vertex( Vertex v ) { /*std::cout << v << " "*/; faces.back().push_back( v ); }
	void finish_face() { /*std::cout << std::endl;*/ }
};

typedef std::vector<std::vector<std::size_t>> custom_embedding;

// A typedef for the storage - a vector of vectors of edge descriptors
typedef std::vector< std::vector< boost::graph_traits<BoostGraph>::edge_descriptor > > embedding_storage_t;
// A typedef for the iterator property map, assuming the graph has an interior vertex index map
typedef boost::iterator_property_map< embedding_storage_t::iterator, boost::property_map<BoostGraph, boost::vertex_index_t>::type> embedding_t;



//TODO: lambda draus machen
template<class It>
boost::iterator_range<It> pair_range( const std::pair<It, It>& p ) {
	return boost::make_iterator_range( p.first, p.second );
}


//enumerates faces of given embedding of a graph
face_enumeration let_boost_enumerate_faces( const BoostGraph& G, const custom_embedding& ext_embedding ) {

	assert( ext_embedding.size() == boost::num_vertices( G ) );
	// Create an instance of the storage and the property map
	embedding_storage_t embedding_storage( num_vertices( G ) );

	//Fill the embedding with edges as given in the adjacency lists
	BGL_FORALL_VERTICES( v, G, BoostGraph ) {
		//std::cout << "Out edges of vertex " << v << ":" << std::endl;
		for ( auto edge_target : ext_embedding[v] ) {
			for ( auto e : pair_range( out_edges( v, G ) ) ) {
				if ( e.m_target == edge_target ) {
					//std::cout << e << std::endl;
					embedding_storage[v].push_back( e );
				}
			}
		}
	}

	//Create an instance of the property map which can be passed to any function expecting a model of PlanarEmbedding.
	embedding_t embedding( embedding_storage.begin(), boost::get( boost::vertex_index, G ) );


	//Initialize the interior edge index map
	typedef boost::graph_traits<BoostGraph>::edge_descriptor EdgeDescriptor;
	typedef std::map<EdgeDescriptor, size_t> EdgeIndexMap;
	EdgeIndexMap mapIndex;
	boost::associative_property_map<EdgeIndexMap> my_edge_index_map( mapIndex );
	boost::graph_traits<BoostGraph>::edges_size_type current_index = 0;
	BGL_FORALL_EDGES( e, G, BoostGraph )
	{
		put( my_edge_index_map, e, current_index++ );
		//std::cout << e << std::endl;
	}


	//Einen "planar_face_traversal_visitor" anlegen, dessen Funktionen aufgerufen werden wenn ein neuer Knoten / eine neue Masche (= "face") / das ende einer Mascche erreicht wird etc.
	output_visitor visitor;

	//Mit der Einbettung, dem visitor und der edge_index_map ab in den planar_face_traversal Algorithmus
	planar_face_traversal( G, embedding, visitor, my_edge_index_map );

	return visitor.faces;
}

/*TODO: warum kann man das nicht nach hier auslagern aus enumerate_embedding_faces?
embedding_storage_t create_boost_embedding(const Graph& G, const custom_embedding& ext_embedding){
assert(ext_embedding.size() == boost::num_vertices(G));
// Create an instance of the storage and the property map
embedding_storage_t embedding_storage( num_vertices( G ) );
//Fill the embedding with edges as given in the adjacency lists
BGL_FORALL_VERTICES( v, G, Graph ) {
std::cout << "Out edges of vertex " << v << ":" << std::endl;
for(auto edge_dst : ext_embedding[v]){
for ( auto e : pair_range( out_edges( v, G ) ) ) {
if(e.m_target == edge_dst){
std::cout << e << std::endl;
embedding_storage[v].push_back( e );
}
}
}
}
return embedding_storage;
}*/

/*
* Sortiert rechts rum von vertex nach (vertex + (0,1) ) aus gesehen
*/
std::vector<std::size_t> sort_by_angle( const embedded_graph& G, std::size_t vertex ) {
	auto adjacent_vertices = G.adjacency_lists[vertex];
	auto coordinates = G.vertex_coordinates;
	auto vertex_coordinates = coordinates[vertex];

	adjacent_vertices = fplus::sort_on( [&coordinates, & vertex_coordinates](std::size_t v) -> double {
		double angle = geom2d::positive_angle( geom2d::Vec2D<double>( 0.0, 1.0 ), coordinates[v] - vertex_coordinates );
		return angle;
	}, adjacent_vertices );

	return adjacent_vertices;
}

custom_embedding compute_order_based_embedding( const embedded_graph& G ) {
	const std::size_t n_vertices = G.adjacency_lists.size();
	custom_embedding embedding( n_vertices );
	for ( std::size_t vertex = 0; vertex < n_vertices; vertex++ ) {
		embedding[vertex] = sort_by_angle( G, vertex );
	}

	return embedding;
}

face_enumeration enumerate_embedding_faces( const embedded_graph& G ) {
	//instantiate the boost graph
	BoostGraph G_boost( G.adjacency_lists.size() );
	// add edges to the boost equivalent of our graph
	for ( std::size_t vertex = 0; vertex < G.adjacency_lists.size(); vertex++ ) {
		for ( auto target : G.adjacency_lists[vertex] ) {

			if ( target < vertex ) {
				continue;
				//break; //Wenn TODO: sicherstellen dass die adjazenzlisten sortiert sind
			}
			assert( target != vertex );
			boost::add_edge( vertex, target, G_boost );
		}
	}

	//compute representation of the embedding that can be translated to a boost embedding
	auto the_custom_embedding = compute_order_based_embedding( G );

	return let_boost_enumerate_faces( G_boost, the_custom_embedding );
}


