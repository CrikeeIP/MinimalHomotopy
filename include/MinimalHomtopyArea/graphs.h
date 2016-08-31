// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/HomotopyArea
// Distributed under the MIT Software License (X11 license)
// (See accompanying file LICENSE)



#pragma once


#include "..\Geometry\include\geometry\geometry2d.h"

#pragma warning( push )
#pragma warning( disable: 4503 )
//Boost functionality
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
//BoostGraph properties
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
//BoostGraph iterators
#include <boost/graph/iteration_macros.hpp>
//BoostGraph implementation
#include <boost/graph/adjacency_list.hpp>
//BoostAlgorithm: ConnectedComponents
#include <boost/graph/connected_components.hpp>
//BoostAlgorithm: PlanarityTest
#include <boost/graph/boyer_myrvold_planar_test.hpp>
//BoostAlgorithm: PlanarFaceTraversal
#include <boost/graph/planar_face_traversal.hpp>
//BoostDarstellung
#include <boost/graph/graphviz.hpp>
#pragma warning( pop )


/*Define the graph type
listS: selects the STL list container to store
the OutEdge list
vecS: selects the STL vector container to store
the vertices
undirectedS: selects directed edges
*/
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> BoostGraph;


/*Define an own embedded graph type, based on an actual embedding in the plane represented in coordinates of the vertices*/
typedef std::vector<std::size_t> adjacency_list_t;
typedef std::vector<adjacency_list_t> adjacency_lists_t;
typedef std::vector<geom2d::Vec2D<double>> vertex_coordinates_t;

typedef std::vector<std::size_t> vertex_indices;
typedef std::vector<vertex_indices> face_enumeration;


struct embedded_graph {
	adjacency_lists_t adjacency_lists;
	vertex_coordinates_t vertex_coordinates;
};


//Prints the given graph to std out
void output_graph( const BoostGraph& G );


//Enumerates the connected components of a the given graph
void connected_components( const BoostGraph& G );


//Tests the given Graph G for planarity. Optionally offers an embedding if G is planar
void planarity_tests( const BoostGraph& G );


//Enumerates the faces of an embedded (and thus planar) graph
face_enumeration enumerate_embedding_faces( const embedded_graph& G );
