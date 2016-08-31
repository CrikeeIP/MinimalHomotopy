// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 licnse).
// (See accompanying file LICENSE)

#include <iostream>
#include <assert.h>

#include "..\Geometry\include\geometry\geometry2d.h"
#include "graphs.h"
#include "homotopy_area_calculator.h"


using namespace geom2d;
using namespace min_ht_area;



struct intersection {
	std::size_t p1_index;
	std::size_t p2_index;
	double t1;
	double t2;
};


point point_on_segment(const point& x1, const point& x2, double t ) {
	point vec = (x2-x1) * t;
	return x1 + vec;
}



/*Stores in field lut[i][j] (i<p1.size()-1, j<p2.size()-1) the parameter t (in [0,1]) of the intersection of segment i = (p1[i],p1[i+1]) =: (a,b) on p1 and j on p2.
* (Here, p1 and p2 refer to the unrefined polygons)
* t == -1.0 if they do not intersect
* otherwise, x = a + t*(b-a) is the intersection of p1 and p2
*/
std::vector<std::pair<intersection, std::size_t>> segment_intersection_lut;
/**Enumerates the vertices which are intersections of p1 and p2, i.e., which are vertices of both (refined) polygons
*  intersection_indices[i] stores the indices of this vertex on p1 (in .first) and on p2 (in .second)
*/
std::vector<std::pair<std::size_t, std::size_t>> intersection_indices;

//Computes the refinement of p2
polygon build_p2(const polygon& p1, const polygon& p2_orig) {
	polygon p2;
	
	//sort by p2-index und time in [0,1]
	std::vector<std::pair<intersection, std::size_t>> intersections = fplus::sort_by( []( const std::pair<intersection, std::size_t>& i1, const std::pair<intersection, std::size_t>& i2 ) -> bool {
		return i1.first.p2_index == i2.first.p2_index ? i1.first.t2 < i2.first.t2 : i1.first.p2_index < i2.first.p2_index;
	}, segment_intersection_lut );
	
	std::size_t intersect_idx = 0;
	//insert all original vertices and the additional intersection vertices into the refinement p2 of p2_orig
	for ( std::size_t p2_idx = 0; p2_idx <= p2_orig.size()-2; p2_idx++ ) {
		p2.push_back( p2_orig[p2_idx] );
		double last_time = -1.0;

		while ( intersect_idx < intersections.size() &&  intersections[intersect_idx].first.p2_index == p2_idx ) {
			const auto& i = intersections[intersect_idx++];
			if ( i.first.t2 == 0.0 ) {
				intersection_indices[i.second].second = p2.size() - 1;
				last_time = i.first.t2;
				continue;  ///do not insert intersections which are already vertices of p1/p2
			}
			else if( i.first.t2 == 1.0 || i.first.t2 == last_time ) {
				last_time = i.first.t2;
				continue;  ///do not insert intersections which are already vertices of p2
			}
			last_time = i.first.t2;
			p2.push_back( p1[intersection_indices[i.second].first] );
			intersection_indices[i.second].second = p2.size() - 1;
		}
	}
	p2.push_back( p2_orig.back() );
	intersection_indices.back().second = p2.size() - 1;

	return p2;
}


std::pair<polygon, polygon> initialize_intersection_lut( const polygon& p1_orig, const polygon& p2_orig ) {
    ///Optional: TODO: Instead of O(n*m) brute force probably some red-blue line/trapezoid sweep? (e.g. Chan http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.44.4227&rep=rep1&type=pdf,  but problem with degeneracies..)
	segment_intersection_lut = std::vector<std::pair<intersection,std::size_t>>();
	intersection_indices = std::vector<std::pair<std::size_t, std::size_t>>();
	polygon p1, p2;

	intersection_indices.push_back( { 0, 0 } );

	//for all segments of p1
	for ( std::size_t p1_i = 0; p1_i <= p1_orig.size() - 2; p1_i++ ) {
		p1.push_back( p1_orig[p1_i] );
		//compute intersections with segments of p2
		std::vector<intersection> intersection_times_p1i;
		///intersection_times.reserve( p2.size() ); ///Optional: TODO: Useful in case of performance issues? Else delete
		for ( std::size_t p2_i = 0; p2_i <= p2_orig.size() - 2; p2_i++ ) {
			auto intersection_time = segment_intersection<double>( { p1_orig[p1_i], p1_orig[p1_i + 1] }, { p2_orig[p2_i], p2_orig[p2_i + 1] } );
			if ( intersection_time.first < 0.0 ) {
				assert( intersection_time.second < 0.0 );
				continue; 
			}
			intersection_times_p1i.push_back( { p1_i, p2_i, intersection_time.first, intersection_time.second } );
		}
		
		//sort by time
		intersection_times_p1i = fplus::sort_by( []( const intersection& i1, const intersection& i2 ) -> bool {
			return (i1.t1 == i2.t1)? i1.p2_index < i2.p2_index : (i1.t1 < i2.t1); ///stable sortation with respect to p2_index in case of equal cutting times
		}, intersection_times_p1i );
		
		//insert all intersections into the luts and the refinement of p1_orig
		double last_time = -1.0;
		
		for ( const auto& i : intersection_times_p1i ) { ///insert vertices of p1_orig, which are also intersections, into the intersection_indices
			if ( i.t1 == 0.0 || i.t1 == last_time ) { ///do not insert intersections which have already been inserted as they appeared with time 1 in the segment before
				last_time = i.t1;
			}
			else if( i.t1 != 1.0) {
				last_time = i.t1;
				p1.push_back( point_on_segment( p1_orig[p1_i], p1_orig[p1_i + 1], i.t1 ) );
				intersection_indices.push_back( { p1.size() - 1, -1 } );
			}
			else {
				intersection_indices.push_back( { p1.size() , -1 } );
				last_time = i.t1;
			}
			segment_intersection_lut.push_back( { i, intersection_indices.size() -1 } );
		}
	}
	p1.push_back( p1_orig.back() );


	//compute the refinement of p2
	p2 = build_p2( p1, p2_orig );
	

	//Sanity checks
	assert( segment_intersection_lut.front().first.t1 == 0.0 && segment_intersection_lut.front().first.t2 == 0.0 ); ///first intersection
	assert( segment_intersection_lut.back().first.t1 == 1.0 && segment_intersection_lut.back().first.t2 == 1.0); ///last intersection
	assert( segment_intersection_lut.front().second == 0 && segment_intersection_lut.back().second == intersection_indices.size()-1 ); ///indices of first and last intersection
	
	return{ p1, p2 };
}


//Determines if a given polygon meets the requirements of the Chambers/Wang algorithm
bool cw_admissible_polygon( const polygon& p) {
	//Keine degenerierten Kanten
	for ( std::size_t p1 = 0; p1 <= p.size() - 2; p1++ ) {
		if ( p[p1] == p[p1 + 1] ) {
			return false;
		}
		point vec = p[p1] - p[p1 + 1];
		if ( norm( vec) < 0.0000000000001 ) {
			return false;
		}
	}

	//Keine Selbstschnitte
	if ( p.size() == 2 ) {
		return true;
	}
	for ( std::size_t p1 = 0; p1 <= p.size() - 3; p1++ ) {
		if ( segment_intersection<double>( { p[p1], p[p1 + 1] }, { p[p1 + 1], p[p1 + 2] } ) != std::pair<double, double>{1.0, 0.0} ) {
			//Aufeinanderfolgende schneiden sich nicht nur in End/Anfangspkt
			return false;
		}
		for ( std::size_t p2 = p1 + 2; p2 <= p.size() - 2; p2++ ) {
			if ( segment_intersection<double>( { p[p1], p[p1 + 1] }, { p[p2], p[p2 + 1] } ) != std::pair<double, double>{ -1.0, -1.0 } ) {
				//Nicht aufeinanderfolgende schneiden sich nicht
				return false;
			}
		}
	}
	return true;
}


/**Computes the graph defined by the sub-polygons p1[intersect1, intersect2], p2[intersect1, intersect2] of p1 and p2.
*  @param intersect1, intersect2 : indices of the intersection at which the sub-polygons start and end, respectively.
* (Note: intersect1 and intersect2 are indices of intersections, not vertices of p1 and p2. Must be translated to the latter using the intersection_indices lut.
* 
*/
embedded_graph create_embedded_graph( const polygon& p1, const polygon& p2, const std::size_t intersect1, const std::size_t intersect2 ) {
	assert( fplus::is_sorted_by( [](std::pair<std::size_t, std::size_t> i1, std::pair<std::size_t, std::size_t> i2) -> bool  {
		return i1.first < i2.first;
	}, intersection_indices ) );

	auto intersection_indices_p1 = fplus::drop_if( [&]( std::pair<std::size_t, std::size_t> intersection_idcs) -> bool {
		auto first_intersect = intersection_indices[intersect1];
		auto last_intersect = intersection_indices[intersect2];
		if ( intersection_idcs.first < first_intersect.first ||
			 intersection_idcs.first > last_intersect.first ||
			 intersection_idcs.second < first_intersect.second ||
			 intersection_idcs.second > last_intersect.second )
		{
			return true;
		}
		return false;
	}, intersection_indices );

	//the graph
	adjacency_lists_t adjacency_lists;
	vertex_coordinates_t vertex_coordinates;
	//lookup which (graph-) indices the intersections have
	std::map<std::pair<std::size_t, std::size_t>, std::size_t> intersection_graph_indices_lut;


	//p1 : just add its vertices beteween intersect1, intersect2 one after another into the graph.
	//Thus the graph index of vertex number i will be i - p1_vertex_index(intersect1)
	std::size_t p1_arc_length = intersection_indices_p1.back().first - intersection_indices_p1.front().first + 1;
	std::size_t next_p1_intersect = 0;
	
	//for each p1 vertex
	for ( std::size_t p1_i = intersection_indices_p1.front().first; p1_i <= intersection_indices_p1.back().first; p1_i++ ) {
		//note coordinates
		vertex_coordinates.push_back( p1[p1_i] );

		//if current vertex is an intersection, add its graph index to the lut
		if ( p1_i == intersection_indices_p1[next_p1_intersect].first ) {
			intersection_graph_indices_lut[intersection_indices_p1[next_p1_intersect] ] = adjacency_lists.size();
			next_p1_intersect++;
		}
		//insert predecessor
		if ( p1_i > intersection_indices_p1.front().first ) {
			adjacency_lists.push_back( { adjacency_lists.size() - 1 } );
		}
		else {
			adjacency_lists.push_back( { } ); ///first vertex has no predecessor
		}
		//insert successor
		if ( p1_i < intersection_indices_p1.back().first ) { ///last vertex has no successor
			adjacency_lists.back().push_back( adjacency_lists.size() );
		}
	}
	

	//p2 : add its vertices beteween intersect1, intersect2 into the graph.
	//Some of them (the intersections with p1) have already been added.
	//Thus, if a vertex, its predecessor or its successor is an intersection, we have to look up its graph index in the intersection_graph_indices_lut.
	auto intersection_indices_p2 = fplus::sort_on( fplus::snd<std::size_t, std::size_t>, intersection_indices_p1 );
	
	//determine the indices of intersect1 and intersect2 along p2
/*	auto intersect1_p2_maybe = fplus::find_first_idx( intersection_indices[intersect1], intersection_indices_p2 );
	assert( intersect1_p2_maybe.is_just() );
	std::size_t first_p2_intersect = intersect1_p2_maybe.unsafe_get_just();
	auto intersect2_p2_maybe = fplus::find_first_idx( intersection_indices[intersect2], intersection_indices_p2 );
	assert( intersect2_p2_maybe.is_just() );
	std::size_t last_p2_intersect = intersect2_p2_maybe.unsafe_get_just();
	*/
	std::size_t next_p2_intersect = 1; ///next intersection which will be encountered

	//tracker
	std::size_t last_p2_graph_index( 0 );
	std::size_t next_p2_graph_index( 0 );

	//for each p2 vertex
	for ( std::size_t p2_i = intersection_indices_p2.front().second; p2_i <= intersection_indices_p2.back().second; p2_i++ ) {
		//determine its graph index
		std::size_t current_graph_index( next_p2_graph_index );	
		//if its not an intersection
		if( current_graph_index == adjacency_lists.size() ){
			//note coordinates
			vertex_coordinates.push_back( p2[p2_i] );
			//Insert adjacency list into which predecessor and successor are inserted
			adjacency_lists.push_back( {} );
		}
		//insert predecessor
		if ( p2_i > intersection_indices_p2.front().second ) { ///first vertex has no predecessor
			adjacency_lists[current_graph_index].push_back( last_p2_graph_index );
		}
		//insert successor
		if ( p2_i < intersection_indices_p2.back().second ) {	///last vertex has no successor
			//successor is an intersection -> look it up
			if ( intersection_indices_p2[next_p2_intersect].second == p2_i + 1 ) {
				std::size_t successor_index = fplus::get_from_map( intersection_graph_indices_lut, intersection_indices_p2[next_p2_intersect] ).unsafe_get_just();
				adjacency_lists[current_graph_index].push_back( successor_index );
				//note next graph index and next intersection
				next_p2_graph_index = successor_index;
				next_p2_intersect++;
			}
			//otherwise
			else {
				adjacency_lists[current_graph_index].push_back( adjacency_lists.size() );
				//note next graph index
				next_p2_graph_index = adjacency_lists.size();
			}
		}
		//note last graph index
		last_p2_graph_index = current_graph_index;
	}


	//sanity checks
	std::size_t p2_arc_length = intersection_indices_p2.back().second - intersection_indices_p2.front().second + 1;
	std::size_t n_vertices = p1_arc_length + p2_arc_length - intersection_indices_p1.size();

	assert( n_vertices == vertex_coordinates.size() );
	assert( n_vertices == adjacency_lists.size() );


	return{ adjacency_lists, vertex_coordinates };
}


//Concats the sub-polygons between intersect1 and intersect2, to a closed polygon p1 ° rev(p2)
polygon p1_rev_p2( const polygon& p1, const polygon& p2, const std::size_t intersect1, const std::size_t intersect2 ) {
	polygon q;
	for ( std::size_t v = intersection_indices[intersect1].first; v <= intersection_indices[intersect2].first; v++ ) {
		q.push_back( p1[v] );
	}
	for ( std::size_t v = intersection_indices[intersect2].second ; v >= intersection_indices[intersect1].second + 1; v-- ) {
		q.push_back( p2[v-1] );
	}
	return q;
}


/**Calculates the winding number of a closed curve with respect to x. 
* x must not be on the curve or too close to an arc, otherwise a nothing will be returned
*/
fplus::maybe<int> winding_number( const polygon& curve, const point& x ) {
	assert( curve.size() >= 2 );
	assert( curve.front() == curve.back() );

	double angle = 0.0;
	point xv1_vec = curve[0] - x;
	for ( std::size_t i = 0; i <= curve.size()-2; i++ ) {
		point xv2_vec = curve[i + 1] - x;

		if ( geom::norm( xv1_vec ) < 0.1 || geom::norm( xv2_vec ) < 0.1 ) { ///here to prevent errors caused by floating point inaccuracy. x should be chosen such that vertices are at a certain distance (here 0.01).
			return{};
		}
		double segment_angle = geom2d::angle( xv1_vec, xv2_vec );
		if ( std::abs(segment_angle) >= 179.9 ) { ///here to prevent errors caused by floating point inaccuracy. x should be chosen such that all angles are measured qualitatively correct. That is, segments are far away from x such that it's clear on which side of x they pass. TODO: Prove this is (un)necessary
			return{};
		}
		angle += segment_angle;

		xv1_vec = xv2_vec;
	}

	//translate angle to a number of winds
	assert( std::abs(angle) < std::numeric_limits<int>::max() - 2 );
	int winds_deg = fplus::round<int>( angle );
	assert( winds_deg % 360 == 0 );
	int winds = winds_deg / 360;

	return{ winds };
}


//Computes a point inside a non self intersecting polygon
point find_point_in_polygon( const polygon & p ) {
	//try centroid
	point centroid = geom2d::centroid( p );
	auto windingnumber_maybe = winding_number( p, centroid );
	int windingnumber;
	if ( windingnumber_maybe.is_nothing() ) {
		windingnumber = 0;
	}
	else {
		windingnumber = windingnumber_maybe.unsafe_get_just();
	}
	assert( std::abs( windingnumber ) <= 1 );

	if ( windingnumber != 0 ) {
		return centroid;
	}

	//for each segment
	for ( std::size_t v = 0; v <= p.size() - 2; v++ ) {
		point segment_vec = p[v + 1] - p[v];
		point segment_center = p[v] + 0.5*(segment_vec);
		point normal = geom2d::normal( segment_vec );
		
		double delta = 1.0;
		for ( std::size_t i = 0 ; i < 5; i++, delta /= 2 ) {
			//try point left of the segment
			point x = segment_center + delta*normal;
			windingnumber_maybe = winding_number( p, x );
			if ( windingnumber_maybe.is_nothing() ) {
				windingnumber = 0;
			}
			else {
				windingnumber = windingnumber_maybe.unsafe_get_just();
			}
			assert( std::abs( windingnumber ) <= 1 );
			if ( windingnumber != 0 ) {
				return x;
			}
			//try point right of the segment
			point y = segment_center - delta*normal;
			windingnumber_maybe = winding_number( p, y );
			if ( windingnumber_maybe.is_nothing() ) {
				windingnumber = 0;
			}
			else {
				windingnumber = windingnumber_maybe.unsafe_get_just();
			}
			assert( std::abs( windingnumber ) <= 1 );
			if ( windingnumber != 0 ) {
				return y;
			}
		}
	}

	assert( false ); ///found no point inside the polygon
	return{ -1,-1 };
}

typedef std::pair<polygon, int> windingnumber;
typedef std::vector<windingnumber> windingnumbers; ///Optional: TODO: Statt dem Polygon kann man auch direkt die jeweilige Flaeche berechnen und speichern? Oder direkt die TW?
/**Calculates the winding numbers of all faces of the graph defined by the union of the sub-polygons (p1[intersect1,intersect2] v p2[intersect1,intersect2]), which share their intersections as vertices
*  @param intersect1, intersect2 : indices of the intersection at which the sub-polygons start and end, respectively.
* (Note: intersect1 and intersect2 are indices of intersections, not vertices of p1 and p2. Must be translated to the latter using the intersection_indices lut.
*/
windingnumbers calculate_winding_numbers( const polygon& p1, const polygon& p2, const std::size_t intersect1, const std::size_t intersect2 ) {
	polygon p1_revp2 = p1_rev_p2( p1, p2, intersect1, intersect2 );
	
	//create graph from sub-polygons
	embedded_graph g = create_embedded_graph( p1, p2, intersect1, intersect2 );
	
	//enumerate its faces
	face_enumeration faces = enumerate_embedding_faces( g );

	//transform graph indices into polygons
	auto face_polygons = fplus::transform( [&g]( adjacency_list_t face_indices ) -> polygon {
		polygon face;
		for ( const auto& i : face_indices ) {
			face.push_back( g.vertex_coordinates[i] );
		}
		face.push_back( face.front() );
		return face;
	}, faces );

	//delete complement of the graph from the face list, which is either the only face with self intersections, or the one with maximal area, because its bouning polygon is the outer contour of the graph
	std::size_t complement_face_index = -1;
	double max_area = -1.0;
	for ( std::size_t i = 0; i < face_polygons.size(); i++ ) {
		polygon vertices = face_polygons[i];
		vertices.pop_back();
		if ( !fplus::all_unique( vertices ) ) {
			complement_face_index = i;
			break;
		}
		double area = geom2d::area( face_polygons[i] );
		if ( area > max_area ) {
			max_area = area;
			complement_face_index = i;
		}
	}
	face_polygons = fplus::drop_idxs( std::vector<std::size_t>{ complement_face_index }, face_polygons );

	//calculate their windingnumbers
	windingnumbers winding_numbers = fplus::transform( [&p1_revp2]( polygon face ) -> std::pair<polygon, int> {
		//point in that polygon
		point x = find_point_in_polygon( face );

		auto windingnumber_maybe = winding_number( p1_revp2, x );
		assert( windingnumber_maybe.is_just() );
		return{ face, windingnumber_maybe.unsafe_get_just() };
	}, face_polygons );

	
	return winding_numbers;
}


//Checks if all winding numbers are ">= 0" or "<= 0"
bool consistent( const windingnumbers& wn ) {
	assert( !wn.empty() );
	int s = 0; ///Keeps track of the sign of the winding numbers
	for ( const auto& w : wn ) {
		if ( w.second == 0 ) { continue; } ///Fits every sign

		int prod = geom::sign( s * w.second );
		
		if ( prod == 0 && w.second != 0 ) { ///Sign was 0, now changes to w.second
			s = geom::sign( w.second );
			continue;
		}
		else if ( prod < 0 ) { ///Only case in which the product is x1*x2 < 0 occurs if xi < 0 and xj > 0
			return  false;
		}
	}
	return true;
}


//Computes the total winding number from a set of faces and their single winding numbers, which is the area of the faces multiplied with their winding number
double total_winding_number( const windingnumbers& winding_numbers ) {
	std::vector<double> multiplied_areas = fplus::transform( []( const std::pair<polygon, int>& wn ) -> double {
		double area = geom2d::area( wn.first );
		return area * std::abs( wn.second );
	}, winding_numbers );

	return fplus::sum( multiplied_areas );
}


//Determines if an intersection is an anchor point 
bool valid_intersection( std::size_t intersect1, std::size_t intersect2) {
	assert( intersection_indices[intersect1] != intersection_indices[intersect2] );
	
	if( (intersection_indices[intersect1].first >= intersection_indices[intersect2].first) ||
		(intersection_indices[intersect1].second >= intersection_indices[intersect2].second) )
	{ return false; }
	return true;
}

double global_min = DBL_MAX;

//Recursive function, which computes the total winding number (TW) of Q := p1°rev(p2) by splitting Q at its anchor points and summing up the TWs of the parts
double decompose_at_anchor_points( const polygon& p1, const polygon& p2, std::size_t intersect1, std::size_t intersect2, double current_homotopy_area) {
	assert(( 0 <= intersect1) && (intersect1 < intersect2) && (intersect2 < intersection_indices.size()) );
	if ( current_homotopy_area > global_min ) { return 0.0; }

	auto winding_numbers = calculate_winding_numbers( p1, p2, intersect1, intersect2 );
	if ( intersect1 == intersect2 - 1 ) { assert( consistent( winding_numbers ) ); }
	if ( consistent( winding_numbers ) ) {
		assert( intersect2 == intersection_indices.size() - 1 );
		double area = total_winding_number( winding_numbers );
		global_min = std::min( global_min, current_homotopy_area + area );
		return area;
	}
	else {
		double min_area = DBL_MAX;
		for ( std::size_t intersect = intersect1 + 1; intersect < intersect2; intersect++ ) {
			if ( !valid_intersection( intersect1, intersect) ) { continue; }
			auto first_half_winding_numbers = calculate_winding_numbers( p1, p2, intersect1, intersect );
			if ( !consistent( first_half_winding_numbers ) ) { continue; }
			double first_half = total_winding_number( first_half_winding_numbers );
			double second_half = decompose_at_anchor_points( p1, p2, intersect, intersect2, current_homotopy_area + first_half );
			double area_if_intersect_is_anchor = first_half + second_half;
			min_area = std::min( min_area, area_if_intersect_is_anchor );
		}
		return min_area;
	}
}


//Calculates the minimum area homotopy as in Chambers/Wang (see http://arxiv.org/abs/1303.7427)
fplus::result<double, std::string> min_ht_area::calculate_min_homotopy_area( const polygon& p1, const polygon& p2, bool check_polygon ) {
	//check boundary conditions
	if ( p1.size() < 2 || p2.size() < 2 ) {
		return fplus::error<double, std::string>( std::string( "Polygons must at least be of length 2. Given lengths: " +
															   std::to_string( p1.size() ) +
															   ", " + std::to_string( p2.size() )
															   ) );
	}
	if ( p1.front() != p2.front() || p1.back() != p2.back() ) {
		return fplus::error<double, std::string>(std::string( "First and last point of the polygons p1 and p2 must be equal! "
															  "(In order to be homotopic modulo {0,1}.)"
															  ));
	}
	//check all necessary polygon restraints
	if ( check_polygon && (!cw_admissible_polygon( p1 ) || !cw_admissible_polygon( p2 ) ) ){
		return fplus::error<double, std::string>( std::string( "One of the given polygons does not fulfill the constraints of the Chambers/Wang algorithm!") ); ///Optional: TODO: Be more specific. Which constraints exist, which one has been hurt?
	}

	//initialize luts and introduce vertices into the polygons at their mutual intersections
	auto splitted_polygons = initialize_intersection_lut( p1, p2 );
	polygon p1_refined = splitted_polygons.first;
	polygon p2_refined = splitted_polygons.second;

	//compute total winding number according to Chambers/Wang
	global_min = DBL_MAX;
	auto result = decompose_at_anchor_points( p1_refined, p2_refined, 0, intersection_indices.size()-1, 0.0 );

	return fplus::ok<double, std::string>( result );
}
