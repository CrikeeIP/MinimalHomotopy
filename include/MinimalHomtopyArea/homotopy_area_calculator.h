// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/Spielwiese
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#pragma once


#include "..\Geometry\include\geometry\geometry2d.h"

#include <fplus/fplus.h>

#include <vector>


namespace min_ht_area {

	typedef geom2d::Vec2D<double> point;

	fplus::result<double, std::string> calculate_min_homotopy_area( const geom2d::polygon& p1, const geom2d::polygon& p2, bool check_polygon = true );

} //namespace min_ht_area