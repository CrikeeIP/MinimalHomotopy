// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)

#pragma once


#include "..\Geometry\include\geometry\geometry2d.h"

#include <fplus/fplus.h>

#include <vector>


namespace min_ht_area {

	typedef geom2d::Vec2D<double> point;

	fplus::result<double, std::string> calculate_min_homotopy_area( const geom2d::polygon& p1, const geom2d::polygon& p2, bool check_polygon = true );

} //namespace min_ht_area
