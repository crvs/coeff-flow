#pragma once

#include <vector>
#include <algorithm>

#include "scomplex/Simplicial_Complex.h"
#include "scomplex/graph_utils.hpp"

// rtree stuff
#include "boost/geometry/index/rtree.hpp"
#include "boost/geometry/geometries/point.hpp"
#include "boost/geometry/geometries/box.hpp"

//#include "boost/geometry/index/parameters.hpp"
//#include "boost/geometry.hpp"

// using namespace boost;

typedef std::vector<double> point_t;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template <size_t Dimension>
using point_template = bg::model::point<double, Dimension, bg::cs::cartesian>;

template <size_t Dimension>
void get_r_tree_from_points<Dimension>(std::list<point_t>& points) {
    typedef point_template<Dimension> point;
    typedef bg::model::box<point> box;
    typedef std::pair<box, unsigned> value;

    // default rtree constructor
    bgi::rtree<value, bgi::quadratic<16>> r_tree;

    /*
        point_t,                  //
        bgi::quadratic,           //
        bgi::indexable<point_t>,  //
        bgi::equal_to<point_t>,   //
        std::allocator<point_t>   //
        > RTree(points.begin(), points.end());
        */
}
