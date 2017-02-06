#pragma once

#include <vector>
#include <algorithm>

#include <CGAL/Kd_tree.h>

#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <boost/tuple/tuple.hpp>

typedef std::vector<double> point_t;

// kd-tree point type
typedef CGAL::Cartesian_d<double>::Point_d point_d;  // d-dimensional point

typedef CGAL::Search_traits_d<CGAL::Cartesian_d<double>> cartesian_traits;

typedef CGAL::Search_traits_adapter<     // define the search tratis
    boost::tuple<point_d, size_t>,       // pairing d-point with its index
    CGAL::Nth_of_tuple_property_map<     // custom kernel
        0,                               // position of the points to compare
        boost::tuple<point_d, size_t>>,  // type tha will be used
    cartesian_traits>                    // traits to use
    traits;                              // name

// type of search tree
typedef CGAL::Orthogonal_k_neighbor_search<traits> neighbor_search_t;

typedef neighbor_search_t::Tree tree_t;

void make_tree(tree_t& tree, std::vector<point_t>& point_list) {
    // tree.reserve(point_list.size());
    for (size_t ind = 0; ind < point_list.size(); ++ind) {
        auto pt = point_list.at(ind);
        point_d pt_d(pt.size(), pt.begin(), pt.end());
        tree.insert(boost::make_tuple(pt_d, ind));
    }
}

// given a point snap it to the tree
point_t nearest_neighbor(tree_t& tree, point_t point) {
    point_d pt(point.size(), point.begin(), point.end());
    // "1" refers to the number of nearest neighbors to find
    neighbor_search_t search(tree, pt, 1);
    point_d result_pt;
    size_t ind;
    boost::tie<point_d, size_t>(result_pt, ind) = search.begin()->first;
    point_t point_vec(result_pt.cartesian_begin(), result_pt.cartesian_end());
    return point_vec;
}

// given a vector of points snap them to the tree
std::vector<point_t> snap_points(tree_t& tree, std::vector<point_t> points) {
    std::vector<point_t> snapped_points;
    for (auto point : points) {
        snapped_points.push_back(nearest_neighbor(tree, point));
    }
    return snapped_points;
}

size_t nearest_neighbor_index(tree_t& tree, point_t point) {
    point_d pt(point.size(), point.begin(), point.end());
    // "1" refers to the number of nearest neighbors to find
    neighbor_search_t search(tree, pt, 1);
    point_d result_pt;
    size_t ind;
    boost::tie<point_d, size_t>(result_pt, ind) = search.begin()->first;
    return ind;
}

std::vector<size_t> snap_points_to_indexes(tree_t& tree,
                                           std::vector<point_t> points) {
    std::vector<size_t> snapped_points;
    for (auto point : points) {
        snapped_points.push_back(nearest_neighbor_index(tree, point));
    }
    return snapped_points;
}
