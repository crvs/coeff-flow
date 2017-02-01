#pragma once

#include <vector>
#include <algorithm>

#include <CGAL/Kd_tree.h>

// #include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/K_neighbor_search.h>
// #include <CGAL/Euclidean_distance.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

typedef std::vector<double> point_t;

typedef CGAL::Cartesian_d<double> cart;
typedef cart::Point_d point_d;
typedef CGAL::Search_traits_d<cart> traits;
typedef CGAL::Orthogonal_k_neighbor_search<traits> neighbor_search_t;
typedef neighbor_search_t::Tree tree_t;

// typedef CGAL::Kd_tree<traits> kd_tree;
// typedef CGAL::K_neighbor_search<traits> k_neighbor;
// typedef CGAL::Euclidean_distance<traits> euclidean;

// template <int Dimension>
// typedef cart::Point_d<Dimension> point_d;

void make_tree(tree_t& tree, std::list<point_t>& point_list) {
    tree.reserve(point_list.size());
    std::list<point_d> insertable_points;
    for (auto pt : point_list) {
        point_d pt_d(pt.size(), pt.begin(), pt.end());
        tree.insert(pt_d);
    }
    //    return tree;
};

point_t nearest_neighbour(tree_t& tree, point_t point) {
    point_d pt(point.size(), point.begin(), point.end());
    neighbor_search_t search(tree, pt, 1);
    point_d result_pt = search.begin()->first;
    point_t point_vec(result_pt.cartesian_begin(), result_pt.cartesian_end());
    return point_vec;
}
