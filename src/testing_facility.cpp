
#include <vector>
#include <iostream>

#include <yaml-cpp/yaml.h>
#include <fstream>
#include <algorithm>

#include "types.hpp"
#include "chains.hpp"
#include "SimpComplex.hpp"
#include "qhull_parsing.hpp"
#include "path_snapper.hpp"
#include "BoundingChain.hpp"
#include "coeff_flow.hpp"
#include "plywriter.hpp"

#include "meshreader.hpp"

#include <time.h>

using namespace std;
using namespace Gsimp;

void call_single_cycle_test(string, vector< size_t >, vector< Point >, size_t,
                            bool);

void call_single_cycle_test(string, vector< Point >, size_t, bool);

void call_single_cycle_test(string, vector< Point >, size_t);

void call_single_cycle_test(string, vector< size_t >, size_t, bool);

void call_single_cycle_test(string, vector< size_t >, size_t);

int main(int argc, char* argv[]) {
    YAML::Node input = YAML::LoadFile(argv[1]);

    string paths_file;
    string complex_file;
    vector< Point  > cycle_points;
    vector< size_t > cycle_index_points;
    bool in_plane;

    if (input["single_cycle_test"]) {
        complex_file =
            input["single_cycle_test"]["complex_file"].as< string >();
        in_plane = input["single_cycle_test"]["in_plane"].as< bool >();

        size_t null_face;
        if (!in_plane)
            null_face = input["single_cycle_test"]["null_face"].as< size_t >();
        else
            null_face = 0;
        if (input["single_cycle_test"]["cycle_index_points"]) {
            cycle_index_points =
                input["single_cycle_test"]["cycle_index_points"]
                    .as< vector< size_t > >();
            call_single_cycle_test(complex_file, cycle_index_points, null_face,
                                   in_plane);
        } else {
            cycle_points = input["single_cycle_test"]["cycle_points"]
                               .as< vector< Point > >();
            cycle_points.push_back(cycle_points[0]);
            call_single_cycle_test(complex_file, cycle_points, null_face,
                                   in_plane);
        }
    }
}

void call_single_cycle_test(string complex_file, vector< Point > cycle_points,
                            size_t null_face, bool in_plane) {
    call_single_cycle_test(complex_file, {}, cycle_points, null_face, in_plane);
}

void call_single_cycle_test(string complex_file, vector< Point > cycle_points,
                            size_t null_face) {
    call_single_cycle_test(complex_file, {}, cycle_points, null_face, false);
}

void call_single_cycle_test(string complex_file, vector< size_t > cycle_indices,
                            size_t null_face, bool in_plane) {
    call_single_cycle_test(complex_file, cycle_indices, {}, null_face,
                           in_plane);
}

void call_single_cycle_test(string complex_file, vector< size_t > cycle_indices,
                            size_t null_face) {
    call_single_cycle_test(complex_file, cycle_indices, {}, null_face, false);
}

//
// mother of all single_cycle_tests
//
void call_single_cycle_test(string complex_file, vector< size_t > cycle_indices,
                            vector< Point > cycle_points, size_t null_face,
                            bool in_plane) {
    cout << "performing the single cycle test with parameters:\n";
    cout << "file containing the mesh: " << complex_file << "\n";
    cout << "points that will form a cycle:\n";
    for (auto p : cycle_points) {
        cout << "  * ";
        for (auto c : p) cout << c << " ";
        cout << "\n";
    }

    cout << "\b\b  \n";
    if (!in_plane) cout << "index of face to be null: " << null_face << "\n\n";

    vector< point_t > points_v;
    vector< Cell > cells_v;

    clock_t t0, t1;
    t0 = clock();
    // reading mesh and shenanigans
    SimpleMesh mesh;


    string meshtype;
    if (meshtype == "ply") {
        mesh = readMesh(complex_file);
    } else {
        throw std::runtime_error("can only handle ply files");
    }

    points_v = mesh.vertices;
    cells_v = mesh.faces;


    // end PCL shenanigans
    t1 = clock();
    cout << "mesh has " << points_v.size() << " vertices and " << cells_v.size()
         << " faces\n";
    cout << "read mesh in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    // watch out for 0 or 1 indexing
    bool zero_index = true;
    for (Cell cell : cells_v) {
        for (size_t vert : cell) {
            if (vert >= points_v.size()) {
                zero_index = false;
                break;
            }
        }
        if (!zero_index) break;
    }

    if (!zero_index) {
        for (int i = 0; i < cells_v.size(); ++i) {
            for (int j = 0; j < cells_v[i].size(); ++j) cells_v[i][j] -= 1;
        }
    }

    t0 = clock();
    shared_ptr< SimpComplex > s_comp(new SimpComplex(points_v, cells_v));
    t1 = clock();
    cout << "created complex in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    cout << "complex comosition:\n";
    cout << "    number of faces: " << s_comp->get_level_size(2) << "\n";
    cout << "    number of edges: " << s_comp->get_level_size(1) << "\n";
    cout << "    number of vertices: " << s_comp->get_level_size(0) << "\n";

    // re sort the triangles according to level 2
    {
        typedef pair< size_t, Cell > sortable;
        vector< sortable > pairing;
        for (auto tri : cells_v) {
            // get index of the triangle
            size_t ind;
            ind = s_comp->cell_to_index(tri);
            pairing.push_back({ind, tri});
        }
        sort(pairing.begin(), pairing.end(),
             [](const sortable& x, const sortable& y) {
                 return x.first < y.first;
             });
        cells_v.clear();
        for (auto el : pairing) {
            cells_v.push_back(el.second);
        }
    }

    // end sorting
    t0 = clock();
    shared_ptr< path_snapper > p_snap(new path_snapper(s_comp));
    t1 = clock();
    cout << "created snapper in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    vector< size_t > snapped;
    if (cycle_indices.size() == 0)
        snapped = p_snap->snap_path_to_indices(cycle_points);
    else {
        snapped = cycle_indices;
    }

    t1 = clock();
    cout << "snapped path in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";
    cout << "computed path contains " << snapped.size() << " points\n";

    t0 = clock();
    chain cycle = p_snap->index_sequence_to_dense_chain(snapped);
    t1 = clock();
    cout << "produced chain from path in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    BoundingChain ch_calc(s_comp);
    t1 = clock();
    cout << "calculated boundary matrices in " << t1 - t0 << " clock cycles\
        " << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    chain cycle2(cycle);
    cycle2.to_sparse();
    t0 = clock();
    chain b_chain_0 = ch_calc.getBoundingChain(cycle2);
    t1 = clock();
    cout << "calculated bounding chain (using Eigen) in " << t1 - t0 << "\
        clock cycles " << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    vector< vector< int > > colors;

    for (size_t i = 0; i < b_chain_0.get_size(); i++) {
        if (abs(b_chain_0[i]) > 10e-3)
            colors.push_back({255, 0, 0});
        else
            colors.push_back({255, 255, 255});
    }

    ofstream my_ply;
    my_ply.open("my_ply.ply");

    vector< vector< size_t > > edges{};
    vector< vector< int > > edge_colors{};

    for (size_t i = 0; i < cycle.get_size(); i++) {
        if (abs(cycle[i]) > 10e-3) edges.push_back(s_comp->index_to_cell(1, i));
        edge_colors.push_back({0, 0, 255});
    }

    make_ply(my_ply, s_comp->points(), cells_v, colors, edges, edge_colors);
    my_ply.close();

    Cell null_cell;
    if (!in_plane) {
        null_cell = s_comp->index_to_cell(s_comp->dimension(), null_face);
        cout << "null cell: ";
        copy(null_cell.begin(), null_cell.end(),
             ostream_iterator< size_t >(cout, " "));
        cout << "\n";
    }

    t0 = clock();
    s_comp->calculate_hasse();
    t1 = clock();
    cout << "calculated the hasse diagram in " << t1 - t0 << " clock cycles "
         << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";
    t0 = clock();

    chain b_chain_1 = s_comp->new_dense_chain(2);
    if (!in_plane)
        b_chain_1 = coeff_flow(*s_comp, cycle, null_cell, 0);
    else
        b_chain_1 = coeff_flow_embedded(*s_comp, cycle);

    t1 = clock();
    cout << "calculated bounding chain using coeff_flow in " << t1 - t0
         << " clock cycles " << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    colors = {};
    for (size_t i = 0; i < b_chain_1.get_size(); i++) {
        if (abs(b_chain_1[i]) > 10e-3)
            colors.push_back({255, 0, 0});
        else
            colors.push_back({255, 255, 255});
    }

    ofstream my_ply2;
    my_ply2.open("my_ply2.ply");

    edges = {};
    edge_colors = {};
    for (size_t i = 0; i < cycle.get_size(); i++) {
        if (abs(cycle[i]) > 10e-3) edges.push_back(s_comp->index_to_cell(1, i));
        edge_colors.push_back({0, 0, 255});
    }

    make_ply(my_ply2, s_comp->points(), cells_v, colors, edges, edge_colors);
    my_ply2.close();
}
