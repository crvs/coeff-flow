#include <pcl/io/obj_io.h>
#include <pcl/io/ply_io.h>

#include <iostream>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <fstream>

#include <algorithm>
#include <scomplex/chain_calc.hpp>
#include <scomplex/coeff_flow.hpp>
#include <scomplex/path_snapper.hpp>
#include <scomplex/plywriter.hpp>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

#include <time.h>

void call_single_cycle_test(std::string, std::vector< size_t >,
                            std::vector< std::vector< double > >, size_t, bool);

void call_single_cycle_test(std::string, std::vector< std::vector< double > >,
                            size_t, bool);

void call_single_cycle_test(std::string, std::vector< std::vector< double > >,
                            size_t);

void call_single_cycle_test(std::string, std::vector< size_t >, size_t, bool);

void call_single_cycle_test(std::string, std::vector< size_t >, size_t);

int main(int argc, char* argv[]) {
    YAML::Node input = YAML::LoadFile(argv[1]);

    std::string paths_file;
    std::string complex_file;
    std::vector< std::vector< double > > cycle_points;
    std::vector< size_t > cycle_index_points;
    bool in_plane;

    if (input["single_cycle_test"]) {
        complex_file =
            input["single_cycle_test"]["complex_file"].as< std::string >();
        in_plane = input["single_cycle_test"]["in_plane"].as< bool >();

        size_t null_face;
        if (!in_plane)
            null_face = input["single_cycle_test"]["null_face"].as< size_t >();
        else
            null_face = 0;
        if (input["single_cycle_test"]["cycle_index_points"]) {
            cycle_index_points =
                input["single_cycle_test"]["cycle_index_points"]
                    .as< std::vector< size_t > >();
            call_single_cycle_test(complex_file, cycle_index_points, null_face,
                                   in_plane);
        } else {
            cycle_points = input["single_cycle_test"]["cycle_points"]
                               .as< std::vector< std::vector< double > > >();
            cycle_points.push_back(cycle_points[0]);
            call_single_cycle_test(complex_file, cycle_points, null_face,
                                   in_plane);
        }
    }
}

void call_single_cycle_test(std::string complex_file,
                            std::vector< std::vector< double > > cycle_points,
                            size_t null_face, bool in_plane) {
    call_single_cycle_test(complex_file, {}, cycle_points, null_face, in_plane);
}

void call_single_cycle_test(std::string complex_file,
                            std::vector< std::vector< double > > cycle_points,
                            size_t null_face) {
    call_single_cycle_test(complex_file, {}, cycle_points, null_face, false);
}

void call_single_cycle_test(std::string complex_file,
                            std::vector< size_t > cycle_indices,
                            size_t null_face, bool in_plane) {
    call_single_cycle_test(complex_file, cycle_indices, {}, null_face,
                           in_plane);
}

void call_single_cycle_test(std::string complex_file,
                            std::vector< size_t > cycle_indices,
                            size_t null_face) {
    call_single_cycle_test(complex_file, cycle_indices, {}, null_face, false);
}

//
// mother of all single_cycle_tests
//
void call_single_cycle_test(std::string complex_file,
                            std::vector< size_t > cycle_indices,
                            std::vector< std::vector< double > > cycle_points,
                            size_t null_face, bool in_plane) {
    std::cout << "performing the single cycle test with parameters:\n";
    std::cout << "file containing the mesh: " << complex_file << "\n";
    std::cout << "points that will form a cycle:\n";
    for (auto p : cycle_points) {
        std::cout << "  * ";
        for (auto c : p) std::cout << c << " ";
        std::cout << "\n";
    }

    std::cout << "\b\b  \n";
    if (!in_plane)
        std::cout << "index of face to be null: " << null_face << "\n\n";

    std::vector< point_t > points_v;
    std::vector< cell_t > cells_v;

    clock_t t0, t1;
    t0 = clock();
    // PCL shenanigans

    pcl::PolygonMesh::Ptr mesh(new pcl::PolygonMesh{});

    std::ifstream file(complex_file);
    std::string meshtype;
    std::getline(file, meshtype);
    if (meshtype == "ply") {
        pcl::PLYReader Reader;
        Reader.read(complex_file, *mesh);
    } else {
        pcl::OBJReader Reader;
        Reader.read(complex_file, *mesh);
    }

    pcl::PointCloud< pcl::PointXYZ >::Ptr cloud(
        new pcl::PointCloud< pcl::PointXYZ >{});
    pcl::fromPCLPointCloud2(mesh->cloud, *cloud);

    for (auto poly : mesh->polygons) {
        cells_v.push_back(
            {poly.vertices[0], poly.vertices[1], poly.vertices[2]});
    }

    for (auto pt : cloud->points) {
        points_v.push_back({pt.x, pt.y, pt.z});
    }

    // end PCL shenanigans
    t1 = clock();
    std::cout << "mesh has " << points_v.size() << " vertices and "
              << cells_v.size() << " faces\n";
    std::cout << "read mesh in " << t1 - t0 << " clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    // watch out for 0 or 1 indexing
    bool zero_index = true;
    for (cell_t cell : cells_v) {
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
    std::shared_ptr< gsimp::simplicial_complex > s_comp =
        std::make_shared< gsimp::simplicial_complex >(points_v, cells_v);
    t1 = clock();
    std::cout << "created complex in " << t1 - t0 << " clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    std::cout << "complex comosition:\n";
    std::cout << "    number of faces: " << s_comp->get_level_size(2) << "\n";
    std::cout << "    number of edges: " << s_comp->get_level_size(1) << "\n";
    std::cout << "    number of vertices: " << s_comp->get_level_size(0)
              << "\n";

    // re sort the triangles according to level 2
    {
        typedef std::pair< size_t, cell_t > sortable;
        std::vector< sortable > pairing;
        for (auto tri : cells_v) {
            // get index of the triangle
            size_t ind;
            ind = s_comp->cell_to_index(tri);
            pairing.push_back({ind, tri});
        }
        std::sort(pairing.begin(), pairing.end(),
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
    std::shared_ptr< gsimp::path_snapper > p_snap =
        std::make_shared< gsimp::path_snapper >(s_comp);
    t1 = clock();
    std::cout << "created snapper in " << t1 - t0 << " clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    std::vector< size_t > snapped;
    if (cycle_indices.size() == 0)
        snapped = p_snap->snap_path_to_indices(cycle_points);
    else {
        snapped = cycle_indices;
    }

    t1 = clock();
    std::cout << "snapped path in " << t1 - t0 << " clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";
    std::cout << "computed path contains " << snapped.size() << " points\n";

    t0 = clock();
    gsimp::chain_v cycle = p_snap->index_sequence_to_v_chain(snapped);
    t1 = clock();
    std::cout << "produced chain from path in " << t1 - t0 << " clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    gsimp::bounding_chain ch_calc(s_comp);
    t1 = clock();
    std::cout << "calculated boundary matrices in " << t1 - t0
              << " clock cycles\
        " << float(t1 - t0) / CLOCKS_PER_SEC
              << " seconds\n";

    gsimp::chain_t cycle2 = p_snap->index_sequence_to_chain(snapped);
    t0 = clock();
    gsimp::chain_t b_chain_0 = ch_calc.get_bounding_chain(cycle2);
    t1 = clock();
    std::cout << "calculated bounding chain (using Eigen) in " << t1 - t0 << "\
        clock cycles "
              << float(t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

    std::vector< std::vector< int > > colors;

    for (size_t i = 0; i < gsimp::chain_size(b_chain_0); i++) {
        if (std::abs(gsimp::chain_val(b_chain_0, i)) > 10e-3)
            colors.push_back({255, 0, 0});
        else
            colors.push_back({255, 255, 255});
    }

    std::ofstream my_ply;
    my_ply.open("my_ply.ply");

    std::vector< std::vector< size_t > > edges{};
    std::vector< std::vector< int > > edge_colors{};

    for (size_t i = 0; i < gsimp::chain_size(cycle); i++) {
        if (std::abs(gsimp::chain_val(cycle, i)) > 10e-3)
            edges.push_back(s_comp->index_to_cell(1, i));
        edge_colors.push_back({0, 0, 255});
    }

    make_ply(my_ply, s_comp->get_points(), cells_v, colors, edges, edge_colors);
    my_ply.close();

    gsimp::cell_t null_cell;
    if (!in_plane) {
        null_cell = s_comp->index_to_cell(s_comp->dimension(), null_face);
        std::cout << "null cell: ";
        std::copy(null_cell.begin(), null_cell.end(),
                  std::ostream_iterator< size_t >(std::cout, " "));
        std::cout << "\n";
    }

    t0 = clock();
    s_comp->calculate_hasse();
    t1 = clock();
    std::cout << "calculated the hasse diagram in " << t1 - t0
              << " clock cycles " << float(t1 - t0) / CLOCKS_PER_SEC
              << " seconds\n";
    t0 = clock();

    gsimp::chain_v b_chain_1;
    if (!in_plane)
        b_chain_1 = gsimp::coeff_flow(*s_comp, cycle, null_cell, 0);
    else
        b_chain_1 = gsimp::coeff_flow_embedded(*s_comp, cycle);

    t1 = clock();
    std::cout << "calculated bounding chain using coeff_flow in " << t1 - t0
              << " clock cycles " << float(t1 - t0) / CLOCKS_PER_SEC
              << " seconds\n";

    colors = {};
    for (size_t i = 0; i < gsimp::chain_size(b_chain_1); i++) {
        if (std::abs(gsimp::chain_val(b_chain_1, i)) > 10e-3)
            colors.push_back({255, 0, 0});
        else
            colors.push_back({255, 255, 255});
    }

    std::ofstream my_ply2;
    my_ply2.open("my_ply2.ply");

    edges = {};
    edge_colors = {};
    for (size_t i = 0; i < gsimp::chain_size(cycle); i++) {
        if (std::abs(gsimp::chain_val(cycle, i)) > 10e-3)
            edges.push_back(s_comp->index_to_cell(1, i));
        edge_colors.push_back({0, 0, 255});
    }

    make_ply(my_ply2, s_comp->get_points(), cells_v, colors, edges,
             edge_colors);
    my_ply2.close();
}
