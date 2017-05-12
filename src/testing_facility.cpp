#include <vector>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include <scomplex/types.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/qhull_parsing.hpp>
#include <scomplex/path_snapper.hpp>
#include <scomplex/chain_calc.hpp>
#include <scomplex/coeff_flow.hpp>

#include <time.h>

void call_single_cycle_test(std::string ,std::vector<std::vector<double>> ,size_t );

int main( int argc , char* argv[] ) {

    YAML::Node input = YAML::LoadFile(argv[1]);

    std::string paths_file;
    std::string complex_file;
    std::vector<std::vector<double>> cycle_points;
    size_t null_face;

    if(input["single_cycle_test"]){
        complex_file = input["single_cycle_test"]["complex_file"].as<std::string>();
        cycle_points = input["single_cycle_test"]["cycle_points"].as<std::vector<std::vector<double>>>();
        cycle_points.push_back(cycle_points[0]);
        null_face = input["single_cycle_test"]["null_face"].as<size_t>();
        call_single_cycle_test(complex_file,cycle_points,null_face);
    }

    /*
    if( input[0]["complex_file"])
        complex_file = input[0]["complex_file"].as<std::string>();

    if( input[0]["paths_file"])
        paths_file = input[0]["paths_file"].as<std::string>();

    if(input[0]["cycle_points"]) {
        cycle_points = input[0]["cycle_points"].as<std::vector<std::vector<double>>>();
        null_face = input[0]["null_face"];
    }
    */
}

void call_single_cycle_test(std::string complex_file,std::vector<std::vector<double>> cycle_points,size_t null_face){


    std::cout << "performing the single cycle test with parameters:\n";
    std::cout << "file containing the mesh: " << complex_file << "\n";
    std::cout << "points that will form a cycle: ";
    for(auto p : cycle_points) {
        for(auto c : p)
            std::cout << c << " ";
        std::cout << "_ ";
    }
    std::cout << "\b\b  \n";
    std::cout << "index of face to be null: " << null_face << "\n\n";

    std::vector<point_t> points_v;
    std::vector<cell_t> cells_v;

    clock_t t0,t1;
    t0 = clock();
    std::tie(points_v, cells_v) = parse_qhull_file(complex_file);
    t1 = clock();
    std::cout << "mesh has " << points_v.size() << " vertices and " << cells_v.size() << " faces\n";
    std::cout << "read mesh in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";

    // watch out for 0 or 1 indexing
    bool zero_index = true;
    for(cell_t cell : cells_v){
        for(size_t vert : cell) {
            if( vert >= points_v.size() ) { zero_index = false; break; }
        }
        if ( ! zero_index ) break;
    }

    if ( ! zero_index ){
        for(int i = 0 ; i < cells_v.size() ; ++i){
            for(int j = 0; j < cells_v[i].size(); ++j)
                cells_v[i][j] -= 1;
        }
    }

    t0 = clock();
    std::shared_ptr<gsimp::simplicial_complex> s_comp(new gsimp::simplicial_complex(points_v,cells_v));
    t1 = clock();
    std::cout << "created complex in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    std::shared_ptr<gsimp::path_snapper> p_snap(new gsimp::path_snapper(s_comp));
    t1 = clock();
    std::cout << "created snapper in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    std::vector<size_t> snapped = p_snap->snap_path_to_indices(cycle_points);
    t1 = clock();
    std::cout << "snapped path in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";
    std::cout << "computed path contains " << snapped.size() << " points\n";

    t0 = clock();
    gsimp::chain_t cycle = p_snap->index_sequence_to_chain(snapped);
    t1 = clock();
    std::cout << "produced chain from path in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";

    /*
    t0 = clock();
    gsimp::bounding_chain ch_calc(s_comp);
    t1 = clock();
    std::cout << "calculated boundary matrices in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";

    t0 = clock();
    gsimp::chain_t b_chain_0 = ch_calc.get_bounding_chain(cycle);
    t1 = clock();
    std::cout << "calculated bounding chain (using Eigen) in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";
    */

    gsimp::cell_t null_cell = s_comp->index_to_cell(s_comp->dimension(),null_face);
    std::cout << "null cell: ";
    std::copy(null_cell.begin(),null_cell.end(),std::ostream_iterator<size_t>(std::cout," "));
    std::cout << "\n";

    t0 = clock();
    gsimp::chain_t b_chain_1 = gsimp::coeff_flow(*s_comp,cycle,null_cell,0);
    t1 = clock();
    std::cout << "calculated bounding chain (using coeff_flow) in " << t1 - t0 << " clock cycles " << float(t1-t0)/CLOCKS_PER_SEC << " seconds\n";





}

