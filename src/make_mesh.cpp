#include <scomplex/qhull_parsing.hpp>
#include <fstream>
#include <scomplex/plywriter.hpp>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    string filename = argv[1];

    pair<vector<vector<double>>, vector<vector<size_t>>> mesh;
    mesh = parse_qhull_file(filename);

    vector<vector<double>> points = get<0>(mesh);
    vector<vector<size_t>> edges = {};
    vector<vector<size_t>> faces = get<1>(mesh);

    ofstream ply_file(filename.append(".ply"));
    make_ply(ply_file, points,faces);
    ply_file.close();
    return 0;}





