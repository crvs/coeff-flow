#pragma once

#include <fstream>

#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <vector>
#include <tuple>

typedef std::vector<double> point_t;
typedef std::vector<size_t> cell_t;

template <typename data_t>
std::vector<data_t> tokenize(std::string str) {
    std::istringstream sstream(str);
    std::vector<data_t> tokens{std::istream_iterator<data_t>{sstream},
                               std::istream_iterator<data_t>{}};
    return tokens;
}

template <typename data_t>
void output_vector(std::vector<data_t> vec) {
    std::copy(vec.begin(), vec.end(),
              std::ostream_iterator<data_t>(std::cout, " "));
    std::cout << '\n';
}

std::pair<std::vector<point_t>, std::vector<cell_t>> parse_qhull_file(
    std::string filename) {
    std::ifstream file_stream(filename);

    if (!file_stream.is_open()) {
        throw std::runtime_error("failed to open \"" + filename + "\"\n");
    }

    std::string line;

    int current_line_number{0};
    int dimensionality;
    //
    int num_points;
    int num_cells;
    //
    int cur_point_num;
    int cur_cell_num;

    int line_number{0};
    int total_points{0};
    int total_cells{0};

    std::vector<point_t> points;
    std::vector<cell_t> cells;

    for (int line_number = 0; line_number < total_cells + total_points + 3;
         ++line_number) {
        std::getline(file_stream, line);
        // case of the first line of the file (get dimensions)
        if (line_number == 0) {
            auto tokens = tokenize<size_t>(line);
            dimensionality = tokens.at(0);
            // std::cout << "found first line! dimensions == " << dimensionality
            //          << '\n';
        }
        // second line of the file (get number of points)
        if (line_number == 1) {
            auto tokens = tokenize<size_t>(line);
            total_points = tokens.at(0);
            // std::cout << "found second line! number of points == "
            //          << total_points << '\n';
        }
        // lines pertaining to points
        if (line_number > 1 && line_number < total_points + 2) {
            auto point = tokenize<double>(line);
            point.pop_back();
            points.push_back(point);
        }
        // first line after points (get number of cells)
        if (line_number == total_points + 2) {
            auto tokens = tokenize<size_t>(line);
            total_cells = tokens.at(0);
            // std::cout << "found cell count! number of cells == " <<
            // total_cells
            //          << '\n';
        }
        // lines pertaining to cells
        if (line_number > total_points + 2 &&
            line_number < total_points + total_cells + 3) {
            auto cell = tokenize<size_t>(line);
            cells.push_back(cell);
        }
    }
    file_stream.close();
    return std::make_pair(points, cells);
}

/*
 * usage example
 *
int main() {
    std::vector<point_t> points;
    std::vector<cell_t> cells;
    std::tie(points, cells) = parse_qhull_file("qh-test.dat");

    std::cout << points.size() << '\n';
    for (point_t point : points) {
        output_vector<double>(point);
    }
    std::cout << cells.size() << '\n';
    for (cell_t cell : cells) {
        output_vector<size_t>(cell);
    }

    return 0;
}
*/
