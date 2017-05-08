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
    std::vector<data_t> tokens;
    data_t val;
    while(sstream >> val) tokens.push_back(val);
    return tokens;
}

template <typename Iterator>
void cout_iterator(Iterator first, Iterator last) {
    std::copy(first, last,
              std::ostream_iterator<
                  typename std::iterator_traits<Iterator>::value_type>(
                  std::cout, " "));
    std::cout << '\n';
}

std::pair<std::vector<point_t>, std::vector<cell_t>> parse_qhull_file(
    std::string filename) {
    std::ifstream file_stream(filename);

    if (!file_stream.is_open()) {
        throw std::runtime_error("failed to open \"" + filename + "\"\n");
    }

    std::string line;

    int dimensionality, total_points{0}, total_cells{0};

    std::vector<point_t> points;
    std::vector<cell_t> cells;

    for (int line_number = 0; line_number < total_cells + total_points + 3;
         ++line_number) {
        std::getline(file_stream, line);
        // case of the first line of the file (get dimensions)
        if (line_number == 0) {
            auto tokens = tokenize<size_t>(line);
            dimensionality = tokens.at(0);
        }
        // second line of the file (get number of points)
        if (line_number == 1) {
            auto tokens = tokenize<size_t>(line);
            total_points = tokens.at(0);
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

std::vector<std::vector<point_t>> parse_track_file(std::string filename) {
    std::ifstream file_stream(filename);

    if (!file_stream.is_open()) {
        throw std::runtime_error("failed to open \"" + filename + "\"\n");
    }

    std::vector<std::vector<point_t>> tracks;

    std::vector<point_t> track;

    for (std::string line; std::getline(file_stream, line);) {

        auto tokens = tokenize<double>(line);

        if (tokens.size() == 1) {
            if (track.size() > 0) {tracks.push_back(track);}
            track = std::vector<point_t>();
        }
        else
            track.push_back(tokens);
    }

    file_stream.close();

    return tracks;
}
