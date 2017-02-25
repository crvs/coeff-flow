#pragma once

#include<limits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>

#include <vector>
#include <string>

namespace gsimp {

class drawing {
    struct impl;
    std::shared_ptr<impl> pimpl;
    std::string filename, temp_filename;
   public:
    drawing(std::string, int, int);
    drawing &operator=(const drawing &other);
    void add_arrow_definition();
    void draw_arrow(std::vector<double>, std::vector<double>, std::string, double);
    void draw_edge(std::vector<double>, std::vector<double>, std::string, double);
    void draw_triangle(std::vector<std::vector<double>>, std::string);
};

struct drawing::impl {
    std::ofstream file;
    double maxx = -std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();
    double minx = std::numeric_limits<double>::infinity();
    double miny = std::numeric_limits<double>::infinity();
    int size_x, size_y;
    std::string element_definitions;
    bool has_arrows;
    std::string filename;
    std::string temp_filename;


    impl(std::string filename, int sizex, int sizey)
        : size_x(sizex), size_y(sizey) {
        this->filename = filename + ".svg";
        temp_filename = filename + "-temp.svg";
        file.open(temp_filename);
    }

    ~impl() {
        file.close();
        std::ifstream temp_file(temp_filename);
        std::ofstream file_finished;
        file_finished.open(filename);

        file_finished << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' viewBox='"    //
                      << 0 << " "            //
                      << 0 << " "            //
                      << maxx - minx << " "  //
                      << maxy - miny << "' height = '500px' width='500px'>\n";
        file_finished << "<g transform='translate(" << -minx << "," << -miny << ")' >\n";
        for (std::string file_line; std::getline(temp_file, file_line);) {
            file_finished << file_line + "\n";
        }
        file_finished << "</g>\n";
        file_finished << "</svg>\n";
        file_finished.close();
        std::remove(temp_filename.c_str());
    }

    void adjust_to_point(std::vector<double> point) {
        minx = fmin(point[0] , minx);
        maxx = fmax(point[0] , maxx);
        miny = fmin(point[1] , miny);
        maxy = fmax(point[1] , maxy);
    }
};

drawing& drawing::operator=(const drawing& other) {
    pimpl = std::move(other.pimpl);
    return *this;
}

drawing::drawing(std::string filename, int sizex, int sizey) {
    pimpl = std::shared_ptr<impl>(new impl(filename, sizex, sizey));
}

void drawing::draw_triangle(std::vector<std::vector<double>> coords,
                            std::string color) {
    if (coords.size() != 3) throw "A triangle should have 3 vertices";

    std::string tri("<polyline points='");
    for (auto coord : coords) {
        double x(coord[0]), y(coord[1]);
        tri += std::to_string(x) + ", " + std::to_string(y) + " ";
        pimpl->adjust_to_point(coord);
    }
    tri += "' stroke='none' fill='" + color + "'/>\n";
    pimpl->file << tri;
}

void drawing::add_arrow_definition() {
    std::string arrow_defn = "<defs>\n";                               //
    arrow_defn += "<marker id='head' orient= auto ";                   //
    arrow_defn += "markerHeight='6' markerWidth='6'",                  //
        arrow_defn += "refX='0' refY='1.5'>\n";                        //
    arrow_defn += "<path d='M 0 0 V 3 L 4 1.5 Z' fill='#ff0000'/>\n";  //
    arrow_defn += "</marker>\n";                                       //
    arrow_defn += "</defs>\n";

    pimpl->file << arrow_defn;
}

void drawing::draw_edge(std::vector<double> src, std::vector<double> trg,
                        std::string color, double stroke) {
    std::string line_s = "<path stroke-linecap='round' ";  //
    line_s += "stroke-width='";                                      //
    line_s += std::to_string(stroke);                               //
    line_s += "' fill='none' stroke='";                              //
    line_s += color;                                                //
    line_s += "' d='M";                                          //

    std::string bgn(""), end(",");
    for (int i = 0; i < src.size(); ++i) {
        bgn += std::to_string(src[i]) + " ";
        end += std::to_string(trg[i]) + " ";
    }
    line_s += bgn + end + "'/>\n";
    pimpl->file << line_s;
}

void drawing::draw_arrow(std::vector<double> src, std::vector<double> trg,
                         std::string color, double stroke) {
    std::string arrow =                                              //
        "<path id = 'arrow' stroke-linecap='round' stroke-width='";  //
    arrow += std::to_string(stroke);                                 //
    arrow += "' fill='none' stroke=";                                //
    arrow += color;                                                  //
    arrow += " d = 'M ";                                             //

    std::string bgn(" "), mid(" L "), end(" L ");
    for (int i = 0; i < src.size(); ++i) {
        bgn += std::to_string(src[i]) + " ";
        mid += std::to_string(trg[i] - src[i] / 2) + " ";
        end += std::to_string(trg[i]) + "' ";
    }
    arrow += bgn + mid + end + " />";
    (pimpl->file) << arrow;
}
};
