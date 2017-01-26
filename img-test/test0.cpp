// generic image libraries
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
// vector arythmetics
#include <Eigen/Dense>
using namespace boost::gil;
#include <iostream>
#include <cmath>
#include <vector>

class Segment {
    Eigen::Vector2d p0, p1;
    Eigen::Vector2d director, normal;

   public:
    Segment(Eigen::Vector2d arg_p0, Eigen::Vector2d arg_p1) {
        p0 = arg_p0;
        p1 = arg_p1;

        assert(length() > 0 && "Points have to be different");

        director = p1 - p0;
        director = (1 / director.norm()) * director;
        normal << director(1), -director(0);
        normal = (1 / normal.norm()) * normal;
    }
    float length() { return (p1 - p0).norm(); }

   private:
    bool is_in_segment(Eigen::Vector2d point, float thickness) {
        bool is_in_segment_v;
        is_in_segment_v = fabs(normal.dot(point - p0)) < thickness / 2  //
                          && director.dot(point - p0) > -thickness / 2  //
                          && -director.dot(point - p1) > -thickness / 2;
        return is_in_segment_v;
    }

   public:
    void draw(rgba8_image_t& img,         //
              const rgba8_pixel_t color,  //
              const float thickness       //
              ) {
        int width = img.width();
        int height = img.height();
        auto v = view(img);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Eigen::Vector2d point;
                point << x, y;
                if (is_in_segment(point, thickness)) {
                    v(x, y) = color;
                }
            }
        }
    }
};

class Triangle {
    Eigen::Vector2d t0, t1, t2;

   private:
    Eigen::Vector2d interior_0, interior_1, interior_2;

   public:
    Triangle(const Eigen::Vector2d arg_t0,  //
             const Eigen::Vector2d arg_t1,  //
             const Eigen::Vector2d arg_t2) {
        t0 = arg_t0;
        t1 = arg_t1;
        t2 = arg_t2;

        Eigen::Vector2d director_0 = (t1 - t0);
        Eigen::Vector2d director_1 = (t2 - t1);
        Eigen::Vector2d director_2 = (t0 - t2);

        interior_0 << director_0(1), -1 * director_0(0);
        if (interior_0.dot(director_1) < 0) {
            interior_0 = -1 * interior_0;
        }
        assert(interior_0.dot(director_1) > 0 &&
               "points are not in general position");
        interior_0 = (1 / interior_0.norm()) * interior_0;

        interior_1 << director_1(1), -1 * director_1(0);
        if (interior_1.dot(director_2) < 0) {
            interior_1 = -1 * interior_1;
        }
        interior_1 = (1 / interior_1.norm()) * interior_1;

        interior_2 << director_2(1), -1 * director_2(0);
        if (interior_2.dot(director_0) < 0) {
            interior_2 = -1 * interior_2;
        }
        interior_2 = (1 / interior_2.norm()) * interior_2;
    }

   private:
    bool is_in_border(const Eigen::Vector2d point, const float thickness) {
        bool is_in_border_v = false;
        // do it all at once so that the compiler can optimize
        is_in_border_v =
            ((fabs(interior_0.dot(point - t0)) < thickness / 2)    //
             && (interior_1.dot(point - t1) > -thickness / 2)      //
             && (interior_2.dot(point - t2) > -thickness / 2)) ||  //
            ((fabs(interior_1.dot(point - t1)) < thickness / 2)    //
             && (interior_2.dot(point - t2) > -thickness / 2)      //
             && (interior_0.dot(point - t0) > -thickness / 2)) ||  //
            ((fabs(interior_2.dot(point - t2)) < thickness / 2)    //
             && (interior_1.dot(point - t1) > -thickness / 2)      //
             && (interior_0.dot(point - t0) > -thickness / 2));
        return is_in_border_v;
    };

    bool is_inside(const Eigen::Vector2d point) {
        bool is_inside_v = (interior_0.dot(point - t0) > 0)     //
                           && (interior_1.dot(point - t1) > 0)  //
                           && (interior_2.dot(point - t2) > 0);
        return is_inside_v;
    }

    void draw_border(rgba8_image_t& img,         //
                     const rgba8_pixel_t color,  //
                     const float thickness       //
                     ) {
        // getting the view
        auto v = view(img);
        // bounds
        int width = img.width();
        int height = img.height();

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Eigen::Vector2d point;
                point << x, y;
                if (is_in_border(point, thickness)) {
                    v(x, y) = color;
                }
            }
        }
    };

    void draw_interior(rgba8_image_t& img,        //
                       const rgba8_pixel_t color  //
                       ) {
        // getting the view
        auto v = view(img);
        // bounds
        int width = img.width();
        int height = img.height();

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Eigen::Vector2d point;
                point << x, y;
                if (is_inside(point)) {
                    v(x, y) = color;
                }
            }
        }
    }

   public:
    void draw(rgba8_image_t& img,                  //
              const rgba8_pixel_t interior_color,  //
              const rgba8_pixel_t border_color,    //
              const float border_thickness)        //
    {
        // first draw the interior of the triangle
        draw_interior(img, interior_color);
        // after the interior add the border
        draw_border(img, border_color, border_thickness);
    };

    float area() {
        Eigen::Vector2d d1{t1 - t0}, d2{t2 - t0};
        return (d1(0) * d2(1) - d1(1) * d2(0)) / 2;
    }
};

int main() {
    Eigen::Vector2d a0(10, 10), a1(100, 10), a2(10, 100);
    Triangle T(a0, a1, a2);
    rgba8_image_t img(512, 512);
    rgba8_pixel_t red(255, 0, 0, 128);
    fill_pixels(view(img), red);
    T.draw(img, {0, 255, 255, 128}, {0, 0, 0, 255}, 2);
    Triangle T2({100, 10}, {10, 100}, {100, 100});
    T2.draw(img, {0, 255, 255, 128}, {0, 0, 0, 255}, 2);

    Segment s({100, 230}, {100, 230});
    std::cout << std::endl
              << T.area() << std::endl;

    png_write_view("redsquare.png", const_view(img));
}
