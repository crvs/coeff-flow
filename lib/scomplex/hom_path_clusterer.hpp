#pragma once

#include <scomplex/simplicial_complex.hpp>

using namespace std;
namespace gsimp {

typedef path_t vector<point_t>;

class value {
    double d_val;
    bool inf;

   public:
    // basic constructors (finite and infinite)
    value(double v) {
        d_val = v;
        inf = false;
    }
    value() { inf = true; }

    // copy constructor
    value(const value& other) {
        d_val = other.d_val;
        inf = other.inf;
    }

    bool is_inf() { return is_inf; }
    double val() { return d_val; }

    void set_val(double d) {
        d_val = d;
        inf = false;
    }
    void mk_infinite() { inf = true; }

    bool operator==(const value& other) {
        if (inf) {
            if (other.inf)
                return true;
            else
                return false;
        } else {
            if (other.inf)
                return false;
            else
                return d_val == other.d_val
        }
    }

    bool operator<(const value& other) {
        if (inf)
            return false;
        else {
            if (other.inf)
                return true;
            else
                return d_val < other.d_val;
        }
    }

    bool operator<=(const value& other) {
        return (*this < other) || (*this == other);
    }

    bool operator>(const value& other) {
        return !(*this < other) && !(*this == other);
    }

    bool operator>=(const value& other) {
        return (*this > other) || (*this == other);
    }

    value operator+(const value& other) {
        value val;
        if (!(inf || other.inf))
            val.set_val(d_val + other.d_val);
        return val;
    }

    value operator-(const value& other) {
        value val;
        if (!(inf || other.inf))
            val.set_val(d_val - other.d_val);
        return val;
    }

    value operator*(const value& other) {
        value val;
        if (!(inf || other.inf))
            val.set_val(d_val * other.d_val);
        return val;
    }

    value operator/(const value& other) {
        value val;
        if (!(inf || other.inf || other.d_val == 0))
            val.set_val(d_val / other.d_val);
        return val;
    }

    value& operator+=(const value& other) {
        if (!(inf || other.inf))
            set_val(d_val + other.d_val);
        else
            mk_infinite();
        return *this;
    }

    value& operator-=(const value& other) {
        if (!(inf || other.inf))
            set_val(d_val - other.d_val);
        else
            mk_infinite();
        return *this;
    }

    value& operator*=(const value& other) {
        if (!(inf || other.inf))
            set_val(d_val * other.d_val);
        else
            mk_infinite();
        return *this;
    }

    value& operator/=(const value& other) {
        if (!(inf || other.inf || other.d_val == 0))
            set_val(d_val / other.d_val);
        else
            mk_infinite();
        return *this;
    }
};

class hom_path_clusterer {
   public:
    hom_path_clusterer(shared_ptr<simplicial_complex>);
    hom_path_clusterer(shared_ptr<simplicial_complex>, vector<path_t>);
    hom_path_clusterer(simplicial_complex&)
        hom_path_clusterer(simplicial_complex&, vector<path_t>)

            hom_path_clusterer(hom_path_clusterer&);

    ~hom_path_clusterer();

    add_path(path_t);
    compute_distance_array();
    vector<double> get_distance_array();
    vector<vector<double>> get_distance_matrix();
};
};
