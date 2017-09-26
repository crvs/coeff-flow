#pragma once
// chains.hpp
// This file presents the utilities to manipulate the chain class. This class
// unifies the approach of using either the Eigen sparse vector class or the
// dense vector, and allows for converting between the two.

#include <vector>
#include <Eigen/Sparse>

using namespace std;

namespace gsimp {

class chain_dimension_mismatch : public std::exception {};
class chain_type_mismatch : public std::exception {};

enum ChainTypes { SPARSE, VECTOR };

class chain {
    int dim;
    ChainTypes type;
    Eigen::SparseVector<double> vs;
    vector<double> vd;

   public:
    chain() {}

    chain(int d, vector<double> vec) {
        dim = d;
        type = VECTOR;
        vd = vec;
    }

    chain(int d, Eigen::SparseVector<double> vec) {
        dim = d;
        type = SPARSE;
        vs = vec;
    }

    chain(const chain& other) {
        dim = other.dim;
        type = other.type;
        if (type == SPARSE) vs = other.vs;
        else  vd = other.vd;

    }

    size_t get_size() {
        if (type == VECTOR)
            return vd.size();
        else
            return vs.innerSize();
    }

    bool is_dense() { return type == VECTOR; }
    bool is_sparse() { return type == SPARSE; }
    vector<double> get_dense() {
        if (type == VECTOR) {
            return vd;
        } else {
            vector<double> vec(vs.innerSize());
            for (auto i = 0; i < vs.innerSize(); i++) {
                if (vs.coeffRef(i) != 0) vec[i] = vs.coeffRef(i);
            }
            return vec;
        }
    }

    Eigen::SparseVector<double> get_sparse() {
        if (type == SPARSE) {
            return vs;
        } else {
            Eigen::SparseVector<double> vec(vd.size());
            for (auto i = 0; i < vd.size(); i++) {
                if (vd[i] != 0) vec.coeffRef(i) = vd[i];
            }
            return vec;
        }
    };

    void to_dense() {
        if (type == SPARSE) {
            vd = get_dense();
            type = VECTOR;
            vs.setZero();
        }
    }

    void to_sparse() {
        if (type == VECTOR) {
            vs = get_sparse();
            type = SPARSE;
            vd.clear();
        }
    }

    double& operator[](const size_t& i) {
        if (type == SPARSE)
            return vs.coeffRef(i);
        else
            return vd.at(i);
    }

    int dimension() { return dim; }

    chain operator+(const chain& other) {
        if (dim != other.dim)
            throw chain_dimension_mismatch();
        else if (type != other.type)
            throw chain_type_mismatch();
        else if (type == SPARSE)
            return chain(dim, vs + other.vs);
        else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = vd[i] + other.vd[i];
            }
            return chain(dim, vec);
        }
    }

    chain operator-(const chain& other) {
        if (dim != other.dim)
            throw chain_dimension_mismatch();
        else if (type != other.type)
            throw chain_type_mismatch();
        else if (type == SPARSE)
            return chain(dim, vs + other.vs);
        else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = vd[i] - other.vd[i];
            }
            return chain(dim, vec);
        }
    }

    chain operator*(const double& d) {
        if (type == SPARSE)
            return chain(dim, d * vs);
        else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = d * vd[i];
            }
            return chain(dim, vec);
        }
    }

    chain abs() {
        if (type == SPARSE) {
            auto vec = vs;
            Eigen::SparseVector<double>::InnerIterator vit(vec);
            while ((bool)vit) {
                vit.valueRef() = std::abs(vit.value());
            }
            return chain(dim, vec);
        } else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = std::abs(vd[i]);
            }
            return chain(dim, vec);
        }
    }

    chain& operator+=(const chain& other) {
        if (dim != other.dim)
            throw chain_dimension_mismatch();
        else if (type != other.type)
            throw chain_type_mismatch();
        else if (type == SPARSE)
            vs += other.vs;
        else {
            for (size_t i = 0; i < vd.size(); i++) 
                vd[i] += other.vd[i];
        }
        return *this;
    }

    chain& operator-=(const chain& other) {
        if (dim != other.dim)
            throw chain_dimension_mismatch();
        else if (type != other.type)
            throw chain_type_mismatch();
        else if (type == SPARSE)
            vs -= other.vs;
        else {
            for (size_t i = 0; i < vd.size(); i++) 
                vd[i] -= other.vd[i];
        }
        return *this;
    }

    double operator^(const chain& other) {
        if (dim != other.dim)
            throw chain_dimension_mismatch();
        else if (type != other.type)
            throw chain_type_mismatch();
        else if (type == VECTOR) {
            double prod = 0;
            for (auto i = 0; i < vd.size(); i++) {
                prod += vd[i] * other.vd[i];
            }
            return prod;
        } else {
            double prod = 0;
            Eigen::SparseVector<double>::InnerIterator it1(vs);
            Eigen::SparseVector<double>::InnerIterator it2(other.vs);
            while ((bool)it1) {
                prod += it1.value() * it2.value();
            }
            return prod;
        }
    }
};
};  // namespace gsimp
