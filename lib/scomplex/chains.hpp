#pragma once
// chains.hpp
// This file presents the utilities to manipulate the chain class. This class
// unifies the approach of using either the Eigen sparse vector class or the
// dense vector, and allows for converting between the two.

// TODO:
// * test
// * integrate into existing code
//      - start by ack-grep chain_[tv]
//      - load into quickfix
//      - check each hit

using namespace std;

namespace gsimp {

enum ChainTypes { SPARSE, VECTOR };

class chain {
    int dim;
    ChainTypes type;
    Eigen::SparseVector<double> vs;
    vector<double> vd;

   public:
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
            vs.SetZero();
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
            return &vd.at(i);
    }

    int dimension() { return dim; }

    chain operator+(const chain& other) {
        if (dim != other.dim)
            throw std::exception("chain dimension mismatch in adition");
        else if (type != other.type)
            throw std::exception("chain type mismatch in adition");
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

    chain operator*(const double& d) {
        else if (type == SPARSE) return chain(dim, d * vs);
        else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = d * vd[i];
            }
            return chain(dim, vec);
        }
    }

    chain abs() {
        else if (type == SPARSE) auto vec = vs;
        auto vit = vec.InnerIterator;
        while ((bool)vit) {
            vit.valueRef() = abs(vit.value());
        }
        return chain(dim, vec);
        else {
            vector<double> vec(vd.size());
            for (size_t i = 0; i < vd.size(); i++) {
                vec[i] = abs(vd[i]);
            }
            return chain(dim, vec);
        }
    }

    double operator<*>(const chain& other) {
        if (dim != other.dim)
            throw std::exception("chain dimension mismatch in inner product");
        else if (type != other.type)
            throw std::exception("chain type mismatch in inner product");
        else if (type == VECTOR) {
            double prod = 0;
            for (auto i = 0; i < vd.size(); i++) {
                prod += vd[i] * other.vd[i];
            }
            return prod;
        } else {
            double prod = 0;
            auto it1 = vs.InnerIterator;
            auto it2 = other.vs.InnerIterator;
            while ((bool)it1) {
                prod += it1.value() * it2.value();
            }
            return prod;
        }
    }
};
};  // namespace gsimp
