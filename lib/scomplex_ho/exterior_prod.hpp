#pragma once

#include <vector>
#include <map>
#include <math.h>
#include <string>

/*
 * this file intends to deal with exterior products, mostly for computing
 * volumes of arbitrarily dimensional simplices
 */

using namespace std;

int factorial(int i) {
    int r = 1;
    for (int j = 1; j <= i; j++) r *= j;
    return r;
}

class ext_vector {
    int dim;
    int algebraic_dim;
    map<vector<int>, double> index_map;

    pair<int, vector<int>> index_merge(int ind_1, vector<int> ind_2) {
        vector<int> l_2 = ind_2;
        int j = 0;
        for (int i = 0; i < ind_2.size(); i++) {
            j = i;
            if (l_2[i] >= ind_1)
                break;
            else if (i == ind_2.size() - 1)
                j++;
        }
        if (l_2[j] == ind_1) {
            pair<int, vector<int>> ret_val(0, l_2);
            return ret_val;
        } else {
            l_2.insert(l_2.begin() + j, ind_1);
            j = pow(-1, j);
            pair<int, vector<int>> ret_val(j, l_2);
            return ret_val;
        }
    }

    pair<int, vector<int>> index_merge(vector<int> ind_1, vector<int> ind_2) {
        vector<int> l_1 = ind_1;
        vector<int> l_2 = ind_2;
        int j = 1;
        while (l_1.size() > 0) {
            int k = l_1.back();
            l_1.pop_back();
            auto merge = index_merge(k, l_2);
            j *= get<0>(merge);
            l_2 = get<1>(merge);
        }

        pair<int, vector<int>> ret_val(j, l_2);
        return ret_val;
    };

   public:
    ext_vector(vector<double> vec) {
        for (int i = 0; i < vec.size(); i++) {
            if (vec[i] != 0)
                index_map.insert(
                    pair<vector<int>, double>(vector<int>({i}), vec[i]));
        }
        dim = vec.size();
        algebraic_dim = 1;
    };

    ext_vector operator^(const ext_vector other) {
        ext_vector prod = ext_vector({});
        prod.dim = dim;
        prod.algebraic_dim = algebraic_dim + other.algebraic_dim;
        for (auto p1 : index_map)
            for (auto p2 : other.index_map) {
                auto m = index_merge(p1.first, p2.first);
                double v = get<0>(m) * p1.second * p2.second;
                if (get<1>(m).size() == prod.algebraic_dim)
                    prod.index_map[get<1>(m)] += v;
            }
        return prod;
    };

    double norm() {
        double norm = 0;
        for (auto p : index_map) norm += pow(p.second, 2);
        return sqrt(norm);
    }

    string to_str() {
        string str;
        for (auto el : index_map) {
            str += '[';
            for (auto i : el.first) str += to_string(i) + ' ';
            str += "\b]:" + to_string(el.second) + ' ';
        }
        str += '\n';
        return str;
    }

    int algebraic_dimension() { return algebraic_dim; }
    int dimension() {return dim;}

    /* TODO
     * poincare dual
     */
};
