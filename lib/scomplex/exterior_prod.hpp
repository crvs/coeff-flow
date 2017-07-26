#pragma once

#include <algorithm>
#include <vector>
#include <map>
#include <math>

/*
 * this file intends to deal with exterior products, mostly for computing
 * volumes of arbitrarily dimensional simplices
 */

using namespace std;

class ext_vector {
    int dimension;
    int algebraic_dimension;
    map<vector<int>, double> index_map;

    pair<int, vector<int>> index_merge(int ind_1, vector<int> ind_2) {
        int j = -1;
        for (int i = 0; i < ind_2.size(); i++) {
            if (ind_1 == ind_2[i]) {
                pair<int, vector<int>> ret_val(0, ind_2);
                return ret_val;
            } else if (ind_1 < ind_2[i]) {
                j *= -1;
            } else {
                insert(ind_2.begin() + i, ind_1);
                pair<int, vector<int>> ret_val(j, ind_2);
                return ret_val;
            }
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

        pair<int, vector<int>> ret_val(j, ind_2);
        return ret_val;
    };

   public:
    ext_vector(vector<double> vec) {
        for (int i = 0; i < vec.size(); i++) {
            index_map.insert(
                pair<vector<int>, double>(vector<int>({i}), vec[i]));
        }
        dimension = vec.size();
        algebraic_dimension = 1;
    };

    ext_vector operator^(const ext_vector other) {
        ext_vector prod = ext_vector({});
        prod.dimension = dimension;
        prod.algebraic_dimension =
            algebraic_dimension + other.algebraic_dimension;
        for (auto p1 : index_map)
            for (auto p2 : other.index_map)
            {
                auto m = index_merge(p1.first,p2.first);
                double v = get<0>(m) * p1.second * p2.second;
                prod.index_map[get<1>(m)] += v;
            }
        return prod;
    };

    double norm() {
        double norm = 0;
        for ( auto p : index_map)
            norm += p.second^2;
        return sqrt(norm);
    }
};
