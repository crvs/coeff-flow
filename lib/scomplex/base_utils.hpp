#pragma once

#include <set>
#include <vector>
#include <list>

std::list<size_t> dedupe_list(std::list<size_t>& list) {
    // O(n log(n))
    std::set<size_t> no_reps;
    std::list<size_t> no_reps_list;
    for (size_t el : list) {
        no_reps.insert(el);
    }
    for (size_t el : no_reps) {
        no_reps_list.push_back(el);
    }
    return no_reps_list;
}

std::vector<size_t> dedupe_vec(std::vector<size_t>& vec) {
    // O(n log(n))
    std::set<size_t> no_reps;
    std::vector<size_t> no_reps_list;
    for (size_t el : vec) {
        no_reps.insert(el);
    }
    for (size_t el : no_reps) {
        no_reps_list.push_back(el);
    }
    return no_reps_list;
}

template <typename T>
std::unique_ptr<T> new_unique(T* t) {
    T* tt = new T(t);
    return std::unique_ptr<T>(tt);
}
