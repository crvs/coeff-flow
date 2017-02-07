#pragma once

#include <vector>
#include <list>
#include <initializer_list>

std::list<int> dedupe_list(std::list<int>& list) {
    // O(n log(n))
    std::set<int> no_reps;
    std::list<int> no_reps_list;
    for (int el : list) {
        no_reps.insert(el);
    }
    for (int el : no_reps) {
        no_reps_list.push_back(el);
    }
    return no_reps_list;
}

std::list<int> dedupe_list(std::vector<int>& list) {
    // O(n log(n))
    std::set<int> no_reps;
    std::list<int> no_reps_list;
    for (int el : list) {
        no_reps.insert(el);
    }
    for (int el : no_reps) {
        no_reps_list.push_back(el);
    }
    return no_reps_list;
}

std::list<int> dedupe_list(std::initializer_list<int>& list) {
    // O(n log(n))
    std::set<int> no_reps;
    std::list<int> no_reps_list;
    for (int el : list) {
        no_reps.insert(el);
    }
    for (int el : no_reps) {
        no_reps_list.push_back(el);
    }
    return no_reps_list;
}
