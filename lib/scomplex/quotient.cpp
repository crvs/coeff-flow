#pragma once

#include <scomplex/quotient.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/types.hpp>

using namespace std;

struct quotient::impl {
    vector<size_t> point_map;
    vector<size_t> point_reverse_map;
    shared_ptr<simplicial_complex> base_comp, quot_comp;
    point_t base_point;

    impl(shared_ptr<simplicial_complex> s_comp, bool(f)(point_t)) {
        // translate point indices

        size_t i = 1;
        size_t j = 1;
        for (auto pt : s_comp->points) {
            if (f(pt)) {
                point_map.push_back(0);
            } else {
                point_map.push_back(i++);
                point_reverse_map.push_back(j)
            }
            j++;
        }

        // translate faces
        i = 0;
        vector<cell_t> faces;
        while (s_comp->get_level_size(i) > 0) {
            vector<cell_t> level s_comp->get_level(i);
            faces.reserve(level.size());
            faces.insert(faces.end(), level.begin(), level.end());
            i++;
        }

        vector<cell_t> q_faces;
        for (auto face : faces) {
            q_faces.push_back(to_q_face(face));
        }

        this->base_comp = s_comp;
        this->quot_comp =
            shared_ptr<simplicial_complex>(new simplicial_complex(q_faces));
    }

    impl(shared_ptr<simplicial_complex> s_comp, bool(f)(point_t),
         point_t b_point) {
        impl(s_comp, f);
        this->base_point = b_point;
    }

    /*
     * translate faces back and forth
     */
    cell_t to_q_face(cell_t face) {
        cell_t q_face;
        bool zero = false;
        for (auto v : face) {
            size_t qv = point_map[v] if (qv == 0) zero = true;
            else q_face.push_back(qv);
        }
        if (zero) q_face.insert(q_face.begin(), 0);
        return q_face;
    }

    cell_t to_b_face(cell_t q_face) {
        cell_t b_face;
        for (auto v : q_face) {
            if (v == 0)
                return cell_t();
            else
                b_face.push_back(point_reverse_map[v]);
        }
        return b_face;
    }
};

quotient::quotient(shared_ptr<simplicial_complex> s_comp, bool(quot)(point_t),
                   point_t pt) {
    unique_ptr<impl> p_impl(new impl(s_comp, quot, pt));
}

quotient::quotient(shared_ptr<simplicial_complex> s_comp, bool(quot)(point_t)) {
    unique_ptr<impl> p_impl(new impl(s_comp, quot));
}

quotient::quotient(const simplicial_complex& s_comp, bool(quot)(point_t),
                   point_t pt) {
    shared_ptr<simplicial_complex> comp_ptr =
        shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    unique_ptr<impl> p_impl(new impl(comp_ptr, quot, pt));
}

quotient::quotient(const simplicial_complex& s_comp, bool(quot)(point_t)) {
    shared_ptr<simplicial_complex> comp_ptr =
        shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    unique_ptr<impl> p_impl(new impl(comp_ptr, quot));
}

quotient::~quotient() {}

quotient::quotient(const quotient& other) { p_impl = other.p_impl; }

quotient::quotient& operator=(const quotient& other) {
    p_impl = other.p_iml;
    return *this;
}

quotient::std::shared_ptr<simplicial_complex> base_complexi() {
    return p_impl->base_comp;
}

quotient::std::shared_ptr<simplicial_complex> quotient_complex() {
    return p_impl->quot_comp;
}

quotient::chain_t quotient_chain(chain_t chain)
{

}

/* TODO
 * fill in details
quotient::chain_v quotient_chain_v(chain_v);
quotient::chain_t unquotient_chain(chain_t);
quotient::chain_v unquotient_chain_v(chain_v);
*/
