#include <boost/function.hpp>

#include <scomplex/quotient.hpp>
#include <scomplex/simplicial_complex.hpp>
#include <scomplex/chains.hpp>

using namespace std;

namespace gsimp {
struct quotient::impl {
    vector<size_t> point_map;
    vector<size_t> point_reverse_map;
    shared_ptr<simplicial_complex> base_comp, quot_comp;
    point_t base_point;

    impl(shared_ptr<simplicial_complex> s_comp, boost::function<bool(point_t)> f) {
        // translate point indices
        size_t i = 1;
        size_t j = 1;
        for (auto pt : s_comp->get_points()) {
            if (f(pt)) {
                point_map.push_back(0);
            } else {
                point_map.push_back(i++);
                point_reverse_map.push_back(j);
            }
            j++;
        }

        // translate faces
        i = 0;
        vector<cell_t> faces;
        while (s_comp->get_level_size(i) > 0) {
            vector<cell_t> level = s_comp->get_level(i);
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

    impl(shared_ptr<simplicial_complex> s_comp, boost::function<bool(point_t)> f,
         point_t b_point) {
        impl(s_comp, f);
        this->base_point = b_point;
    }

    impl& operator=(const impl other) {
        point_map = other.point_map;
        point_reverse_map = other.point_reverse_map;
        base_comp = other.base_comp;
        quot_comp = other.quot_comp;
        base_point = other.base_point;
        return *this;
    }

    /*
     * translate faces back and forth
     */
    cell_t to_q_face(cell_t face) {
        cell_t q_face;
        bool zero = false;
        for (auto v : face) {
            size_t qv = point_map[v];
            if (qv == 0)
                zero = true;
            else
                q_face.push_back(qv);
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

quotient::quotient(shared_ptr<simplicial_complex> s_comp, boost::function<bool(point_t)> quot,
                   point_t pt) {
    unique_ptr<impl> p_impl(new impl(s_comp, quot, pt));
}

quotient::quotient(shared_ptr<simplicial_complex> s_comp, boost::function<bool(point_t)> quot) {
    unique_ptr<impl> p_impl(new impl(s_comp, quot));
}

quotient::quotient(const simplicial_complex& s_comp, boost::function<bool(point_t)> quot,
                   point_t pt) {
    shared_ptr<simplicial_complex> comp_ptr =
        shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    unique_ptr<impl> p_impl(new impl(comp_ptr, quot, pt));
}

quotient::quotient(const simplicial_complex& s_comp, boost::function<bool(point_t)> quot) {
    shared_ptr<simplicial_complex> comp_ptr =
        shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    unique_ptr<impl> p_impl(new impl(comp_ptr, quot));
}

quotient::~quotient() {}

quotient::quotient(const quotient& other) {
    p_impl = unique_ptr<quotient::impl>(new impl(*(other.p_impl)));
}

quotient& quotient::operator=(const quotient& other) {
    p_impl = unique_ptr<quotient::impl>(new impl(*(other.p_impl)));
    return *this;
}

shared_ptr<simplicial_complex> quotient::base_complex() {
    return p_impl->base_comp;
}

shared_ptr<simplicial_complex> quotient::quotient_complex() {
    return p_impl->quot_comp;
}

chain quotient::quotient_chain(chain rep) {
    chain q_rep = p_impl->quot_comp->new_chain(rep.dimension());
    if (rep.is_dense()) q_rep.to_dense();

    for (cell_t cell : p_impl->base_comp->get_level(rep.dimension())) {
        cell_t q_cell = p_impl->to_q_face(cell);
        if (cell.size() == q_cell.size()) {
            size_t base_ind = p_impl->base_comp->cell_to_index(cell);
            size_t quot_ind = p_impl->quot_comp->cell_to_index(q_cell);
            q_rep[quot_ind] = rep[base_ind];
        }
    }
    return q_rep;
}

chain quotient::unquotient_chain(chain q_rep) {
    chain rep = p_impl->base_comp->new_chain(q_rep.dimension());
    if (q_rep.is_dense()) rep.to_dense();

    for (cell_t q_cell : p_impl->quot_comp->get_level(q_rep.dimension())) {
        cell_t cell = p_impl->to_b_face(q_cell);
        size_t base_ind = p_impl->base_comp->cell_to_index(cell);
        size_t quot_ind = p_impl->quot_comp->cell_to_index(q_cell);
        rep[base_ind] = q_rep[quot_ind];
    }
    return q_rep;
}
};
