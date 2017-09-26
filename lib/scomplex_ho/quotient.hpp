#pragma once

#include <boost/function.hpp>

#include <scomplex_ho/simplicial_complex.hpp>
#include <scomplex_ho/chains.hpp>

using namespace std;

namespace gsimp {
class quotient {
    struct impl;
    std::shared_ptr<impl> p_impl;

   public:
    quotient(const simplicial_complex&, boost::function<bool(point_t)>);
    quotient(const simplicial_complex&, boost::function<bool(point_t)>,
             point_t);
    quotient(std::shared_ptr<simplicial_complex>,
             boost::function<bool(point_t)>);
    quotient(std::shared_ptr<simplicial_complex>,
             boost::function<bool(point_t)>, point_t);
    quotient(const quotient&);
    quotient& operator=(const quotient&);

    ~quotient();

    std::shared_ptr<simplicial_complex> base_complex();
    std::shared_ptr<simplicial_complex> quotient_complex();

    // push chains forward through quotient
    chain quotient_chain(chain);

    // pull chains backward through quotient
    chain unquotient_chain(chain);

    // faces and indices
    cell_t quotient_face(cell_t);
    cell_t unquotient_face(cell_t);
    size_t base_index(cell_t);
    size_t quotient_index(cell_t);

};  // quotient

struct quotient::impl {
    vector<size_t> point_map;
    vector<size_t> point_reverse_map;
    std::shared_ptr<simplicial_complex> base_comp;
    std::shared_ptr<simplicial_complex> quot_comp;
    bool has_base_point;
    point_t base_point;

    impl(std::shared_ptr<simplicial_complex> s_comp,
         boost::function<bool(point_t)> f)
        : base_comp{s_comp} {
        // translate point indices
        std::vector<point_t> points = base_comp->get_points();
        point_reverse_map.push_back(points.size());
        size_t q_i = 1;
        for (size_t b_i = 0; b_i < points.size(); b_i++) {
            if (f(points[b_i])) {
                point_map.push_back(0);
            } else {
                point_map.push_back(q_i);
                point_reverse_map.push_back(b_i);
                q_i++;
            }
        }

        // translate faces
        size_t i = 0;
        vector<cell_t> faces;
        while (this->base_comp->get_level_size(i) > 0) {
            vector<cell_t> level = this->base_comp->get_level(i);
            faces.reserve(level.size());
            faces.insert(faces.end(), level.begin(), level.end());
            i++;
        }

        vector<cell_t> q_faces;
        for (auto face : faces) {
            q_faces.push_back(to_q_face(face));
        }

        this->quot_comp = std::shared_ptr<simplicial_complex>(
            new simplicial_complex(q_faces));
    }

    impl(std::shared_ptr<simplicial_complex> s_comp,
         boost::function<bool(point_t)> f, point_t b_point) {
        *this = impl(s_comp, f);
        this->base_point = b_point;
        this->has_base_point = false;
    }

    impl& operator=(const impl other) {
        this->point_map = other.point_map;
        this->point_reverse_map = other.point_reverse_map;
        this->base_comp = other.base_comp;
        this->quot_comp = other.quot_comp;
        this->base_point = other.base_point;
        this->has_base_point = other.has_base_point;
        return *this;
    }

    // translate faces back and forth

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
        if (zero)
            q_face.insert(q_face.begin(), 0);
        return q_face;
    }

    int calculate_factor(cell_t cell_i) {
        cell_t cell = cell_i;
        std::sort(cell.begin(),cell.end());
        for (int i = 0 ; i < cell.size() ; i++ ) {
            if ( point_map[cell[i]] == 0 )
                return (i%2 == 0?1:-1);
        }
        return 1;
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

quotient::quotient(std::shared_ptr<simplicial_complex> s_comp,
                   boost::function<bool(point_t)> quot, point_t pt)
    : p_impl{std::shared_ptr<impl>(new impl(s_comp, quot, pt))} {}

quotient::quotient(std::shared_ptr<simplicial_complex> s_comp,
                   boost::function<bool(point_t)> quot)
    : p_impl{std::shared_ptr<impl>(new impl(s_comp, quot))} {}

quotient::quotient(const simplicial_complex& s_comp,
                   boost::function<bool(point_t)> quot, point_t pt) {
    std::shared_ptr<simplicial_complex> comp_ptr =
        std::shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    p_impl = std::shared_ptr<impl>(new impl(comp_ptr, quot, pt));
}

quotient::quotient(const simplicial_complex& s_comp,
                   boost::function<bool(point_t)> quot) {
    std::shared_ptr<simplicial_complex> comp_ptr =
        std::shared_ptr<simplicial_complex>(new simplicial_complex(s_comp));
    p_impl = std::shared_ptr<impl>(new impl(comp_ptr, quot));
}

quotient::~quotient() {}

quotient::quotient(const quotient& other)
    : p_impl{std::shared_ptr<quotient::impl>(new impl(*(other.p_impl)))} {}

quotient& quotient::operator=(const quotient& other) {
    p_impl = std::shared_ptr<quotient::impl>(new impl(*(other.p_impl)));
    return *this;
}

std::shared_ptr<simplicial_complex> quotient::base_complex() {
    return this->p_impl->base_comp;
}

std::shared_ptr<simplicial_complex> quotient::quotient_complex() {
    return this->p_impl->quot_comp;
}

chain quotient::quotient_chain(chain rep) {
    int dim = rep.dimension();
    chain q_rep = this->p_impl->quot_comp->new_dense_chain(dim);
    std::vector<cell_t> cells = this->p_impl->base_comp->get_level(dim);
    for (cell_t cell : cells) {
        cell_t q_cell = this->p_impl->to_q_face(cell);
        if (cell.size() == q_cell.size()) {
            int factor = p_impl->calculate_factor(cell);
            size_t base_ind = this->p_impl->base_comp->cell_to_index(cell);
            size_t quot_ind = this->p_impl->quot_comp->cell_to_index(q_cell);
            q_rep[quot_ind] += factor * rep[base_ind];

            /*
            if( rep[base_ind] != 0 || q_rep[quot_ind] != 0) {
                std::cout << "b " << base_ind << ": ";
                for (auto i : cell) std::cout << i << "[" << p_impl->point_map[i] << "] ";
                std::cout << ": " << rep[base_ind] << " | ";
                std::cout << "q " << quot_ind << ": ";
                for (auto i : q_cell) std::cout << i << " ";
                std::cout << ": " << q_rep[quot_ind] << " | ";
                std::cout << factor << '\n';
            }
            */
        }
    }

    /*
    chain bdry = this->p_impl->quot_comp->new_dense_chain(0);
    for (size_t i = 0 ; i < q_rep.get_size() ; i++) {
        cell_t q_cell = this->p_impl->quot_comp->index_to_cell(1,i);
        std::sort(q_cell.begin(),q_cell.end());
        bdry[q_cell[0]] += -1 * q_rep[i];
        bdry[q_cell[1]] += q_rep[i];
    }
    std::cout << "boundary: \n";
    for (size_t i = 0 ; i < bdry.get_size() ; i++ ) {
        if (bdry[i] != 0)
            std::cout << i << ": " << bdry[i] << "\n";
    }
    */

    return q_rep;
}

chain quotient::unquotient_chain(chain q_rep) {
    chain rep = this->p_impl->base_comp->new_dense_chain(q_rep.dimension());

    for (cell_t q_cell :
         this->p_impl->quot_comp->get_level(q_rep.dimension())) {
        size_t quot_ind = this->p_impl->quot_comp->cell_to_index(q_cell);
        if (q_rep[quot_ind] == 0) continue;
        cell_t cell = this->p_impl->to_b_face(q_cell);
        if (cell.size() == 0) continue;
        size_t base_ind = this->p_impl->base_comp->cell_to_index(cell);
        rep[base_ind] = q_rep[quot_ind];
    }
    return rep;
}

cell_t quotient::quotient_face(cell_t cell) {
    return p_impl->to_q_face(cell);
}

cell_t quotient::unquotient_face(cell_t cell) {
    return p_impl->to_b_face(cell);
}

size_t quotient::base_index(cell_t cell) {
    return p_impl->base_comp->cell_to_index(cell);
}

size_t quotient::quotient_index(cell_t cell) {
    return p_impl->quot_comp->cell_to_index(cell);
}

};
