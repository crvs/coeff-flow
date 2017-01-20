// here we code the trace_ objects they keep a path
#ifndef TRACE_H
#define TRACE_H

#include <vector>
#include <initializer_list>
#include <string>

namespace trace {
typedef std::vector<double> point_t;
class trace {
   private:
    size_t dimension_;
    size_t length_;
    std::vector<point_t> trace_;

   public:
    trace(std::vector<point_t>& list) {
        trace_ = std::vector<point_t>(list);

        // error if we make an empty trace
        assert(trace_.size() > 0);

        dimension_ = trace_.at(0).size();
        length_ = trace_.size();
    }

    trace(trace& other) {
        trace_ = std::vector<point_t>(other.trace_);
        dimension_ = trace_.at(0).size();
        length_ = trace_.size();
    }

    size_t dimension() { return dimension_; };
    size_t size() { return length_; }
    std::vector<point_t> get_trace() { return trace_; }

    std::string print() {
        std::string out_str;
        out_str = "dim: ";
        out_str += std::to_string(dimension_);
        out_str += " size: ";
        out_str += std::to_string(length_);
        out_str += " points";
        for (point_t pt : trace_) {
            out_str += ":";
            for (double coord : pt) {
                out_str += std::to_string(coord) + " ";
            }
        }
        out_str += std::to_string(length_);
    }
};
};
#endif
