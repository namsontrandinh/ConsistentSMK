// src/algs/result.h
#ifndef ALGS_RESULT_H
#define ALGS_RESULT_H

#include <cstdint>
#include <string>

#include "sfunctions.h"

namespace algs {

struct Result {
    std::string algo;
    std::string constraint;

    subm::Solution inS;

    double f_value = 0.0;
    std::uint64_t queries = 0;
    double time_sec = 0.0;
    double mem_mb = 0.0;
};

} // namespace algs

#endif // ALGS_RESULT_H
