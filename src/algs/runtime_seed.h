// src/algs/runtime_seed.h
#ifndef ALGS_RUNTIME_SEED_H
#define ALGS_RUNTIME_SEED_H

#include <cstdint>
#include <chrono>
#include <thread>
#include <functional>
#include <random>

namespace algs {

inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}

inline std::uint64_t runtime_seed64() {
    std::uint64_t s = 0;

    s ^= static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );

    s ^= static_cast<std::uint64_t>(
        std::hash<std::thread::id>{}(std::this_thread::get_id())
    );

    try {
        std::random_device rd;
        const std::uint64_t a = (static_cast<std::uint64_t>(rd()) << 32) ^
                                static_cast<std::uint64_t>(rd());
        s ^= a;
    } catch (...) {
        // ignore
    }

    return splitmix64(s);
}

} // namespace algs

#endif // ALGS_RUNTIME_SEED_H
