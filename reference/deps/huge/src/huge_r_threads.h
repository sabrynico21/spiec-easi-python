// R-only OpenMP budget guard.
//
// The standalone C++ core is also consumed by Python and CMake clients, so
// the package-owned OpenMP budget is bounded at the R adapter boundary.
#pragma once

#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace huge {

class RThreadLimitGuard {
public:
    RThreadLimitGuard() noexcept
    {
#ifdef _OPENMP
        previous_max_threads_ = omp_get_max_threads();
        const char* marker = std::getenv("HUGE_R_FORK_WORKER");
        const bool fork_worker =
            marker != nullptr && marker[0] == '1' && marker[1] == '\0';
        const int limit = fork_worker ? 1 : 2;
        if (previous_max_threads_ > limit)
            omp_set_num_threads(limit);
#endif
    }

    ~RThreadLimitGuard() noexcept
    {
#ifdef _OPENMP
        omp_set_num_threads(previous_max_threads_);
#endif
    }

    RThreadLimitGuard(const RThreadLimitGuard&) = delete;
    RThreadLimitGuard& operator=(const RThreadLimitGuard&) = delete;

private:
#ifdef _OPENMP
    int previous_max_threads_ = 1;
#endif
};

} // namespace huge
