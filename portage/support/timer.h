/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include <chrono>
#include <cstdio>
#include <fstream>

/* Wrapper for high precision time point */ 
namespace timer {

/*!
 @brief Get current time
 */ 
inline std::chrono::high_resolution_clock::time_point now() { 
  return std::chrono::high_resolution_clock::now(); 
}

/*!
 @brief Get elapsed time in seconds.
 @param[in] tic start time  
 */ 
inline float elapsed(
  std::chrono::high_resolution_clock::time_point& tic, bool reset = false
) {
  auto const toc = now();
  auto const timing = static_cast<float>(
    std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()
  ) / 1.E3;

  if (reset) { tic = now(); }
  return timing;
}

/*!
 @brief Dump current time and reset it afterwards.
 @param[in] tic start time  
 */ 
inline void reset(
  std::chrono::high_resolution_clock::time_point& tic, float* cumul = nullptr
) {
  if (cumul != nullptr) {
    *cumul += elapsed(tic);
  }
  tic = now();
}

} // namespace 'Portage::timer'

class Profiler {

public:
  // time counters
  struct {
    float mesh_init   = 0;
    float redistrib   = 0;
    float interface   = 0;
    float search      = 0;
    float intersect   = 0;
    float gradient    = 0;
    float interpolate = 0;
    float mismatch    = 0;
    float remap       = 0;
    float total       = 0;
  } time;

  // run parameters
  struct {
    int ranks    = 1;
    int threads  = 1;
    int dim      = 2;
    int nsource  = 0;
    int ntarget  = 0;
    int nmats    = 0;
    int order    = 1;
    std::string output = "";
  } params;

  // constructors
  Profiler() = default;
  Profiler(Profiler const&) = delete;
  Profiler(Profiler&&) noexcept = default;
  ~Profiler() = default;

  // reset all counters
  inline void reset() {
    time.mesh_init    = 0;
    time.redistrib    = 0;
    time.interface    = 0;
    time.search       = 0;
    time.intersect    = 0;
    time.gradient     = 0;
    time.interpolate  = 0;
    time.mismatch     = 0;
    time.remap        = 0;
    time.total        = 0;
    params.ranks      = 1;
    params.threads    = 1;
    params.dim        = 2;
    params.nsource    = 0;
    params.ntarget    = 0;
    params.nmats      = 0;
    params.order      = 1;
    params.output     = "";
  }

  // dump all data to file
  inline bool dump() {

    // save timing for each step
    constexpr int const nsteps = 10;
    constexpr float const time_eps = 1.E-4;

    float const elap[nsteps] = {
      std::max(time_eps, time.mesh_init),
      std::max(time_eps, time.redistrib),
      std::max(time_eps, time.interface),
      std::max(time_eps, time.search),
      std::max(time_eps, time.intersect),
      std::max(time_eps, time.gradient),
      std::max(time_eps, time.interpolate),
      std::max(time_eps, time.mismatch),
      std::max(time_eps, time.remap),
      std::max(time_eps, time.total - time.mesh_init - time.interface - time.remap)
    };

    int time_ratio[nsteps];
    int max_elap = 0;

    for (int i = 0; i < nsteps; ++i) {
      time_ratio[i] = static_cast<int>(100 * (elap[i] / time.total));
      max_elap = std::max(max_elap, static_cast<int>(elap[i]));
    }

    // for number formatting
    int const n_dec = 3;
    int const n_tot = 1 + n_dec + (max_elap > 0 ? ((int) std::floor(std::log10(max_elap))) + 1 : 0);

    std::string const& path = params.output;
    std::printf("\nRecap: total elapsed time %.3f s\n", time.total);
    std::printf(" \u2022 %2d %% generate mesh       \e[32m(%*.3f s)\e[0m.\n", time_ratio[0], n_tot, elap[0]);
    std::printf(" \u2022 %2d %% redistribute        \e[32m(%*.3f s)\e[0m.\n", time_ratio[1], n_tot, elap[1]);
    std::printf(" \u2022 %2d %% interface reconst.  \e[32m(%*.3f s)\e[0m.\n", time_ratio[2], n_tot, elap[2]);
    std::printf(" \u2022 %2d %% search              \e[32m(%*.3f s)\e[0m.\n", time_ratio[3], n_tot, elap[3]);
    std::printf(" \u2022 %2d %% intersect           \e[32m(%*.3f s)\e[0m.\n", time_ratio[4], n_tot, elap[4]);
    std::printf(" \u2022 %2d %% gradient            \e[32m(%*.3f s)\e[0m.\n", time_ratio[5], n_tot, elap[5]);
    std::printf(" \u2022 %2d %% interpolate         \e[32m(%*.3f s)\e[0m.\n", time_ratio[6], n_tot, elap[6]);
    std::printf(" \u2022 %2d %% mismatch            \e[32m(%*.3f s)\e[0m.\n", time_ratio[7], n_tot, elap[7]);
    std::printf(" \u2022 %2d %% remap               \e[32m(%*.3f s)\e[0m.\n", time_ratio[8], n_tot, elap[8]);
    std::printf(" \u2022 %2d %% post-process        \e[32m(%*.3f s)\e[0m.\n", time_ratio[9], n_tot, elap[9]);
    std::fflush(stdout);

    std::printf("Exporting stats to '%s' ... ", path.data());
    std::fflush(stdout);

    auto tic = timer::now();

    std::ofstream file(path, std::ios::out|std::ios::app);
    if (not file.good()) {
      std::fprintf(stderr, "Could not open file :%s\n", path.data());
      return false;
    }

    // add header if file is empty
    std::ifstream checkfile(path);
    assert(checkfile.good());
    checkfile.seekg(0, std::ios::end);
    bool const is_empty = checkfile.tellg() == 0;
    checkfile.close();

    // generate headers if required
    if (is_empty) {
      file << "# Profiling data for t-junction app" << std::endl;
      file << "#"                                   << std::endl;
      file << "# Fields"                            << std::endl;
      file << "#  1. number of ranks"               << std::endl;
      file << "#  2. number of threads"             << std::endl;
      file << "#  3. initialization time"           << std::endl;
      file << "#  4. redistribution time"           << std::endl;
      file << "#  5. interface reconstruction time" << std::endl;
      file << "#  6. search time"                   << std::endl;
      file << "#  7. intersection time"             << std::endl;
      file << "#  8. gradient time"                 << std::endl;
      file << "#  9. interpolation time"            << std::endl;
      file << "# 10. mismatch time"                 << std::endl;
      file << "# 11. remapping time"                << std::endl;
      file << "# 12. total elapsed time"            << std::endl;
      file << "# 13. mesh dimension"                << std::endl;
      file << "# 14. source cells count"            << std::endl;
      file << "# 15. target cells count"            << std::endl;
      file << "# 16. materials count"               << std::endl;
      file << "# 17. remap order"      << std::endl << std::endl;
    }

    file << params.ranks      << "\t"
         << params.threads    << "\t"
         << time.mesh_init    << "\t"
         << time.redistrib    << "\t"
         << time.interface    << "\t"
         << time.search       << "\t"
         << time.intersect    << "\t"
         << time.gradient     << "\t"
         << time.interpolate  << "\t"
         << time.gradient     << "\t"
         << time.remap        << "\t"
         << time.total        << "\t"
         << params.dim        << "\t"
         << params.nsource    << "\t"
         << params.ntarget    << "\t"
         << params.nmats      << "\t"
         << params.order      << std::endl;

    file.close();
    std::printf("done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
    std::fflush(stdout);
    return true;
  }
};
