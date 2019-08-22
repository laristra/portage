/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include <chrono>

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
inline int reset(
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
    time.mesh_init = 0;
    time.redistrib = 0;
    time.interface = 0;
    time.remap     = 0;
    time.total     = 0;
    params.ranks   = 1;
    params.threads = 1;
    params.dim     = 2;
    params.nsource = 0;
    params.ntarget = 0;
    params.nmats   = 0;
    params.order   = 1;
    params.output  = "";
  }

  // dump all data to file
  inline bool dump() {

    // save timing for each step
    constexpr int const nsteps = 3;
    constexpr float const time_eps = 1.E-4;

    float const elap[nsteps] = {
      std::max(time_eps, time.mesh_init),
      std::max(time_eps, time.interface),
      std::max(time_eps, time.remap)
    };

    int time_ratio[nsteps];

    for (int i = 0; i < nsteps; ++i) {
      time_ratio[i] = static_cast<int>(elap[i] * 100 / time.total);
    }

    std::string const& path = params.output;

    std::printf("\nRecap: total elapsed time %.3f s\n", time.total);
    std::printf("= %2d %% mesh initialization     (%6.1f s).\n", time_ratio[0], elap[0]);
    std::printf("= %2d %% interface recontruction (%6.1f s).\n", time_ratio[1], elap[1]);
    std::printf("= %2d %% remapping               (%6.1f s).\n", time_ratio[2], elap[2]);
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
      file << "#  6. remapping time"                << std::endl;
      file << "#  7. total elapsed time"            << std::endl;
      file << "#  8. mesh dimension"                << std::endl;
      file << "#  9. source cells count"            << std::endl;
      file << "# 10. target cells count"            << std::endl;
      file << "# 11. materials count"               << std::endl;
      file << "# 12. remap order"      << std::endl << std::endl;
    }

    file << params.ranks   << "\t"
         << params.threads << "\t"
         << time.mesh_init << "\t"
         << time.redistrib << "\t"
         << time.interface << "\t"
         << time.remap     << "\t"
         << time.total     << "\t"
         << params.dim     << "\t"
         << params.nsource << "\t"
         << params.ntarget << "\t"
         << params.nmats   << "\t"
         << params.order   << std::endl;

    file.close();
    std::printf("done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
    std::fflush(stdout);
    return true;
  }
};