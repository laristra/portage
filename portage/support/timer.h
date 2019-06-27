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
inline float elapsed(std::chrono::high_resolution_clock::time_point& tic, bool reset = false) {
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
