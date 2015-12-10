/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef PORTAGE_H
#define PORTAGE_H

#ifdef THRUST

#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#else // no thrust

#include <vector>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

#endif

namespace Portage {

#ifdef THRUST
  
  template<typename T>
    using vector = thrust::device_vector<T>;
  
  typedef thrust::counting_iterator<int> counting_iterator;
  counting_iterator make_counting_iterator(int const i) {
    return thrust::make_counting_iterator(i);
  }

  template<typename InputIterator, typename OutputIterator, typename UnaryFunction>
    OutputIterator transform(InputIterator first, InputIterator last,
                             OutputIterator result, UnaryFunction op) {
    return thrust::transform(first, last, result, op);
  }
  template<typename InputIterator1, typename InputIterator2, typename OutputIterator, typename BinaryFunction>
    OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                             InputIterator2 first2, OutputIterator result, BinaryFunction op) {
    return thrust::transform(first1, last1, first2, result, op);
  }


#else // no thrust
  template<typename T>
    using vector = std::vector<T>;
  
  typedef boost::counting_iterator<int> counting_iterator;
  counting_iterator make_counting_iterator(int const i) {
    return boost::make_counting_iterator<int>(i);
  }

  template<typename InputIterator, typename OutputIterator, typename UnaryFunction>
    OutputIterator transform(InputIterator first, InputIterator last,
		             OutputIterator result, UnaryFunction op) {
    return std::transform(first, last, result, op);
  }
  template<typename InputIterator1, typename InputIterator2, typename OutputIterator, typename BinaryFunction>
    OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
			     InputIterator2 first2, OutputIterator result, BinaryFunction op) {
    return std::transform(first1, last1, first2, result, op);
  }

#endif
  
}

#endif // PORTAGE_H
