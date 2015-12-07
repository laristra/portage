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

  template<typename UnaryOp>
    counting_iterator transform(counting_iterator first1, counting_iterator last1,
				counting_iterator result, UnaryOp op) {
    return thrust::transform(first1, last1, result, op);
  }
  template<typename BinaryOp>
    counting_iterator transform(counting_iterator first1, counting_iterator last1,
				counting_iterator first2,
				counting_iterator result, BinaryOp op) {
    return thrust::transform(first1, last1, first2, result, op);
  }

#else // no thrust
  
  template<typename T>
    using vector = std::vector<T>;
  
  typedef boost::counting_iterator<int> counting_iterator;
  counting_iterator make_counting_iterator(int const i) {
    return boost::make_counting_iterator<int>(i);
  }

  template<typename UnaryOp>
    counting_iterator transform(counting_iterator first1, counting_iterator last1,
				counting_iterator result, UnaryOp op) {
    return std::transform(first1, last1, result, op);
  }
  template<typename BinaryOp>
    counting_iterator transform(counting_iterator first1, counting_iterator last1,
				counting_iterator first2,
				counting_iterator result, BinaryOp op) {
    return std::transform(first1, last1, first2, result, op);
  }

#endif
  
}

#endif // PORTAGE_H
