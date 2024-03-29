/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SUPPORT_BASIS_H_
#define PORTAGE_SUPPORT_BASIS_H_

#include <cassert>
#include <cmath>
#include <array>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

// portage includes
#include "portage/support/portage.h"

namespace Portage { namespace Meshfree { namespace basis {

enum Type { Unitary, Linear, Quadratic, LastBasis };

constexpr std::array<size_t, 4> const quadratic_sizes = {0, 3, 6, 10};

////////////////////////////////////////////////////////////////////////////////
// Traits
////////////////////////////////////////////////////////////////////////////////

template<Type type, int dim>
class Traits {};

template<int dim>
class Traits<Unitary, dim> {
public:
  static constexpr size_t function_size = 1;
  static constexpr std::array<size_t, 2> jet_size = {1, 1};
  using array_t  = std::array<double, function_size>;
  using matrix_t = std::array<std::array<double, function_size>, function_size>;
};

// this definition completes declaration of the static member
template<int dim>
constexpr std::array<size_t, 2> Traits<Unitary, dim>::jet_size;

template<int dim>
class Traits<Linear, dim> {
public:
  static constexpr size_t function_size = dim + 1;
  static constexpr std::array<size_t, 2> jet_size = {dim + 1, dim + 1};
  using array_t  = std::array<double, function_size>;
  using matrix_t = std::array<std::array<double, function_size>, function_size>;
};

// this definition completes declaration of the static member
template<int dim>
constexpr std::array<size_t, 2> Traits<Linear, dim>::jet_size;

template<int dim>
class Traits<Quadratic, dim> {
public:
  static constexpr size_t function_size = quadratic_sizes[dim];
  static constexpr std::array<size_t, 2> jet_size = {function_size, function_size};
  using array_t  = std::array<double, function_size>;
  using matrix_t = std::array<std::array<double, function_size>, function_size>;
};

// this definition completes declaration of the static member, see e.g.
// stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
// this should be revisited after switching to C++17
template<int dim>
constexpr std::array<size_t, 2> Traits<Quadratic, dim>::jet_size;

template<int dim>
size_t function_size(Type btype) {
  switch (btype) {
    case Unitary:   return Traits<Unitary, dim>::function_size;
    case Linear:    return Traits<Linear, dim>::function_size;
    case Quadratic: return Traits<Quadratic, dim>::function_size;
    default: return 0;
  }
}

template<int dim>
std::array<size_t, 2> jet_size(Type btype) {
  switch (btype) {
    case Unitary:   return Traits<Unitary, dim>::jet_size;
    case Linear:    return Traits<Linear, dim>::jet_size;
    case Quadratic: return Traits<Quadratic, dim>::jet_size;
    default: return {0, 0};
  }
}

////////////////////////////////////////////////////////////////////////////////
// Templates
////////////////////////////////////////////////////////////////////////////////

template<Type type, size_t dim>
std::array<double, Traits<type, dim>::function_size>
function(Wonton::Point<dim> x) { return {}; }

template<Type type, size_t dim>
std::array<
  std::array<double, Traits<type, dim>::function_size>,
  Traits<type, dim>::function_size
>
jet(Wonton::Point<dim> x) { return {}; }

template<Type type, size_t dim>
std::array<
  std::array<double, Traits<type, dim>::function_size>,
  Traits<type, dim>::function_size
>
inverse_jet(Wonton::Point<dim> x) { return {}; }

template<Type type, size_t dim>
std::array<double, Traits<type, dim>::function_size>
shift(Wonton::Point<dim> x, Wonton::Point<dim> y) {

  auto ijet = inverse_jet<type, dim>(x);
  auto by   = function<type, dim>(y);
  int const n = Traits<type, dim>::function_size;
  std::array<double, Traits<type, dim>::function_size> r {};

  for (int i = 0; i < n; i++) {
    r[i] = 0.;
    for (int j = 0; j < n; j++) {
      r[i] += ijet[i][j] * by[j];
    }
  }
  return r;
}

template<Type type, size_t dim>
typename Traits<type, dim>::matrix_t transfactor(const Wonton::Point<dim>& c) {
  typename Traits<type, dim>::matrix_t result;
  assert(false);
  return result;
}

template<size_t dim>
std::vector<double> function(Type type, Wonton::Point<dim> x) {
  size_t nbasis = function_size<dim>(type);
  std::vector<double> result(nbasis);
  switch (type) {
    case Unitary: {
      auto resulta = function<Unitary, dim>(x);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    case Linear: {
      auto resulta = function<Linear, dim>(x);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    case Quadratic: {
      auto resulta = function<Quadratic, dim>(x);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
std::vector<double> shift(Type type, Wonton::Point<dim> x, Wonton::Point<dim> y) {
  size_t nbasis = function_size<dim>(type);
  std::vector<double> result(nbasis);
  switch (type) {
    case Unitary: {
      auto resulta = shift<Unitary, dim>(x, y);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    case Linear: {
      auto resulta = shift<Linear, dim>(x, y);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    case Quadratic: {
      auto resulta = shift<Quadratic, dim>(x, y);
      for (size_t i = 0; i < nbasis; i++)
        result[i] = resulta[i];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
std::vector<std::vector<double>> jet(Type type, Wonton::Point<dim> x) {
  auto njet = jet_size<dim>(type);
  std::vector<std::vector<double>> result(njet[0], std::vector<double>(njet[1], 0.));
  switch (type) {
    case Unitary: {
      auto resulta = jet<Unitary, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    case Linear: {
      auto resulta = jet<Linear, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    case Quadratic: {
      auto resulta = jet<Quadratic, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
std::vector<std::vector<double>> inverse_jet(Type type, Wonton::Point<dim> x) {
  auto njet = jet_size<dim>(type);
  std::vector<std::vector<double>> result(njet[0], std::vector<double>(njet[1], 0.));
  switch (type) {
    case Unitary: {
      auto resulta = inverse_jet<Unitary, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    case Linear: {
      auto resulta = inverse_jet<Linear, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    case Quadratic: {
      auto resulta = inverse_jet<Quadratic, dim>(x);
      for (size_t i = 0; i < njet[0]; i++)
        for (size_t j = 0; j < njet[1]; j++)
          result[i][j] = resulta[i][j];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
std::vector<std::vector<double>> transfactor(const Type type, const Wonton::Point<dim>& c) {
  size_t nbasis = function_size<dim>(type);
  std::vector<std::vector<double>> result(nbasis, std::vector<double>(nbasis, 0.));
  switch (type) {
    case Unitary: {
      auto tf = transfactor<Unitary, dim>(c);
      for (size_t i = 0; i < nbasis; i++)
        for (size_t j = 0; j < nbasis; j++)
          result[i][j] = tf[i][j];
      break;
    }
    case Linear: {
      auto tf = transfactor<Linear, dim>(c);
      for (size_t i = 0; i < nbasis; i++)
        for (size_t j = 0; j < nbasis; j++)
          result[i][j] = tf[i][j];
      break;
    }
    case Quadratic: {
      auto tf = transfactor<Quadratic, dim>(c);
      for (size_t i = 0; i < nbasis; i++)
        for (size_t j = 0; j < nbasis; j++)
          result[i][j] = tf[i][j];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Unitary
////////////////////////////////////////////////////////////////////////////////

template<>
inline
std::array<double, Traits<Unitary, 1>::function_size>
function<Unitary, 1>(Wonton::Point<1> x) {
  std::array<double, Traits<Unitary, 1>::function_size> r {};
  r[0] = 1.;
  return r;
}

template<>
inline
std::array<double, Traits<Unitary, 2>::function_size>
function<Unitary, 2>(Wonton::Point<2> x) {
  std::array<double, Traits<Unitary, 2>::function_size> r {};
  r[0] = 1.;
  return r;
}

template<>
inline
std::array<double, Traits<Unitary, 3>::function_size>
function<Unitary, 3>(Wonton::Point<3> x) {
  std::array<double, Traits<Unitary, 3>::function_size> r {};
  r[0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 1>::function_size>,
  Traits<Unitary, 1>::function_size
>
jet<Unitary, 1>(Wonton::Point<1> x) {
  std::array<
    std::array<double, Traits<Unitary, 1>::function_size>,
    Traits<Unitary, 1>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 2>::function_size>,
  Traits<Unitary, 2>::function_size
>
jet<Unitary, 2>(Wonton::Point<2> x) {
  std::array<
    std::array<double, Traits<Unitary, 2>::function_size>,
    Traits<Unitary, 2>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 3>::function_size>,
  Traits<Unitary, 3>::function_size
>
jet<Unitary, 3>(Wonton::Point<3> x) {
  std::array<
    std::array<double, Traits<Unitary, 3>::function_size>,
    Traits<Unitary, 3>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 1>::function_size>,
  Traits<Unitary, 1>::function_size
>
inverse_jet<Unitary, 1>(Wonton::Point<1> x) {
  std::array<
    std::array<double, Traits<Unitary, 1>::function_size>,
    Traits<Unitary, 1>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 2>::function_size>,
  Traits<Unitary, 2>::function_size
>
inverse_jet<Unitary, 2>(Wonton::Point<2> x) {
  std::array<
    std::array<double, Traits<Unitary, 2>::function_size>,
    Traits<Unitary, 2>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Unitary, 3>::function_size>,
  Traits<Unitary, 3>::function_size
>
inverse_jet<Unitary, 3>(Wonton::Point<3> x) {
  std::array<
    std::array<double, Traits<Unitary, 3>::function_size>,
    Traits<Unitary, 3>::function_size
  > r {};
  r.fill({}); // zero-init
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<double, Traits<Unitary, 1>::function_size>
shift<Unitary, 1>(Wonton::Point<1> x, Wonton::Point<1> y) {
  std::array<double, Traits<Unitary, 1>::function_size> r {};
  r.fill({}); // zero-init
  r[0] = 1.;
  return r;
}

template<>
inline
std::array<double, Traits<Unitary, 2>::function_size>
shift<Unitary, 2>(Wonton::Point<2> x, Wonton::Point<2> y) {
  std::array<double, Traits<Unitary, 2>::function_size> r {};
  r.fill({}); // zero-init
  r[0] = 1.;
  return r;
}

template<>
inline
std::array<double, Traits<Unitary, 3>::function_size>
shift<Unitary, 3>(Wonton::Point<3> x, Wonton::Point<3> y) {
  std::array<double, Traits<Unitary, 3>::function_size> r {};
  r.fill({}); // zero-init
  r[0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
typename Traits<Unitary, 1>::matrix_t transfactor<Unitary, 1>(const Wonton::Point<1>& c) {
  typename Traits<Unitary, 1>::matrix_t tf;
  tf[0][0] = 1.;
  return tf;
}

template<>
inline
typename Traits<Unitary, 2>::matrix_t transfactor<Unitary, 2>(const Wonton::Point<2>& c) {
  typename Traits<Unitary, 2>::matrix_t tf;
  tf[0][0] = 1.;
  return tf;
}

template<>
inline
typename Traits<Unitary, 3>::matrix_t transfactor<Unitary, 3>(const Wonton::Point<3>& c) {
  typename Traits<Unitary, 3>::matrix_t tf;
  tf[0][0] = 1.;
  return tf;
}

////////////////////////////////////////////////////////////////////////////////
// Linear
////////////////////////////////////////////////////////////////////////////////

template<>
inline
std::array<double, Traits<Linear, 1>::function_size>
function<Linear, 1>(Wonton::Point<1> x) {
  std::array<double, Traits<Linear, 1>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  return r;
}

template<>
inline
std::array<double, Traits<Linear, 2>::function_size>
function<Linear, 2>(Wonton::Point<2> x) {
  std::array<double, Traits<Linear, 2>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  return r;
}

template<>
inline
std::array<double, Traits<Linear, 3>::function_size>
function<Linear, 3>(Wonton::Point<3> x) {
  std::array<double, Traits<Linear, 3>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = x[2];
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Linear, 1>::function_size>,
  Traits<Linear, 1>::function_size
>
jet<Linear, 1>(Wonton::Point<1> x) {
  std::array<
    std::array<double, Traits<Linear, 1>::function_size>,
    Traits<Linear, 1>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Linear, 2>::function_size>,
  Traits<Linear, 2>::function_size
>
jet<Linear, 2>(Wonton::Point<2> x) {
  std::array<
    std::array<double, Traits<Linear, 2>::function_size>,
    Traits<Linear, 2>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[2][0] = x[1];
  r[0][1] = 0.;
  r[1][1] = 1.;
  r[2][1] = 0.;
  r[0][2] = 0.;
  r[1][2] = 0.;
  r[2][2] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Linear, 3>::function_size>,
  Traits<Linear, 3>::function_size
>
jet<Linear, 3>(Wonton::Point<3> x) {
  std::array<
    std::array<double, Traits<Linear, 3>::function_size>,
    Traits<Linear, 3>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[2][0] = x[1];
  r[3][0] = x[2];
  r[0][1] = 0.;
  r[1][1] = 1.;
  r[2][1] = 0.;
  r[3][1] = 0.;
  r[0][2] = 0.;
  r[1][2] = 0.;
  r[2][2] = 1.;
  r[3][2] = 0.;
  r[0][3] = 0.;
  r[1][3] = 0.;
  r[2][3] = 0.;
  r[3][3] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Linear, 1>::function_size>,
  Traits<Linear, 1>::function_size
>
inverse_jet<Linear, 1>(Wonton::Point<1> x) {
  std::array<
    std::array<double, Traits<Linear, 1>::function_size>,
    Traits<Linear, 1>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = -x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Linear, 2>::function_size>,
  Traits<Linear, 2>::function_size
>
inverse_jet<Linear, 2>(Wonton::Point<2> x) {
  std::array<
    std::array<double, Traits<Linear, 2>::function_size>,
    Traits<Linear, 2>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = -x[0];
  r[2][0] = -x[1];
  r[0][1] = 0.;
  r[1][1] = 1.;
  r[2][1] = 0.;
  r[0][2] = 0.;
  r[1][2] = 0.;
  r[2][2] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Linear, 3>::function_size>,
  Traits<Linear, 3>::function_size
>
inverse_jet<Linear, 3>(Wonton::Point<3> x) {
  std::array<
    std::array<double, Traits<Linear, 3>::function_size>,
    Traits<Linear, 3>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = -x[0];
  r[2][0] = -x[1];
  r[3][0] = -x[2];
  r[0][1] = 0.;
  r[1][1] = 1.;
  r[2][1] = 0.;
  r[3][1] = 0.;
  r[0][2] = 0.;
  r[1][2] = 0.;
  r[2][2] = 1.;
  r[3][2] = 0.;
  r[0][3] = 0.;
  r[1][3] = 0.;
  r[2][3] = 0.;
  r[3][3] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<double, Traits<Linear, 1>::function_size>
shift<Linear, 1>(Wonton::Point<1> x, Wonton::Point<1> y) {
  std::array<double, Traits<Linear, 1>::function_size> r {};
  r[0] = 1.;
  r[1] = y[0] - x[0];
  return r;
}

template<>
inline
std::array<double, Traits<Linear, 2>::function_size>
shift<Linear, 2>(Wonton::Point<2> x, Wonton::Point<2> y) {
  std::array<double, Traits<Linear, 2>::function_size> r {};
  r[0] = 1.;
  r[1] = y[0] - x[0];
  r[2] = y[1] - x[1];
  return r;
}

template<>
inline
std::array<double, Traits<Linear, 3>::function_size>
shift<Linear, 3>(Wonton::Point<3> x, Wonton::Point<3> y) {
  std::array<double, Traits<Linear, 3>::function_size> r {};
  r[0] = 1.;
  r[1] = y[0] - x[0];
  r[2] = y[1] - x[1];
  r[3] = y[2] - x[2];
  return r;
}

//---------------------------------------------------------------

template<>
inline
typename Traits<Linear, 1>::matrix_t transfactor<Linear, 1>(const Wonton::Point<1>& c) {
  typename Traits<Linear, 1>::matrix_t tf;
  tf[0][0] = 1;
  tf[0][1] = 0;
  tf[1][0] = c[0];
  tf[1][1] = 1;
  return tf;
}

template<>
inline
typename Traits<Linear, 2>::matrix_t transfactor<Linear, 2>(const Wonton::Point<2>& c) {
  typename Traits<Linear, 2>::matrix_t tf;
  tf[0][0] = 1;
  tf[0][1] = 0;
  tf[0][2] = 0;
  tf[1][0] = c[0];
  tf[1][1] = 1;
  tf[1][2] = 0;
  tf[2][0] = c[1];
  tf[2][1] = 0;
  tf[2][2] = 1;
  return tf;
}

template<>
inline
typename Traits<Linear, 3>::matrix_t transfactor<Linear, 3>(const Wonton::Point<3>& c) {
  typename Traits<Linear, 3>::matrix_t tf;
  tf[0][0] = 1;
  tf[0][1] = 0;
  tf[0][2] = 0;
  tf[0][3] = 0;
  tf[1][0] = c[0];
  tf[1][1] = 1;
  tf[1][2] = 0;
  tf[1][3] = 0;
  tf[2][0] = c[1];
  tf[2][1] = 0;
  tf[2][2] = 1;
  tf[2][3] = 0;
  tf[3][0] = c[2];
  tf[3][1] = 0;
  tf[3][2] = 0;
  tf[3][3] = 1;
  return tf;
}

////////////////////////////////////////////////////////////////////////////////
// Quadratic
////////////////////////////////////////////////////////////////////////////////

template<>
inline
std::array<double, Traits<Quadratic, 1>::function_size>
function<Quadratic, 1>(Wonton::Point<1> x) {
  std::array<double, Traits<Quadratic, 1>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = .5 * x[0] * x[0];
  return r;
}

template<>
inline
std::array<double, Traits<Quadratic, 2>::function_size>
function<Quadratic, 2>(Wonton::Point<2> x) {
  std::array<double, Traits<Quadratic, 2>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = .5 * x[0] * x[0];
  r[4] = x[0] * x[1];
  r[5] = .5 * x[1] * x[1];
  return r;
}

template<>
inline
std::array<double, Traits<Quadratic, 3>::function_size>
function<Quadratic, 3>(Wonton::Point<3> x) {
  std::array<double, Traits<Quadratic, 3>::function_size> r {};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = x[2];
  r[4] = .5 * x[0] * x[0];
  r[5] = x[0] * x[1];
  r[6] = x[0] * x[2];
  r[7] = .5 * x[1] * x[1];
  r[8] = x[1] * x[2];
  r[9] = .5 * x[2] * x[2];
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 1>::function_size>,
  Traits<Quadratic, 1>::function_size
>
jet<Quadratic, 1>(Wonton::Point<1> x) {
  std::array<
    std::array<double, Traits<Quadratic, 1>::function_size>,
    Traits<Quadratic, 1>::function_size
  > r {};
  r[0][0] = 1.;

  r[1][0] = x[0];
  r[1][1] = 1.;

  r[2][0] = x[0] * x[0] * .5;
  r[2][1] = x[0];
  r[2][2] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 2>::function_size>,
  Traits<Quadratic, 2>::function_size
>
jet<Quadratic, 2>(Wonton::Point<2> x) {
  std::array<
    std::array<double, Traits<Quadratic, 2>::function_size>,
    Traits<Quadratic, 2>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[2][0] = x[1];
  r[3][0] = .5 * x[0] * x[0];
  r[4][0] = x[0] * x[1];
  r[5][0] = .5 * x[1] * x[1];

  r[1][1] = 1.;
  r[3][1] = x[0];
  r[4][1] = x[1];

  r[2][2] = 1.;
  r[4][2] = x[0];
  r[5][2] = x[1];

  for (size_t i = 3; i < 6; i++)
    r[i][i] = 1.;
  return r;
}

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 3>::function_size>,
  Traits<Quadratic, 3>::function_size
>
jet<Quadratic, 3>(Wonton::Point<3> x) {
  std::array<
    std::array<double, Traits<Quadratic, 3>::function_size>,
    Traits<Quadratic, 3>::function_size
  > r {};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[2][0] = x[1];
  r[3][0] = x[2];
  r[4][0] = .5 * x[0] * x[0];
  r[5][0] = x[0] * x[1];
  r[6][0] = x[0] * x[2];
  r[7][0] = .5 * x[1] * x[1];
  r[8][0] = x[1] * x[2];
  r[9][0] = .5 * x[2] * x[2];

  r[1][1] = 1.;
  r[4][1] = x[0];
  r[5][1] = x[1];
  r[6][1] = x[2];

  r[2][2] = 1.;
  r[5][2] = x[0];
  r[7][2] = x[1];
  r[8][2] = x[2];

  r[3][3] = 1.;
  r[6][3] = x[0];
  r[8][3] = x[1];
  r[9][3] = x[2];

  for (size_t i = 4; i < 10; i++)
    r[i][i] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 1>::function_size>,
  Traits<Quadratic, 1>::function_size
>
inverse_jet<Quadratic, 1>(Wonton::Point<1> x) {
  Wonton::Point<1> mx;
  for (size_t i = 0; i < 1; i++) mx[i] = -x[i];
  return jet<Quadratic, 1>(mx);
}

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 2>::function_size>,
  Traits<Quadratic, 2>::function_size
>
inverse_jet<Quadratic, 2>(Wonton::Point<2> x) {
  Wonton::Point<2> mx;
  for (size_t i = 0; i < 2; i++) mx[i] = -x[i];
  return jet<Quadratic, 2>(mx);
}

template<>
inline
std::array<
  std::array<double, Traits<Quadratic, 3>::function_size>,
  Traits<Quadratic, 3>::function_size
>
inverse_jet<Quadratic, 3>(Wonton::Point<3> x) {
  Wonton::Point<3> mx;
  for (size_t i = 0; i < 3; i++) mx[i] = -x[i];
  return jet<Quadratic, 3>(mx);
}

//---------------------------------------------------------------

template<>
inline
std::array<double, Traits<Quadratic, 1>::function_size>
shift<Quadratic, 1>(Wonton::Point<1> x, Wonton::Point<1> y) {
  Wonton::Point<1> d;
  for (size_t i = 0; i < 1; i++) d[i] = y[i] - x[i];
  return function<Quadratic, 1>(d);
}

template<>
inline
std::array<double, Traits<Quadratic, 2>::function_size>
shift<Quadratic, 2>(Wonton::Point<2> x, Wonton::Point<2> y) {
  Wonton::Point<2> d;
  for (size_t i = 0; i < 2; i++) d[i] = y[i] - x[i];
  return function<Quadratic, 2>(d);
}

template<>
inline
std::array<double, Traits<Quadratic, 3>::function_size>
shift<Quadratic, 3>(Wonton::Point<3> x, Wonton::Point<3> y) {
  Wonton::Point<3> d;
  for (size_t i = 0; i < 3; i++) d[i] = y[i] - x[i];
  return function<Quadratic, 3>(d);
}

//---------------------------------------------------------------

template<>
inline
typename Traits<Quadratic, 1>::matrix_t transfactor<Quadratic, 1>(const Wonton::Point<1>& c) {
  typename Traits<Quadratic, 1>::matrix_t tf;
  tf[0] = { 1, 0, 0 };
  tf[1] = { c[0], 1, 0 };
  tf[2] = { 0.5 * c[0] * c[0], c[0], 1 };
  return tf;
}

template<>
inline
typename Traits<Quadratic, 2>::matrix_t transfactor<Quadratic, 2>(const Wonton::Point<2>& c) {
  typename Traits<Quadratic, 2>::matrix_t tf;
  tf[0] = { 1, 0, 0, 0, 0, 0 };
  tf[1] = { c[0], 1, 0, 0, 0, 0 };
  tf[2] = { c[1], 0, 1, 0, 0, 0 };
  tf[3] = { 0.5 * c[0] * c[0], c[0], 0, 1, 0, 0 };
  tf[4] = { c[0] * c[1], c[1], c[0], 0, 1, 0 };
  tf[5] = { 0.5 * c[1] * c[1], 0, c[1], 0, 0, 1 };
  return tf;
}

template<>
inline
typename Traits<Quadratic, 3>::matrix_t transfactor<Quadratic, 3>(const Wonton::Point<3>& c) {
  typename Traits<Quadratic, 3>::matrix_t tf;
  tf[0] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  tf[1] = { c[0], 1, 0, 0, 0, 0, 0, 0, 0, 0 };
  tf[2] = { c[1], 0, 1, 0, 0, 0, 0, 0, 0, 0 };
  tf[3] = { c[2], 0, 0, 1, 0, 0, 0, 0, 0, 0 };
  tf[4] = { 0.5 * c[0] * c[0], c[0], 0, 0, 1, 0, 0, 0, 0, 0 };
  tf[5] = { c[0] * c[1], c[1], c[0], 0, 0, 1, 0, 0, 0, 0 };
  tf[6] = { c[0] * c[2], c[2], 0, c[0], 0, 0, 1, 0, 0, 0 };
  tf[7] = { 0.5 * c[1] * c[1], 0, c[1], 0, 0, 0, 0, 1, 0, 0 };
  tf[8] = { c[1] * c[2], 0, c[2], c[1], 0, 0, 0, 0, 1, 0 };
  tf[9] = { 0.5 * c[2] * c[2], 0, 0, c[2], 0, 0, 0, 0, 0, 1 };
  return tf;
}

}}}  // namespace Portage::Meshfree::Basis

#endif  // PORTAGE_SUPPORT_BASIS_H_
