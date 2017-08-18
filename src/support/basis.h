/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#define BASIS_H_INC_

#include <cassert>
#include <cmath>
#include <array>
#include <vector>

#include "portage/support/Point.h"

namespace Portage {
namespace Meshfree {
namespace Basis {

using std::array;
using std::vector;

enum Type {Unitary, Linear, Quadratic, LastBasis};

constexpr size_t quadratic_sizes[4]={0,3,6,10};

////////////////////////////////////////////////////////////////////////////////
// Traits
////////////////////////////////////////////////////////////////////////////////

template<Type type, size_t dim>
class Traits{};

template<size_t dim>
class Traits<Unitary, dim>
{
  public:
  static constexpr size_t function_size=1;
  static constexpr array<size_t,2> jet_size={1,1};
};

template<size_t dim>
class Traits<Linear, dim>
{
  public:
  static constexpr size_t function_size=dim+1;
  static constexpr array<size_t,2> jet_size={dim+1,dim+1};
};

template<size_t dim>
class Traits<Quadratic, dim>
{
  public:
  static constexpr size_t function_size=quadratic_sizes[dim];
  static constexpr array<size_t,2>
    jet_size={quadratic_sizes[dim], quadratic_sizes[dim]};
};

template<size_t dim>
size_t function_size(Type btype){
  size_t result;
  switch (btype) {
    case Unitary:
      result = Traits<Unitary, dim>::function_size;
      break;
    case Linear:
      result = Traits<Linear, dim>::function_size;
      break;
    case Quadratic:
      result = Traits<Quadratic, dim>::function_size;
      break;
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
array<size_t,2> jet_size(Type btype) {
  array<size_t,2> result;
  switch(btype) {
    case Unitary:
      result = Traits<Unitary, dim>::jet_size;
      break;
    case Linear:
      result = Traits<Linear, dim>::jet_size;
      break;
    case Quadratic:
      result = Traits<Quadratic, dim>::jet_size;
      break;
    default:
      assert(false);
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Templates
////////////////////////////////////////////////////////////////////////////////

template<Type type, size_t dim>
array<double, Traits<type, dim>::function_size>
function(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<
array<double, Traits<type, dim>::function_size>,
Traits<type, dim>::function_size
>
jet(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<
array<double, Traits<type, dim>::function_size>,
Traits<type, dim>::function_size
>
inverse_jet(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<double, Traits<type, dim>::function_size>
shift(Point<dim> x, Point<dim> y){
  auto ijet = inverse_jet<type,dim>(x);
  auto by = function<type,dim>(y);
  array<double, Traits<type, dim>::function_size> r{0.0};

  size_t n = Traits<type, dim>::function_size;
  for (size_t i=0; i<n; i++) r[i] = 0.;
  for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) {
    r[i] += ijet[i][j]*by[j];
  }

  return r;
}

template<size_t dim>
vector<double> function(Type type, Point<dim> x){
  size_t nbasis = function_size<dim>(type);
  vector<double> result(nbasis);
  switch (type) {
    case Unitary: {
      auto resulta = function<Unitary, dim>(x);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    case Linear: {
      auto resulta = function<Linear, dim>(x);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    case Quadratic: {
      auto resulta = function<Quadratic, dim>(x);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
vector<double> shift(Type type, Point<dim> x, Point<dim> y){
  size_t nbasis = function_size<dim>(type);
  vector<double> result(nbasis);
  switch (type) {
    case Unitary: {
      auto resulta = shift<Unitary, dim>(x,y);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    case Linear: {
      auto resulta = shift<Linear, dim>(x,y);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    case Quadratic: {
      auto resulta = shift<Quadratic, dim>(x,y);
      for (size_t i=0; i<nbasis; i++) result[i] = resulta[i];
      break;
    }
    default:
      assert(false);
  }
  return result;
}

template<size_t dim>
vector<vector<double>> jet(Type type, Point<dim> x){
  auto njet = jet_size<dim>(type);
  vector<vector<double>> result(njet[0],vector<double>(njet[1],0.));
  switch (type) {
    case Unitary: {
      auto resulta = jet<Unitary, dim>(x);
      for (size_t i=0; i<njet[0]; i++) for (size_t j=0; j<njet[1]; j++)
        result[i][j] = resulta[i][j];
      break;
    }
    case Linear: {
      auto resulta = jet<Linear, dim>(x);
      for (size_t i=0; i<njet[0]; i++) for (size_t j=0; j<njet[1]; j++)
        result[i][j] = resulta[i][j];
      break;
    }
    case Quadratic: {
      auto resulta = jet<Quadratic, dim>(x);
      for (size_t i=0; i<njet[0]; i++) for (size_t j=0; j<njet[1]; j++)
        result[i][j] = resulta[i][j];
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
array<double, Traits<Unitary, 1>::function_size>
function<Unitary, 1>(Point<1> x) {
  array<double, Traits<Unitary, 1>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, Traits<Unitary, 2>::function_size>
function<Unitary, 2>(Point<2> x) {
  array<double, Traits<Unitary, 2>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, Traits<Unitary, 3>::function_size>
function<Unitary, 3>(Point<3> x) {
  array<double, Traits<Unitary, 3>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, Traits<Unitary, 1>::function_size>,
Traits<Unitary, 1>::function_size
>
jet<Unitary,1>(Point<1> x){
  array<
  array<double, Traits<Unitary, 1>::function_size>,
  Traits<Unitary, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Unitary, 2>::function_size>,
Traits<Unitary, 2>::function_size
>
jet<Unitary,2>(Point<2> x){
  array<
  array<double, Traits<Unitary, 2>::function_size>,
  Traits<Unitary, 2>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Unitary, 3>::function_size>,
Traits<Unitary, 3>::function_size
>
jet<Unitary,3>(Point<3> x){
  array<
  array<double, Traits<Unitary, 3>::function_size>,
  Traits<Unitary, 3>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, Traits<Unitary, 1>::function_size>,
Traits<Unitary, 1>::function_size
>
inverse_jet<Unitary,1>(Point<1> x){
  array<
  array<double, Traits<Unitary, 1>::function_size>,
  Traits<Unitary, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Unitary, 2>::function_size>,
Traits<Unitary, 2>::function_size
>
inverse_jet<Unitary,2>(Point<2> x){
  array<
  array<double, Traits<Unitary, 2>::function_size>,
  Traits<Unitary, 2>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Unitary, 3>::function_size>,
Traits<Unitary, 3>::function_size
>
inverse_jet<Unitary,3>(Point<3> x){
  array<
  array<double, Traits<Unitary, 3>::function_size>,
  Traits<Unitary, 3>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<double, Traits<Unitary, 1>::function_size>
shift<Unitary, 1>(Point<1> x, Point<1> y) {
  array<double, Traits<Unitary, 1>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, Traits<Unitary, 2>::function_size>
shift<Unitary, 2>(Point<2> x, Point<2> y) {
  array<double, Traits<Unitary, 2>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, Traits<Unitary, 3>::function_size>
shift<Unitary, 3>(Point<3> x, Point<3> y) {
  array<double, Traits<Unitary, 3>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Linear
////////////////////////////////////////////////////////////////////////////////

template<>
array<double, Traits<Linear, 1>::function_size>
function<Linear, 1>(Point<1> x) {
  array<double, Traits<Linear, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  return r;
}

template<>
array<double, Traits<Linear, 2>::function_size>
function<Linear, 2>(Point<2> x) {
  array<double, Traits<Linear, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  return r;
}

template<>
array<double, Traits<Linear, 3>::function_size>
function<Linear, 3>(Point<3> x) {
  array<double, Traits<Linear, 3>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = x[2];
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, Traits<Linear, 1>::function_size>,
Traits<Linear, 1>::function_size
>
jet<Linear,1>(Point<1> x){
  array<
  array<double, Traits<Linear, 1>::function_size>,
  Traits<Linear, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Linear, 2>::function_size>,
Traits<Linear, 2>::function_size
>
jet<Linear,2>(Point<2> x){
  array<
  array<double, Traits<Linear, 2>::function_size>,
  Traits<Linear, 2>::function_size
    > r{0.0};
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
array<
array<double, Traits<Linear, 3>::function_size>,
Traits<Linear, 3>::function_size
>
jet<Linear,3>(Point<3> x){
  array<
  array<double, Traits<Linear, 3>::function_size>,
  Traits<Linear, 3>::function_size
    > r{0.0};
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
array<
array<double, Traits<Linear, 1>::function_size>,
Traits<Linear, 1>::function_size
>
inverse_jet<Linear,1>(Point<1> x){
  array<
  array<double, Traits<Linear, 1>::function_size>,
  Traits<Linear, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  r[1][0] = -x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Linear, 2>::function_size>,
Traits<Linear, 2>::function_size
>
inverse_jet<Linear,2>(Point<2> x){
  array<
  array<double, Traits<Linear, 2>::function_size>,
  Traits<Linear, 2>::function_size
    > r{0.0};
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
array<
array<double, Traits<Linear, 3>::function_size>,
Traits<Linear, 3>::function_size
>
inverse_jet<Linear,3>(Point<3> x){
  array<
  array<double, Traits<Linear, 3>::function_size>,
  Traits<Linear, 3>::function_size
    > r{0.0};
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
array<double, Traits<Linear, 1>::function_size>
shift<Linear, 1>(Point<1> x, Point<1> y) {
  array<double, Traits<Linear, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = y[0]-x[0];
  return r;
}

template<>
array<double, Traits<Linear, 2>::function_size>
shift<Linear, 2>(Point<2> x, Point<2> y) {
  array<double, Traits<Linear, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = y[0]-x[0];
  r[2] = y[1]-x[1];
  return r;
}

template<>
array<double, Traits<Linear, 3>::function_size>
shift<Linear, 3>(Point<3> x, Point<3> y) {
  array<double, Traits<Linear, 3>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = y[0]-x[0];
  r[2] = y[1]-x[1];
  r[3] = y[2]-x[2];
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Quadratic
////////////////////////////////////////////////////////////////////////////////

template<>
array<double, Traits<Quadratic, 1>::function_size>
function<Quadratic, 1>(Point<1> x) {
  array<double, Traits<Quadratic, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = .5 * x[0] * x[0];
  return r;
}

template<>
array<double, Traits<Quadratic, 2>::function_size>
function<Quadratic, 2>(Point<2> x) {
  array<double, Traits<Quadratic, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = .5 * x[0] * x[0];
  r[4] = x[0] * x[1];
  r[5] = .5 * x[1] * x[1];
  return r;
}

template<>
array<double, Traits<Quadratic, 3>::function_size>
function<Quadratic, 3>(Point<3> x) {
  array<double, Traits<Quadratic, 3>::function_size> r{0.0};
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
array<
array<double, Traits<Quadratic, 1>::function_size>,
Traits<Quadratic, 1>::function_size
>
jet<Quadratic,1>(Point<1> x){
  array<
  array<double, Traits<Quadratic, 1>::function_size>,
  Traits<Quadratic, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;

  r[1][0] = x[0];
  r[1][1] = 1.;

  r[2][0] = x[0] * x[0] * .5;
  r[2][1] = x[0];
  r[2][2] = 1.;
  return r;
}

template<>
array<
array<double, Traits<Quadratic, 2>::function_size>,
Traits<Quadratic, 2>::function_size
>
jet<Quadratic,2>(Point<2> x){
  array<
  array<double, Traits<Quadratic, 2>::function_size>,
  Traits<Quadratic, 2>::function_size
    > r{0.0};
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
array<
array<double, Traits<Quadratic, 3>::function_size>,
Traits<Quadratic, 3>::function_size
>
jet<Quadratic,3>(Point<3> x){
  array<
  array<double, Traits<Quadratic, 3>::function_size>,
  Traits<Quadratic, 3>::function_size
    > r{0.0};
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
array<
array<double, Traits<Quadratic, 1>::function_size>,
Traits<Quadratic, 1>::function_size
>
inverse_jet<Quadratic,1>(Point<1> x){
  Point<1> mx;
  for (size_t i=0; i<1; i++) mx[i] = -x[i];
  return jet<Quadratic,1>(mx);
}

template<>
array<
array<double, Traits<Quadratic, 2>::function_size>,
Traits<Quadratic, 2>::function_size
>
inverse_jet<Quadratic,2>(Point<2> x){
  Point<2> mx;
  for (size_t i=0; i<2; i++) mx[i] = -x[i];
  return jet<Quadratic,2>(mx);
}

template<>
array<
array<double, Traits<Quadratic, 3>::function_size>,
Traits<Quadratic, 3>::function_size
>
inverse_jet<Quadratic,3>(Point<3> x){
  Point<3> mx;
  for (size_t i=0; i<3; i++) mx[i] = -x[i];
  return jet<Quadratic,3>(mx);
}

//---------------------------------------------------------------

template<>
array<double, Traits<Quadratic, 1>::function_size>
shift<Quadratic, 1>(Point<1> x, Point<1> y) {
  Point<1> d;
  for (size_t i=0; i<1; i++) d[i] = y[i]-x[i];
  return function<Quadratic,1>(d);
}

template<>
array<double, Traits<Quadratic, 2>::function_size>
shift<Quadratic, 2>(Point<2> x, Point<2> y) {
  Point<2> d;
  for (size_t i=0; i<2; i++) d[i] = y[i]-x[i];
  return function<Quadratic,2>(d);
}

template<>
array<double, Traits<Quadratic, 3>::function_size>
shift<Quadratic, 3>(Point<3> x, Point<3> y) {
  Point<3> d;
  for (size_t i=0; i<3; i++) d[i] = y[i]-x[i];
  return function<Quadratic,3>(d);
}

}
}
}

#endif // BASIS_H_INC
