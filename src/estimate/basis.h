/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef BASIS_H_INC_
#define BASIS_H_INC_

#include <cassert>
#include <cmath>
#include <array>

#include "portage/support/Point.h"

namespace Portage {
namespace Meshfree {
namespace Basis {

using std::array;

enum Type {Unitary, Linear, Quadratic, LastBasis};

constexpr size_t quadratic_sizes[4]={0,3,6,10};

////////////////////////////////////////////////////////////////////////////////
// Traits
////////////////////////////////////////////////////////////////////////////////

template<Type type, size_t dim>
class BasisTraits{};

template<size_t dim>
class BasisTraits<Unitary, dim>
{
  public:
  static constexpr size_t function_size=1;
  static constexpr array<size_t,2> jet_size={1,1};
};

template<size_t dim>
class BasisTraits<Linear, dim>
{
  public:
  static constexpr size_t function_size=dim+1;
  static constexpr array<size_t,2> jet_size={dim+1,dim+1};
};

template<size_t dim>
class BasisTraits<Quadratic, dim>
{
  public:
  static constexpr size_t function_size=quadratic_sizes[dim];
  static constexpr array<size_t,2>
    jet_size={quadratic_sizes[dim], quadratic_sizes[dim]};
};

////////////////////////////////////////////////////////////////////////////////
// Templates
////////////////////////////////////////////////////////////////////////////////

template<Type type, size_t dim>
array<double, BasisTraits<type, dim>::function_size>
basis_function(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<
array<double, BasisTraits<type, dim>::function_size>,
BasisTraits<type, dim>::function_size
>
basis_jet(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<
array<double, BasisTraits<type, dim>::function_size>,
BasisTraits<type, dim>::function_size
>
basis_inverse_jet(Point<dim> x){assert(false);}

template<Type type, size_t dim>
array<double, BasisTraits<type, dim>::function_size>
basis_shift(Point<dim> x, Point<dim> y){
  auto ijet = basis_inverse_jet<type,dim>(x);
  auto by = basis_function<type,dim>(y);
  array<double, BasisTraits<type, dim>::function_size> r{0.0};

  size_t n = BasisTraits<type, dim>::function_size;
  for (size_t i=0; i<n; i++) r[i] = 0.;
  for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) {
    r[i] += ijet[i][j]*by[j];
  }

  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Unitary
////////////////////////////////////////////////////////////////////////////////

template<>
array<double, BasisTraits<Unitary, 1>::function_size>
basis_function<Unitary, 1>(Point<1> x) {
  array<double, BasisTraits<Unitary, 1>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, BasisTraits<Unitary, 2>::function_size>
basis_function<Unitary, 2>(Point<2> x) {
  array<double, BasisTraits<Unitary, 2>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, BasisTraits<Unitary, 3>::function_size>
basis_function<Unitary, 3>(Point<3> x) {
  array<double, BasisTraits<Unitary, 3>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, BasisTraits<Unitary, 1>::function_size>,
BasisTraits<Unitary, 1>::function_size
>
basis_jet<Unitary,1>(Point<1> x){
  array<
  array<double, BasisTraits<Unitary, 1>::function_size>,
  BasisTraits<Unitary, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Unitary, 2>::function_size>,
BasisTraits<Unitary, 2>::function_size
>
basis_jet<Unitary,2>(Point<2> x){
  array<
  array<double, BasisTraits<Unitary, 2>::function_size>,
  BasisTraits<Unitary, 2>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Unitary, 3>::function_size>,
BasisTraits<Unitary, 3>::function_size
>
basis_jet<Unitary,3>(Point<3> x){
  array<
  array<double, BasisTraits<Unitary, 3>::function_size>,
  BasisTraits<Unitary, 3>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, BasisTraits<Unitary, 1>::function_size>,
BasisTraits<Unitary, 1>::function_size
>
basis_inverse_jet<Unitary,1>(Point<1> x){
  array<
  array<double, BasisTraits<Unitary, 1>::function_size>,
  BasisTraits<Unitary, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Unitary, 2>::function_size>,
BasisTraits<Unitary, 2>::function_size
>
basis_inverse_jet<Unitary,2>(Point<2> x){
  array<
  array<double, BasisTraits<Unitary, 2>::function_size>,
  BasisTraits<Unitary, 2>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Unitary, 3>::function_size>,
BasisTraits<Unitary, 3>::function_size
>
basis_inverse_jet<Unitary,3>(Point<3> x){
  array<
  array<double, BasisTraits<Unitary, 3>::function_size>,
  BasisTraits<Unitary, 3>::function_size
    > r{0.0};
  r[0][0] = 1.;
  return r;
}

//---------------------------------------------------------------

template<>
array<double, BasisTraits<Unitary, 1>::function_size>
basis_shift<Unitary, 1>(Point<1> x, Point<1> y) {
  array<double, BasisTraits<Unitary, 1>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, BasisTraits<Unitary, 2>::function_size>
basis_shift<Unitary, 2>(Point<2> x, Point<2> y) {
  array<double, BasisTraits<Unitary, 2>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

template<>
array<double, BasisTraits<Unitary, 3>::function_size>
basis_shift<Unitary, 3>(Point<3> x, Point<3> y) {
  array<double, BasisTraits<Unitary, 3>::function_size> r{0.0};
  r[0] = 1.;
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Linear
////////////////////////////////////////////////////////////////////////////////

template<>
array<double, BasisTraits<Linear, 1>::function_size>
basis_function<Linear, 1>(Point<1> x) {
  array<double, BasisTraits<Linear, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  return r;
}

template<>
array<double, BasisTraits<Linear, 2>::function_size>
basis_function<Linear, 2>(Point<2> x) {
  array<double, BasisTraits<Linear, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  return r;
}

template<>
array<double, BasisTraits<Linear, 3>::function_size>
basis_function<Linear, 3>(Point<3> x) {
  array<double, BasisTraits<Linear, 3>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = x[2];
  return r;
}

//---------------------------------------------------------------

template<>
array<
array<double, BasisTraits<Linear, 1>::function_size>,
BasisTraits<Linear, 1>::function_size
>
basis_jet<Linear,1>(Point<1> x){
  array<
  array<double, BasisTraits<Linear, 1>::function_size>,
  BasisTraits<Linear, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  r[1][0] = x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Linear, 2>::function_size>,
BasisTraits<Linear, 2>::function_size
>
basis_jet<Linear,2>(Point<2> x){
  array<
  array<double, BasisTraits<Linear, 2>::function_size>,
  BasisTraits<Linear, 2>::function_size
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
array<double, BasisTraits<Linear, 3>::function_size>,
BasisTraits<Linear, 3>::function_size
>
basis_jet<Linear,3>(Point<3> x){
  array<
  array<double, BasisTraits<Linear, 3>::function_size>,
  BasisTraits<Linear, 3>::function_size
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
array<double, BasisTraits<Linear, 1>::function_size>,
BasisTraits<Linear, 1>::function_size
>
basis_inverse_jet<Linear,1>(Point<1> x){
  array<
  array<double, BasisTraits<Linear, 1>::function_size>,
  BasisTraits<Linear, 1>::function_size
    > r{0.0};
  r[0][0] = 1.;
  r[1][0] = -x[0];
  r[0][1] = 0.;
  r[1][1] = 1.;
  return r;
}

template<>
array<
array<double, BasisTraits<Linear, 2>::function_size>,
BasisTraits<Linear, 2>::function_size
>
basis_inverse_jet<Linear,2>(Point<2> x){
  array<
  array<double, BasisTraits<Linear, 2>::function_size>,
  BasisTraits<Linear, 2>::function_size
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
array<double, BasisTraits<Linear, 3>::function_size>,
BasisTraits<Linear, 3>::function_size
>
basis_inverse_jet<Linear,3>(Point<3> x){
  array<
  array<double, BasisTraits<Linear, 3>::function_size>,
  BasisTraits<Linear, 3>::function_size
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
array<double, BasisTraits<Linear, 1>::function_size>
basis_shift<Linear, 1>(Point<1> x, Point<1> y) {
  array<double, BasisTraits<Linear, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = y[0]-x[0];
  return r;
}

template<>
array<double, BasisTraits<Linear, 2>::function_size>
basis_shift<Linear, 2>(Point<2> x, Point<2> y) {
  array<double, BasisTraits<Linear, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = y[0]-x[0];
  r[2] = y[1]-x[1];
  return r;
}

template<>
array<double, BasisTraits<Linear, 3>::function_size>
basis_shift<Linear, 3>(Point<3> x, Point<3> y) {
  array<double, BasisTraits<Linear, 3>::function_size> r{0.0};
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
array<double, BasisTraits<Quadratic, 1>::function_size>
basis_function<Quadratic, 1>(Point<1> x) {
  array<double, BasisTraits<Quadratic, 1>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = .5 * x[0] * x[0];
  return r;
}

template<>
array<double, BasisTraits<Quadratic, 2>::function_size>
basis_function<Quadratic, 2>(Point<2> x) {
  array<double, BasisTraits<Quadratic, 2>::function_size> r{0.0};
  r[0] = 1.;
  r[1] = x[0];
  r[2] = x[1];
  r[3] = .5 * x[0] * x[0];
  r[4] = x[0] * x[1];
  r[5] = .5 * x[1] * x[1];
  return r;
}

template<>
array<double, BasisTraits<Quadratic, 3>::function_size>
basis_function<Quadratic, 3>(Point<3> x) {
  array<double, BasisTraits<Quadratic, 3>::function_size> r{0.0};
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
array<double, BasisTraits<Quadratic, 1>::function_size>,
BasisTraits<Quadratic, 1>::function_size
>
basis_jet<Quadratic,1>(Point<1> x){
  array<
  array<double, BasisTraits<Quadratic, 1>::function_size>,
  BasisTraits<Quadratic, 1>::function_size
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
array<double, BasisTraits<Quadratic, 2>::function_size>,
BasisTraits<Quadratic, 2>::function_size
>
basis_jet<Quadratic,2>(Point<2> x){
  array<
  array<double, BasisTraits<Quadratic, 2>::function_size>,
  BasisTraits<Quadratic, 2>::function_size
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
array<double, BasisTraits<Quadratic, 3>::function_size>,
BasisTraits<Quadratic, 3>::function_size
>
basis_jet<Quadratic,3>(Point<3> x){
  array<
  array<double, BasisTraits<Quadratic, 3>::function_size>,
  BasisTraits<Quadratic, 3>::function_size
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
array<double, BasisTraits<Quadratic, 1>::function_size>,
BasisTraits<Quadratic, 1>::function_size
>
basis_inverse_jet<Quadratic,1>(Point<1> x){
  Point<1> mx;
  for (size_t i=0; i<1; i++) mx[i] = -x[i];
  return basis_jet<Quadratic,1>(mx);
}

template<>
array<
array<double, BasisTraits<Quadratic, 2>::function_size>,
BasisTraits<Quadratic, 2>::function_size
>
basis_inverse_jet<Quadratic,2>(Point<2> x){
  Point<2> mx;
  for (size_t i=0; i<2; i++) mx[i] = -x[i];
  return basis_jet<Quadratic,2>(mx);
}

template<>
array<
array<double, BasisTraits<Quadratic, 3>::function_size>,
BasisTraits<Quadratic, 3>::function_size
>
basis_inverse_jet<Quadratic,3>(Point<3> x){
  Point<3> mx;
  for (size_t i=0; i<3; i++) mx[i] = -x[i];
  return basis_jet<Quadratic,3>(mx);
}

//---------------------------------------------------------------

template<>
array<double, BasisTraits<Quadratic, 1>::function_size>
basis_shift<Quadratic, 1>(Point<1> x, Point<1> y) {
  Point<1> d;
  for (size_t i=0; i<1; i++) d[i] = y[i]-x[i];
  return basis_function<Quadratic,1>(d);
}

template<>
array<double, BasisTraits<Quadratic, 2>::function_size>
basis_shift<Quadratic, 2>(Point<2> x, Point<2> y) {
  Point<2> d;
  for (size_t i=0; i<2; i++) d[i] = y[i]-x[i];
  return basis_function<Quadratic,2>(d);
}

template<>
array<double, BasisTraits<Quadratic, 3>::function_size>
basis_shift<Quadratic, 3>(Point<3> x, Point<3> y) {
  Point<3> d;
  for (size_t i=0; i<3; i++) d[i] = y[i]-x[i];
  return basis_function<Quadratic,3>(d);
}

}
}
}

#endif // BASIS_H_INC
