/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGEAPP_USER_FIELD_H_
#define PORTAGEAPP_USER_FIELD_H_

// ExprTk package for expression parsing
// See http://www.partow.net/programming/exprtk/ for source and examples

#include "exprtk.hpp"
#include "wonton/support/Point.h"
#ifdef HAVE_JALI
  #include "Point.hh"
#endif

// This functor initializes a general field from a string expression
// and returns its value at any given point

class user_field_t {
public:
  double x = 0.;
  double y = 0.;
  double z = 0.;
  int dim_ = 0;

  exprtk::symbol_table<double> symbol_table {};
  exprtk::expression<double> expression {};
  exprtk::parser<double> parser {};

  user_field_t() = default;
  ~user_field_t() = default;

  int initialize(int dim, std::string const& expression_str) {
    dim_ = dim;
    assert(dim_ >= 1 && dim_ <= 3);
    if (dim_ > 0) symbol_table.add_variable("x", x);
    if (dim_ > 1) symbol_table.add_variable("y", y);
    if (dim_ > 2) symbol_table.add_variable("z", z);

    x = y = z = 0.0;

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(expression_str, expression)) {
      std::cerr << "Could not parse field expression " << expression_str <<
          std::endl;
      std::cerr << "One culprit could be that the variables in the expression " 
                << "do not match the specified dimension of the problem" 
                << std::endl;
      std::cerr << "Also check the page http://www.partow.net/programming/exprtk/ for syntax description" << std::endl;
      return 0;
    }
    
    return 1;
  }

  template<int D>
  double operator()(Wonton::Point<D> const& c) {
    if (dim_ == D) {
      if (dim_ > 0) x = c[0];
      if (dim_ > 1) y = c[1];
      if (dim_ > 2) z = c[2];
      double const val = expression.value();
      return val;
    } else
      throw std::runtime_error("incompatible dimensions");
  }

#ifdef HAVE_JALI
  double operator()(JaliGeometry::Point const& c) {
    if (dim_ == c.dim()) {
      if (dim_ > 0) x = c[0];
      if (dim_ > 1) y = c[1];
      if (dim_ > 2) z = c[2];
      double const val = expression.value();
      return val;
    } else
      throw std::runtime_error("incompatible dimensions");
  }
#endif
};



#endif  // PORTAGEAPP_USER_FIELD_H_
