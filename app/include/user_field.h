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

// This functor initializes a general field from a string expression
// and returns its value at any given point

typedef struct user_field {
  double x, y, z;
  int dim_;
  exprtk::symbol_table<double> symbol_table;
  exprtk::expression<double> expression;
  exprtk::parser<double> parser;

  user_field() : dim_(0) {}

  user_field(int dim, std::string const& expression_str) : dim_(dim) {
    initialize(dim, expression_str);
  }

  void initialize(int dim, std::string const& expression_str) {
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
      exit(-1);  // Exit the program gracefully
    }
  }

  template<class P> double operator()(const P &c) {
    if (dim_ > 0) x = c[0];
    if (dim_ > 1) y = c[1];
    if (dim_ > 2) z = c[2];
    double val = expression.value();
    return val;
  }
} user_field_t;



#endif  // PORTAGEAPP_USER_FIELD_H_
