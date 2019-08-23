  /*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

// ExprTk package for expression parsing
// See http://www.partow.net/programming/exprtk/ for source and examples

#include "exprtk.hpp"

// This functor initializes a general field from a string expression
// and returns its value at any given point

class filter_t {

public:

   filter_t() = default;
  ~filter_t() = default;

  bool initialize(int dim, std::string const& expression_str) {
    assert(dim >= 1 and dim <= 3);

    // register variables and symbol table
    if (dim > 0) symbol_table.add_variable("x", x);
    if (dim > 1) symbol_table.add_variable("y", y);
    if (dim > 2) symbol_table.add_variable("z", z);
    expression.register_symbol_table(symbol_table);

    // parse expression
    if (!parser.compile(expression_str, expression)) {
      std::fprintf(stderr,
        "Could not parse the point filter expression '%s'\n"
        "One culprit could be that one or more variables"
        "do not match the specified dimension of the problem, "
        "or the expression itself does not return a boolean.\n"
        "Please check http://www.partow.net/programming/exprtk "
        "for more details on syntax description\n",
          expression_str.data()
      );
      return false;
    }
    dim_ = dim;
    return true;
  }

  template<class Point>
  bool operator()(Point const& c) {
    if (dim_ > 0) x = c[0];
    if (dim_ > 1) y = c[1];
    if (dim_ > 2) z = c[2];
    return (bool) expression.value();
  }

private:
  double x = 0;
  double y = 0;
  double z = 0;
  int dim_ = 0;

  exprtk::symbol_table<double> symbol_table;
  exprtk::expression<double> expression;
  exprtk::parser<double> parser;
};