/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include "exprtk.hpp"
#include "portage/support/portage.h"

/**
 * @brief Assess analytical formula.
 *
 */
class Formula {
public:

  /**
   * @brief Initialize the parser.
   *
   * @param formula: the analytical expression to parse.
   * @return true if everything correctly set, false otherwise.
   */
  bool initialize(std::string const& formula) {
    // register variables and symbol table
    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    expression.register_symbol_table(symbol_table);

    // parse expression
    if (not parser.compile(formula, expression)) {
      std::fprintf(stderr,
        "Could not parse the point filter expression '%s'\n"
        "Please check http://www.partow.net/programming/exprtk "
        "for more details on syntax description\n",
        formula.data()
      );
      return false;
    }
    return true;
  }

  /**
   * @brief Evalute the function on the given point.
   *
   * @param p: the point.
   * @return the value of the function at p.
   */
  double operator()(Wonton::Point<2> const& p) {
    x = p[0];
    y = p[1];
    return expression.value();
  }

private:
  double x = 0.0;                             /** coordinate in x-axis */
  double y = 0.0;                             /** coordinate in y-axis */
  exprtk::symbol_table<double> symbol_table;  /** table of variables */
  exprtk::expression<double> expression;      /** the expression to parse */
  exprtk::parser<double> parser;              /** the actual parser */
};
