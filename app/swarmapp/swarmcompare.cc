/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iostream>

namespace {
// numerical tolerance
const double epsilon = 1.e-12;

/**
 * @brief Parsed data.
 */
class Data {
public:
  void add(std::string const& x, std::string const& y, std::string const& v) {
    coords.emplace_back(std::vector<double>{std::stod(x), std::stod(y)});
    values.emplace_back(std::stod(v));
    if (not dim)
      dim = 2;
  }

  void add(std::string const& x, std::string const& y,
           std::string const& z, std::string const& v) {
    coords.emplace_back(std::vector<double>{std::stod(x), std::stod(y), std::stod(z)});
    values.emplace_back(std::stod(v));
    if (not dim)
      dim = 3;
  }

  int size() const { return values.size(); }

  std::vector<std::vector<double>> coords {};
  std::vector<double> values {};
  int dim = 0;
};

/**
 * @brief Parse a single line of a csv file.
 *
 * @param file: input file
 * @return a list of field values
 */
std::vector<std::string> parse_line(std::ifstream& file) {

  std::vector<std::string> result;
  std::string line;
  std::string field;
  std::stringstream buffer;

  std::getline(file, line);
  if (line.empty())
    return result;

  //std::stringstream buffer(line, std::ios_base::in);
  buffer << line;

  while (std::getline(buffer, field, ',')) {
    result.emplace_back(field);
  }
  // check for a trailing comma with no data after it.
  if (not buffer and field.empty())
    result.emplace_back("");
  return result;
}

/**
 * @brief Parse the file and retrieve field values.
 *
 * @param file
 * @return
 */
Data parse(const char* path) {

  assert(path != nullptr);
  std::ifstream file(path);
  assert(file.good());

  Data data;
  bool stop = false;
  parse_line(file); // skip file header

  while (not stop) {
    auto cols = parse_line(file);
    int nb_col = cols.size();

    switch (nb_col) {
      case 0: stop = true; break;
      case 3: data.add(cols[0], cols[1], cols[2]); break;
      case 4: data.add(cols[0], cols[1], cols[2], cols[3]); break;
      default: throw std::runtime_error("uneven columns in file");
    }
  }
  return data;
}

} // end anonymous namespace

/**
 * @brief Main method.
 *
 * @param argc: number of arguments.
 * @param argv: array of arguments.
 * @return status code, 0 if ok.
 */
int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr << "Error: not enough arguments" << std::endl;
    return EXIT_FAILURE;
  }

  auto const gold = parse(argv[1]);
  auto const test = parse(argv[2]);

  if (gold.size() != test.size()) {
    std::cerr << "Error: mismatch on gold and test data sizes" << std::endl;
    return EXIT_FAILURE;
  }

  // check if two real values are equal
  auto equals = [](double const& u, double const& v) -> bool {
    return std::abs(u - v) < epsilon;
  };

  int const dim = gold.dim;

  for (int i = 0; i < gold.size(); i++) {
    bool coords_match = true;
    for (int j = 0; j < dim; ++j)
      coords_match &= equals(gold.coords[i][j], test.coords[i][j]);

    if (not coords_match) {
      std::cerr << "Error: mismatched coordinates at line " << i << std::endl;
      std::cerr << "\tgold: [";
      for (int j = 0; j < dim; ++j) {
        std::cerr << gold.coords[i][j];
        if (j < dim - 1)
          std::cerr << ", ";
      }
      std::cerr << "], test: [";
      for (int j = 0; j < dim; ++j) {
        std::cerr << test.coords[i][j];
        if (j < dim - 1)
          std::cerr << ", ";
      }
      std::cerr << "]." << std::endl;
      return EXIT_FAILURE;
    }

    if (not equals(gold.values[i], test.values[i])) {
      std::cerr << "Error: mismatched values at line " << i << std::endl;
      std::cerr << "\tgold: " << gold.values[i] <<", test: " << test.values[i] << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "files agree" << std::endl;
  return EXIT_SUCCESS;
}
