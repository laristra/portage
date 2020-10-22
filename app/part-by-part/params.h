/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include <fstream>
#include <memory>
#include <map>
#include "mpi.h" // mandatory

// parsers
#include "json.h"
#include "user_field.h"
#include "filter.h"

// jali
#include "MeshFactory.hh"
#include "LabeledSetRegion.hh"
#include "GeometricModel.hh"

// meshes and states
#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// remap kernels and driver
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/driver/coredriver.h"
#include "portage/support/mpi_collate.h"
#include "portage/support/timer.h"

/**
 * @brief Encaspulate parameters, parsing methods and input checks.
 *
 */
class Params {

  /**
   * @brief Structure to store part pair meta-data.
   *
   */
  struct Part {
    /** source part: whether expression or block indices */
    struct { std::string expr; std::set<int> blocks; } source;
    /** target part: whether expression or block indices */
    struct { std::string expr; std::set<int> blocks; } target;
  };

public:
  /** MPI rank */
  int rank = 0;
  /** number of MPI ranks */
  int nb_ranks = 1;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_WORLD;
  /** mesh dimension */
  int dimension = 2;
  /** is mesh conformal */
  bool conformal = true;
  /** export mesh and states */
  bool dump = false;
  /** is source field internal */
  bool internal = false;
  /** source mesh file */
  std::string source;
  /** target mesh file */
  std::string target;
  /** remap order */
  int order = 1;
  /** numerical tolerance for remap */
  double tolerance = Portage::DEFAULT_NUMERIC_TOLERANCES<2>.relative_conservation_eps;
  /** upper bound for remap */
  double const upper_bound = std::numeric_limits<double>::max();
  /** lower bound for remap */
  double const lower_bound = -upper_bound;
  /** field kind: nodal or cell-centered */
  Wonton::Entity_kind kind = Wonton::CELL;
  /** name and expression of each source field */
  std::map<std::string, std::string> fields {};
  /** name and expression of each part pair */
  std::map<std::string, std::vector<Part>> parts {};
  /** number of fixup iterations */
  int fix_iter = 5;
  /** gradient limiter to use for internal entities */
  Portage::Limiter_type limiter = Portage::NOLIMITER;
  /** gradient limiter to use for boundary entities */
  Portage::Boundary_Limiter_type bnd_limiter = Portage::BND_NOLIMITER;
  /** empty cells fixup policy */
  Portage::Empty_fixup_type empty_fixup = Portage::EXTRAPOLATE;
  /** partially overlapped cells fixup policy */
  Portage::Partial_fixup_type partial_fixup = Portage::GLOBALLY_CONSERVATIVE;

  /**
   * @brief Parse and store parameters.
   *
   * @param argc: number of app arguments
   * @param argv: list of app arguments
   * @return parsing status flag.
   */
  int parse(int argc, char* argv[]) {
    if (argc != 2 or argv == nullptr)
      return abort("wrong arguments");

    if (std::string(argv[1]) == "-h")
      return abort();

    std::ifstream file(argv[1]);
    if (not file.good())
      return abort("unable to open input file", false);

    nlohmann::json json;

    try {

      file >> json;
      file.close();

      // check if every parameter is specified
      if (not json.count("mesh"))
        return abort("missing mesh attributes");
      else {
        if (not json["mesh"].count("dimension"))
          return abort("unspecified mesh dimension");

        if (not json["mesh"].count("source"))
          return abort("unspecified source mesh");

        if (not json["mesh"].count("target"))
          return abort("unspecified target mesh");

        if (not json["mesh"].count("conformal"))
          return abort("must precise if conformal mesh");

        if (not json["mesh"].count("export"))
          return abort("must precise if export results or not");
      }

      if(not json.count("remap"))
        return abort("missing remap attributes");
      else {
        if (not json["remap"].count("kind"))
          return abort("unspecified entity kind for remap");

        if (not json["remap"].count("order"))
          return abort("unspecified order of accuracy for remap");

        if (not json["remap"].count("limiter"))
          return abort("unspecified default gradient limiter");

        if (not json["remap"].count("bnd_limiter"))
          return abort("unspecified default gradient boundary limiter");

        if (not json["remap"].count("fixup"))
          return abort("unspecified mismatch fixup parameters");
        else {
          if (not json["remap"]["fixup"].count("partial"))
            return abort("unspecified partially filled cells fixup scheme");
          if (not json["remap"]["fixup"].count("empty"))
            return abort("unspecified empty cells fixup scheme");
          if (not json["remap"]["fixup"].count("max-iter"))
            return abort("unspecified maximum number of fixup iterations");
        }

        if (not json["remap"].count("fields"))
          return abort("unspecified material fields");
        else {
          for (auto&& field : json["remap"]["fields"]) {
            if (not field.count("name")) { return abort("unknown field name"); }
            if (not field.count("internal") or not field["internal"]) {
              if (not field.count("expr")) { return abort("no field expression"); }
            }
          }
        }
      }

      if (not json.count("parts"))
        return abort("no parts field in parameter file");
      else {

        for (auto&& entry : json["parts"]) {
          if (not entry.count("field"))
            return abort("no field name for part");

          if (entry.count("expr")) {
            for (auto&& pair : entry["expr"]) {
              if (not pair.count("source"))
                return abort("no source part");
              if (not pair.count("target"))
                return abort("no target part");
            }
          } else if (entry.count("blocks")) {
            for (auto&& pair : entry["blocks"]) {
              if (not pair.count("source"))
                return abort("no source part");
              if (not pair.count("target"))
                return abort("no target part");
            }
          } else
            return abort("undefined parts");
        }
      }

      // then store them
      dimension = json["mesh"]["dimension"];
      conformal = json["mesh"]["conformal"];
      dump      = json["mesh"]["export"];
      source    = json["mesh"]["source"];
      target    = json["mesh"]["target"];
      order     = json["remap"]["order"];
      fix_iter  = json["remap"]["fixup"]["max-iter"];

      bool const use_limiter  = json["remap"]["limiter"];
      bool const use_bnd_limiter  = json["remap"]["bnd_limiter"];
      std::string remap_kind  = json["remap"]["kind"];
      std::string parts_fixup = json["remap"]["fixup"]["partial"];
      std::string empts_fixup = json["remap"]["fixup"]["empty"];

      if (remap_kind == "cell")
        kind = Wonton::Entity_kind::CELL;
      else
        kind = Wonton::Entity_kind::NODE;

      if (use_limiter)
        limiter = Portage::Limiter_type::BARTH_JESPERSEN;
      else
        limiter = Portage::Limiter_type::NOLIMITER;

      if (use_bnd_limiter)
        bnd_limiter = Portage::Boundary_Limiter_type::BND_ZERO_GRADIENT;
      else
        bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;

      if (parts_fixup == "locally_conservative")
        partial_fixup = Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE;
      else if (parts_fixup == "constant")
        partial_fixup = Portage::Partial_fixup_type::CONSTANT;
      else
        partial_fixup = Portage::Partial_fixup_type::GLOBALLY_CONSERVATIVE;

      if (empts_fixup == "leave_empty")
        empty_fixup = Portage::Empty_fixup_type::LEAVE_EMPTY;
      else
        empty_fixup = Portage::Empty_fixup_type::EXTRAPOLATE;

      /* parts field */
      for (auto&& scalar : json["remap"]["fields"]) {
        bool is_internal = scalar.count("internal") and scalar["internal"];
        fields[scalar["name"]] = is_internal ? "internal" : scalar["expr"];
      }

      for (auto&& entry : json["parts"]) {
        std::string field = entry["field"];
        Part part;
        if (entry.count("expr")) {
          for (auto&& expr : entry["expr"]) {
            part.source.expr = expr["source"];
            part.target.expr = expr["target"];
            parts[field].emplace_back(part);
          }
        } else {
          for (auto&& block : entry["blocks"]) {
            part.source.blocks.clear();
            part.target.blocks.clear();
            for (int c : block["source"]) { part.source.blocks.insert(c); }
            for (int c : block["target"]) { part.target.blocks.insert(c); }
            parts[field].emplace_back(part);
          }
        }
        assert(not parts[field].empty());
      }

      // check their validity eventually
      file.open(source);
      if (not file.good())
        return abort("unable to read source mesh file", false);
      file.close();

      file.open(target);
      if (not file.good())
        return abort("unable to read target mesh file", false);

      if (dimension < 2 or dimension > 3)
        return abort("invalid mesh dimension [2|3]", false);

      if (kind == Wonton::NODE)
        return abort("multi-part node remap is not supported", false);

      if (order < 1 or order > 2)
        return abort("invalid order of accuracy for remap [1|2]", false);

      // check that fields and parts have same keys
      auto same_keys = [](auto const& a, auto const& b) { return a.first == b.first; };

      bool have_same_size = fields.size() == parts.size();
      bool have_same_keys = std::equal(fields.begin(), fields.end(),
                                       parts.begin(), same_keys);

      if (not have_same_size or not have_same_keys)
        return abort("numerical fields and per-part fields mismatch");


    } catch(nlohmann::json::parse_error& e) {
      return abort(e.what());
    }

    // everything was ok
    if (rank == 0)
      print(json);

    return EXIT_SUCCESS;
  }

  /**
   * @brief Handle runtime errors.
   *
   * @param message: a message to be displayed if any.
   * @param show_usage: hint for usage instructions print.
   * @return status
   */
  int abort(std::string const& message = "", bool show_usage = true) const {
    if (rank == 0) {
      if (show_usage) { print_usage(); }

      if (not message.empty()) {
        std::fprintf(stderr,
          " ---------------------------------------------------------- \n"
          " \e[31m Error: %s. \e[0m                                    \n"
          " ---------------------------------------------------------- \n",
         message.data()
        );
      }
    }

    MPI_Finalize();
    return EXIT_FAILURE;
  }

  /**
   * @brief Get number of digits of the given scalar.
   *
   * @tparam: scalar type
   * @param number: scalar value
   * @return its number of digits for print.
   */
  template<typename T>
  int get_number_digit(T number) const {
    return (number > 0 ? static_cast<int>(std::floor(std::log10(number))) + 1 : 0);
  }

private:
  /**
   * @brief Display run instructions and input file format.
   *
   */
  static void print_usage() {
    std::printf(
      " Usage: mpirun -np [nranks] ./part-remap [file.json]                \n\n"
      " [file.json]:                                                         \n"
      " {                                                                    \n"
      "  \e[32m'mesh'\e[0m: {                                                \n"
      "    \e[32m'dimension'\e[0m: <2|3>,                                    \n"
      "    \e[32m'source'\e[0m: '/path/to/source/mesh.exo',                  \n"
      "    \e[32m'target'\e[0m: '/path/to/target/mesh.exo',                  \n"
      "    \e[32m'conformal'\e[0m: <boolean>                                 \n"
      "    \e[32m'export'\e[0m: <boolean>                                    \n"
      "  },                                                                  \n"
      "  \e[32m'remap'\e[0m: {                                               \n"
      "    \e[32m'kind'\e[0m: 'cell',                                        \n"
      "    \e[32m'order'\e[0m: <1|2>,                                        \n"
      "    \e[32m'limiter'\e[0m: <boolean>,                                  \n"
      "    \e[32m'bnd_limiter'\e[0m: <boolean>,                              \n"
      "    \e[32m'fixup'\e[0m: {                                             \n"
      "      \e[32m'partial'\e[0m: 'constant|<locally|shifted>_conservative',\n"
      "      \e[32m'empty'\e[0m: '<leave_empty|extrapolate>',                \n"
      "      \e[32m'max-iter'\e[0m: <unsigned integer>,                      \n"
      "    },                                                                \n"
      "    \e[32m'fields'\e[0m: [                                            \n"
      "      { \e[32m'name'\e[0m:'density',\e[32m'expr'\e[0m: '<math>' }     \n"
      "      { \e[32m'name'\e[0m:'temperature', 'internal': true }           \n"
      "    ]                                                                 \n"
      "  },                                                                  \n"
      "  \e[32m'parts'\e[0m: [                                               \n"
      "    {                                                                 \n"
      "      \e[32m'field'\e[0m: 'density',                                  \n"
      "      \e[32m'expr'\e[0m: [                                            \n"
      "        { \e[32m'source'\e[0m: <math>, \e[32m'target'\e[0m: <math> }, \n"
      "        { \e[32m'source'\e[0m: <math>, \e[32m'target'\e[0m: <math> }  \n"
      "      ]                                                               \n"
      "    },                                                                \n"
      "    {                                                                 \n"
      "      \e[32m'field'\e[0m: 'temperature',                              \n"
      "      \e[32m'blocks'\e[0m: [                                          \n"
      "        { \e[32m'source'\e[0m: [1,2,3], \e[32m'target'\e[0m: [5] },   \n"
      "        { \e[32m'source'\e[0m: [2,4], \e[32m'target'\e[0m: [7,8] }    \n"
      "      ]                                                               \n"
      "    }                                                                 \n"
      "  ]                                                                   \n"
      " }                                                                    \n"
    );
  }

  /**
   * @brief Display input parameters.
   *
   * @param json: the JSON structure that contains the input.
   */
  void print(nlohmann::json const& json) const {

    // extract file path base name
    auto base = [](const std::string& path) -> std::string {
      std::string b = path;
      size_t i = b.find_last_not_of('/');
      if (i == std::string::npos) {
        if (b[0] == '/') b.erase(1);
        return b;
      }

      b.erase(i + 1, b.length() - i - 1);
      i = b.find_last_of('/');
      if (i != std::string::npos)
        b.erase(0, i + 1);
      return b;
    };

    auto display = [](auto const& list) -> std::string {
      int const size = list.size();
      std::stringstream buffer;
      int i = 0;
      for (auto const& id : list) {
        buffer << id;
        if (++i < size)
          buffer << ",";
      }
      return buffer.str();
    };

    std::string remap_kind  = json["remap"]["kind"];
    std::string dump_result = json["mesh"]["export"] ? "yes": "no";
    std::string use_limiter = json["remap"]["limiter"] ? "yes": "no";
    std::string use_bnd_limiter = json["remap"]["bnd_limiter"] ? "yes": "no";
    std::string partial_fix = json["remap"]["fixup"]["partial"];
    std::string empty_fix   = json["remap"]["fixup"]["empty"];

    std::printf("Parameters: \n");
    std::printf(" \u2022 MPI ranks: \e[32m%d\e[0m\n", nb_ranks);
    std::printf(" \u2022 dimension: \e[32m%d\e[0m\n", dimension);
    std::printf(" \u2022 source mesh: '\e[32m%s\e[0m'\n", base(source).data());
    std::printf(" \u2022 target mesh: '\e[32m%s\e[0m'\n", base(target).data());
    std::printf(" \u2022 dump results: \e[32m%s\e[0m\n", dump_result.data());
    std::printf(" \u2022 remap order: \e[32m%d\e[0m\n", order);
    std::printf(" \u2022 remap kind: \e[32m%s\e[0m\n", remap_kind.data());
    std::printf(" \u2022 use limiter: \e[32m%s\e[0m\n", use_limiter.data());
    std::printf(" \u2022 use boundary limiter: \e[32m%s\e[0m\n", use_bnd_limiter.data());
    std::printf(" \u2022 partial filled fixup: \e[32m%s\e[0m\n", partial_fix.data());
    std::printf(" \u2022 empty cells fixup: \e[32m%s\e[0m\n", empty_fix.data());
    std::printf(" \u2022 fixup iterations: \e[32m%d\e[0m\n", fix_iter);
    std::printf(" \u2022 numerical fields: \n");

    for (auto&& field: fields) {
      std::printf("   - %s: '\e[32m%s\e[0m'\n", field.first.data(), field.second.data());

      int id = 1;
      for (auto&& part : parts.at(field.first)) {
        std::printf("   \e[32m%d\e[0m.", id++);
        if (not part.source.expr.empty() and not part.target.expr.empty()) {
          std::printf(" source: '\e[32m%s\e[0m',", part.source.expr.data());
          std::printf(" target: '\e[32m%s\e[0m'\n", part.target.expr.data());
        } else {
          assert(not part.source.blocks.empty() and not part.target.blocks.empty());
          std::printf(" source: '\e[32m%s\e[0m', target: '\e[32m%s\e[0m'\n",
                      display(part.source.blocks).data(),
                      display(part.target.blocks).data());
        }
      }
      std::printf("\n");
    }
  }
};
