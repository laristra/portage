/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <fstream>
#include <memory>
#ifdef PORTAGE_ENABLE_MPI
  #include "mpi.h"
#endif

// jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_mm_wrapper.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/driver/coredriver.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"

// parse input
#include "json.h"

namespace app {

  // get rid of long namespaces within the app when no possible confusion.
  constexpr auto CELL = Wonton::Entity_kind::CELL;
  constexpr auto NODE = Wonton::Entity_kind::NODE;
  constexpr auto ALL  = Wonton::Entity_type::ALL;
  constexpr auto NO_LIMITER = Portage::Limiter_type::NOLIMITER;
  constexpr auto BARTH_JESPERSEN = Portage::Limiter_type::BARTH_JESPERSEN;
  constexpr unsigned const default_heap_size = 512;

  // a compact object to store app parameters.
  struct Params {
    int dimension;                                      // meshes dimension
    bool is_conform;                                    // conformal meshes?
    bool do_export;                                     // dump results?
    std::string source_file;                            // source mesh file
    std::string target_file;                            // target mesh file

    int order;                                          // remap accuracy order
    Portage::Limiter_type limiter;                      // gradient limiter
    Wonton::Entity_kind entity_kind;                    // node|cell-based remap

    struct Entities {                                   // parts pair data
      unsigned id;                                      // index of part pair
      std::vector<int> source;                          // related source entities
      std::vector<int> target;                          // related target entities
    };
    std::map<std::string, std::string> fields;          // scalar fields
    std::map<std::string, std::vector<Entities>> parts; // meshes parts to remap
  };
} // namespace 'app'

/**
 * @brief Display run instructions and input file format.
 *
 */
void print_usage() {

  std::printf(" Usage: ./part-remap [params.json]                          \n");
  std::printf(" [params.json]:                                             \n");
  std::printf(" {                                                          \n");
  std::printf("  'mesh': {                                                 \n");
  std::printf("    'dimension': <2|3>,                                     \n");
  std::printf("    'source': '/path/to/source/mesh.exo',                   \n");
  std::printf("    'target': '/path/to/target/mesh.exo',                   \n");
  std::printf("    'conformal': <boolean>                                  \n");
  std::printf("    'export': <boolean>                                     \n");
  std::printf("  },                                                        \n");
  std::printf("  'remap': {                                                \n");
  std::printf("    'kind': 'cell',                                         \n");
  std::printf("    'order': <1|2>,                                         \n");
  std::printf("    'limiter': <0-2>,                                       \n");
  std::printf("    'fields': [                                             \n");
  std::printf("      { 'name': 'density', 'expr': '<math>' }               \n");
  std::printf("      { 'name': 'temperature', 'expr': '<math>' }           \n");
  std::printf("    ]                                                       \n");
  std::printf("  },                                                        \n");
  std::printf("  'parts': [                                                \n");
  std::printf("    {                                                       \n");
  std::printf("      'field': 'density',                                   \n");
  std::printf("      'entities': [                                         \n");
  std::printf("        { 'uid': 1, 'source': [0,2,5], 'target': [0,3]   }, \n");
  std::printf("        { 'uid': 2, 'source': [1,3,4], 'target': [1,2,4] }  \n");
  std::printf("      ]                                                     \n");
  std::printf("    }                                                       \n");
  std::printf("  ]                                                         \n");
  std::printf(" }                                                          \n");
  std::printf(" ---------------------------------------------------------- \n");
}

static app::Params params;

/**
 * @brief Handle errors during argument parsing.
 *
 * @param 'message'    a message to be displayed if any.
 * @param 'show_usage' hint for usage instructions print.
 * @return false
 */
bool abort(std::string message, bool show_usage = true) {
  if (show_usage)
    print_usage();

  std::fprintf(stderr, "** Error: %s. **\n", message.data());
  return false;
}

/**
 * @brief Retrieve and store parameters.
 *
 * @param 'path' the JSON parameter file path.
 * @return true if no parsing errors, false otherwise.
 */
bool parse_params(std::string path) {

  std::ifstream file(path);
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

      if (not json["remap"].count("fields"))
        return abort("unspecified material fields");
      else {
        for (auto&& field : json["remap"]["fields"]) {
          if (not field.count("name")) { return abort("unknown field name"); }
          if (not field.count("expr")) { return abort("no field expression"); }
        }
      }
    }

    if (not json.count("parts"))
      return abort("no parts field in parameter file");
    else {
      // lookup table to check uid
      std::set<int> helper;

      for (auto&& entry : json["parts"]) {
        if (not entry.count("field"))
          return abort("no field name for part");
        if (not entry.count("pairs"))
          return abort("no given entities for part");
        else {
          for (auto&& pair : entry["pairs"]) {
            // check uid
            if (not pair.count("uid"))
              return abort("no given unique id for part pair");
            else {
              int const uid = pair["uid"];
              if (helper.count(uid))
                return abort("already used uid for part pair");
              else
                helper.insert(uid);
            }
            // check source list
            if (not pair.count("source"))
              return abort("no source data for part pair");
            else if (pair["source"].empty())
              return abort("empty source entities list for part pair");
            // check target list
            if (not pair.count("target"))
              return abort("no target entities for part pair");
            else if (pair["target"].empty())
              return abort("empty target entities list for part pair");
          }
        }
      }
      helper.clear();
    }

    // then store them
    params.dimension   = json["mesh"]["dimension"];
    params.is_conform  = json["mesh"]["conformal"];
    params.do_export   = json["mesh"]["export"];
    params.source_file = json["mesh"]["source"];
    params.target_file = json["mesh"]["target"];
    params.order       = json["remap"]["order"];
    params.entity_kind = json["remap"]["kind"] == "cell" ? app::CELL : app::NODE;
    params.limiter     = json["remap"]["limiter"] ? app::NO_LIMITER
                                                  : app::BARTH_JESPERSEN;

    for (auto&& scalar : json["remap"]["fields"])
      params.fields[scalar["name"]] = scalar["expr"];

    for (auto&& entry : json["parts"]) {
      auto const field = entry["field"];
      params.parts[field].reserve(app::default_heap_size);

      for (auto&& pair : entry["entities"]) {
        app::Params::Entities parts;
        parts.id = pair["uid"];
        for (auto&& s : pair["source"]) { parts.source.push_back(s); }
        for (auto&& t : pair["target"]) { parts.target.push_back(t); }
        params.parts[field].emplace_back(parts);
      }
      params.parts[field].shrink_to_fit();
    }

    // check their validity eventually
    file.open(params.source_file);
    if (not file.good())
      return abort("unable to read source mesh file", false);
    else
      file.close();

    file.open(params.target_file);
    if (not file.good())
      return abort("unable to read target mesh file", false);

    if (params.dimension < 2 or params.dimension > 3)
      return abort("invalid mesh dimension [2|3]", false);

    if (params.order != 1 and params.order != 2)
      return abort("invalid order of accuracy for remap [0-1]", false);

  } catch(nlohmann::json::parse_error& e) {
    return abort(e.what());
  }
  // everything was ok
  return true;
}


int main(int argc, char* argv[]) {

  std::printf(" ---------------------------------------------------------- \n");
  std::printf("  Demonstration app for multi-part field interpolation.     \n");
  std::printf("  It handles pure cells remap of multi-material meshes and  \n");
  std::printf("  non-smoothed remap of fields with sharp discontinuities.  \n");
  std::printf(" ---------------------------------------------------------- \n");

  if (parse_params(argv[1])) {
    std::printf("Everything was fine\n");
  }
  return EXIT_SUCCESS;
}