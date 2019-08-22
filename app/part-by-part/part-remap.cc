/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <fstream>
#include <memory>
#include <chrono>
#include <map>
#include "mpi.h" // cannot work without MPI

// meshes and states
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "MeshFactory.hh"

// remap kernels and driver
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/driver/coredriver.h"

// parsers
#include "json.h"
#include "user_field.h"

namespace app {
  namespace entity {
    auto const cell = Wonton::Entity_kind::CELL;
    auto const node = Wonton::Entity_kind::NODE;
    // entities parts pair data
    struct part {
      int id;
      std::string source;
      std::string target;
    };
  }
  namespace limiter {
    auto const none  = Portage::Limiter_type::NOLIMITER;
    auto const barth = Portage::Limiter_type::BARTH_JESPERSEN;
  }

  // a compact object to store app parameters.
  struct Params {
    // mesh
    int dimension;                                          // meshes dimension
    bool conformal;                                         // conformal meshes?
    bool dump;                                              // export results?
    std::string source;                                     // source mesh file
    std::string target;                                     // target mesh file

    // remap
    int order;                                              // accuracy order
    Portage::Limiter_type limiter;                          // gradient limiter
    Wonton::Entity_kind   kind;                             // node|cell-based
    std::map<std::string, std::string> fields;              // fields expression
    std::map<std::string, std::vector<entity::part>> parts; // per-field parts
  };
}

namespace timer {
  // get rid of long namespaces
  using time_t = std::chrono::high_resolution_clock::time_point;

  // get current time point
  inline time_t now() { return std::chrono::high_resolution_clock::now(); }

  // retrieve elapsed time in seconds.
  inline float elapsed(time_t& tic, bool reset = false) {
    auto const secs = static_cast<float>(
      std::chrono::duration_cast<std::chrono::seconds>(now() - tic).count()
    );

    if (reset) tic = now();
    return secs;
  }
}


/**
 * @brief Display run instructions and input file format.
 *
 */
void print_usage() {

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
    "    \e[32m'limiter'\e[0m: <boolean>                                   \n"
    "    \e[32m'fields'\e[0m: [                                            \n"
    "      { \e[32m'name'\e[0m:'density',\e[32m'expr'\e[0m: '<math>' }     \n"
    "      { \e[32m'name'\e[0m:'temperature', \e[32m'expr'\e[0m: '<math>' }\n"
    "    ]                                                                 \n"
    "  },                                                                  \n"
    "  \e[32m'parts'\e[0m: [                                               \n"
    "    {                                                                 \n"
    "      \e[32m'field'\e[0m: 'density',                                  \n"
    "      \e[32m'pairs'\e[0m: [                                           \n"
    "        {                                                             \n"
    "          \e[32m'uid'\e[0m: 1,                                        \n"
    "          \e[32m'source'\e[0m: <math>,                                \n"
    "          \e[32m'target'\e[0m: <math>                                 \n"
    "        },                                                            \n"
    "        {                                                             \n"
    "          \e[32m'uid'\e[0m: 2,                                        \n"
    "          \e[32m'source'\e[0m: <math>,                                \n"
    "          \e[32m'target'\e[0m: <math>                                 \n"
    "        }                                                             \n"
    "      ]                                                               \n"
    "    }                                                                 \n"
    "  ]                                                                   \n"
    " }                                                                    \n"
  );
}

/**
 * @brief Handle errors during argument parsing.
 *
 * @param message: a message to be displayed if any.
 * @param show_usage: hint for usage instructions print.
 * @return status
 */
bool abort(std::string message, bool show_usage = true) {

  if (show_usage)
    print_usage();

  std::fprintf(stderr,
    " ---------------------------------------------------------- \n"
    " \e[31m Error: %s. \e[0m                                    \n"
    " ---------------------------------------------------------- \n",
    message.data()
  );
  return false;
}

static app::Params params;

/**
 * @brief Parse and store app parameters.
 *
 * @param path: the JSON parameter file path.
 * @return parsing status flag.
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
              return abort("no source entities expression for part");
            // check target list
            if (not pair.count("target"))
              return abort("no target entities expression for part");
          }
        }
        helper.clear();
      }
    }

    using namespace app;

    // then store them
    params.dimension = json["mesh"]["dimension"];
    params.conformal = json["mesh"]["conformal"];
    params.dump      = json["mesh"]["export"];
    params.source    = json["mesh"]["source"];
    params.target    = json["mesh"]["target"];
    params.order     = json["remap"]["order"];
    params.kind      = json["remap"]["kind"] == "cell" ? entity::cell : entity::node;
    params.limiter   = json["remap"]["limiter"] ? limiter::none : limiter::barth;

    for (auto&& scalar : json["remap"]["fields"])
      params.fields[scalar["name"]] = scalar["expr"];

    for (auto&& entry : json["parts"]) {
      auto const field = entry["field"];
      for (auto&& pair : entry["entities"]) {
        params.parts[field].push_back(
          { pair["uid"], pair["source"], pair["target"] }
        );
      }
    }

    // check their validity eventually
    file.open(params.source);
    if (not file.good())
      return abort("unable to read source mesh file", false);
    file.close();

    file.open(params.target);
    if (not file.good())
      return abort("unable to read target mesh file", false);

    if (params.dimension < 2 or params.dimension > 3)
      return abort("invalid mesh dimension [2|3]", false);

    if (params.order < 1 or params.order > 2)
      return abort("invalid order of accuracy for remap [1|2]", false);

  } catch(nlohmann::json::parse_error& e) {
    return abort(e.what());
  }
  // everything was ok
  return true;
}

/**
 * @brief Check application parameters.
 *
 * @param argc arguments count
 * @param argv arguments values
 * @param my_rank current rank
 * @return status
 */
bool valid(int argc, char* argv[], int my_rank) {

  if (argc != 2) {
    if (my_rank == 0)
      print_usage();
    return false;
  }

  return parse_params(argv[1]);
}


template <int dim>
void remap(Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper,
           Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper,
           Wonton::Jali_State_Wrapper source_state_wrapper,
           Wonton::Jali_State_Wrapper target_state_wrapper) {

  // useful shortcuts
  using Remapper = Portage::CoreDriver<dim, Wonton::Entity_kind::CELL,
                                            Wonton::Jali_Mesh_Wrapper,
                                            Wonton::Jali_State_Wrapper>;

  using Parts = Portage::Parts<dim, Wonton::Entity_kind::CELL,
                                    Wonton::Jali_Mesh_Wrapper,
                                    Wonton::Jali_State_Wrapper>;

  // perform kernels
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto source_weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);


}


/**
 * @brief Run the application.
 *
 * @param argc: arguments count
 * @param argv: arguments values
 * @return status code
 */
int main(int argc, char* argv[]) {

  using namespace app;

  std::printf(" ---------------------------------------------------------- \n");
  std::printf("  Demonstration app for multi-part field interpolation.     \n");
  std::printf("  It handles pure cells remap of multi-material meshes and  \n");
  std::printf("  non-smoothed remap of fields with sharp discontinuities.  \n");
  std::printf(" ---------------------------------------------------------- \n");

  auto tic = timer::now();

  int threading = 0;
  int my_rank = 0;
  int nb_ranks = 0;
  MPI_Comm comm = MPI_COMM_WORLD;

  // init MPI
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threading);
  MPI_Comm_size(comm, &nb_ranks);
  MPI_Comm_rank(comm, &my_rank);

  // check and parse parameters
  if (not valid(argc, argv, my_rank)) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("Initializing distributed meshes ... ");

  // source and target meshes and states
  std::shared_ptr<Jali::Mesh>  source_mesh;
  std::shared_ptr<Jali::Mesh>  target_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::State> target_state;

  // load both distributed meshes
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.included_entities({Jali::Entity_kind::ALL_KIND});
  mesh_factory.partitioner(Jali::Partitioner_type::METIS);

  source_mesh  = mesh_factory(params.source);
  target_mesh  = mesh_factory(params.target);
  source_state = Jali::State::create(source_mesh);
  target_state = Jali::State::create(target_mesh);

  // interfaces with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  // ensure that source and target mesh have the same dimension,
  // and that it corresponds to the one specified by the user.
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(params.dimension == source_mesh->space_dimension());

  // retrieve mesh resolutions
  int const nb_source_cells = source_mesh_wrapper.num_owned_cells()
                            + source_mesh_wrapper.num_ghost_cells();

  int const nb_target_cells = source_mesh_wrapper.num_owned_cells()
                            + source_mesh_wrapper.num_ghost_cells();

  auto init_time = timer::elapsed(tic, true);
  if (my_rank == 0)
    std::printf(" done. \e[32m(%.3f s)\e[0m\n", init_time);

  MPI_Barrier(comm);

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("Running part-by-part remap ... \n");

  assert(not params.fields.empty());
  assert(not params.parts.empty());

  // print some infos for the user
  if (my_rank == 0) {
    std::printf(" = source mesh has %d cells.\n", nb_source_cells);
    std::printf(" = target mesh has %d cells.\n", nb_target_cells);
    std::printf(" = specified numerical fields: \n");
    for (auto&& field : params.fields)
      std::printf("   - %s: %s\n", field.first.data(), field.second.data());

    std::printf("\n");
  }

  // assign scalar fields
  for (auto&& field : params.fields) {
    double* field_data = nullptr;
    // add scalar field to both meshes
    source_state_wrapper.mesh_add_data<double>(entity::cell, field.first, 0.);
    target_state_wrapper.mesh_add_data<double>(entity::cell, field.first, 0.);
    source_state_wrapper.mesh_get_data(entity::cell, field.first, &field_data);

    // evaluate field expression and assign it
    user_field_t source_field;
    assert(source_field.initialize(params.dimension, field.second));
    for (int c = 0; c < nb_source_cells; c++) {
      field_data[c] = source_field(source_mesh->cell_centroid(c));
    }
  }

  // process remap








  MPI_Finalize();
  return EXIT_SUCCESS;
}