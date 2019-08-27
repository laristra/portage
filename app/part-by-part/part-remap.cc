/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <memory>
#include <chrono>
#include <map>
#include "mpi.h" // mandatory

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
#include "filter.h"
/* -------------------------------------------------------------------------- */
namespace entity {
  auto const cell = Wonton::Entity_kind::CELL;
  auto const node = Wonton::Entity_kind::NODE;

  // entities parts pair data
  struct part {
    int id;
    std::string source;
    std::string target;

    part(int in_id, std::string in_source, std::string in_target)
      : id(in_id), source(in_source), target(in_target)
    {}
  };
}
/* -------------------------------------------------------------------------- */
namespace limiter {
  auto const none  = Portage::Limiter_type::NOLIMITER;
  auto const barth = Portage::Limiter_type::BARTH_JESPERSEN;
}

/* -------------------------------------------------------------------------- */
namespace timer {
  // get rid of long namespaces
  using time_t = std::chrono::high_resolution_clock::time_point;

  // get current time point
  inline time_t now() { return std::chrono::high_resolution_clock::now(); }

  // retrieve elapsed time in seconds.
  inline float elapsed(time_t& tic, bool reset = false) {
    auto secs = static_cast<float>(
      std::chrono::duration_cast<std::chrono::seconds>(now() - tic).count()
    );

    if (reset)
      tic = now();

    return secs;
  }
}

/* -------------------------------------------------------------------------- */
struct Params {

  /* mesh */
  int dimension = 2;         // meshes dimension
  bool conformal = true;     // conformal meshes?
  bool dump = false;         // export results?
  std::string source {};     // source mesh file
  std::string target {};     // target mesh file

  /* remap */
  int order = 1;
  double tolerance = Portage::DEFAULT_CONSERVATION_TOL;
  Wonton::Entity_kind kind = entity::cell;
  Portage::Limiter_type limiter = limiter::none;
  std::map<std::string, std::string> fields {};    // fields expression
  std::map<std::string, std::vector<entity::part>> parts {}; // per-field

  /* fixups */
  int max_fix_iter = 5;
  Portage::Partial_fixup_type partial_fixup = Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE;
  Portage::Empty_fixup_type empty_fixup = Portage::Empty_fixup_type::EXTRAPOLATE;
};

/* -------------------------------------------------------------------------- */
// app parameters
static Params params;

// MPI variables
static int threading = 0;
static int my_rank   = 0;
static int nb_ranks  = 0;
static MPI_Comm comm = MPI_COMM_WORLD;

// constants for interpolation
static double const upper_bound = std::numeric_limits<double>::max();
static double const lower_bound = -upper_bound;
static double const epsilon = 1.E-10;

/* -------------------------------------------------------------------------- */
/**
 * @brief Display run instructions and input file format.
 *
 */
void print_usage();

/**
 * @brief Handle runtime errors.
 *
 * @param message: a message to be displayed if any.
 * @param show_usage: hint for usage instructions print.
 * @return status code
 */
int abort(std::string message, bool show_usage = true);

/**
 * @brief Parse and store app parameters.
 *
 * @param path: the JSON parameter file path.
 * @return parsing status flag.
 */
int parse(int argc, char* argv[]);

/**
 * @brief Process part-by-part remapping for given field.
 *
 * @tparam dim:
 * @param field
 * @param nb_parts
 * @param source_mesh
 * @param target_mesh
 * @param source_mesh_wrapper
 * @param target_mesh_wrapper
 * @param source_state_wrapper
 * @param target_state_wrapper
 * @param executor
 * @param source_cells
 * @param target_cells
 * @return status
 */
template <int dim>
bool remap(std::string field, int nb_parts,
           std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
           Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
           Wonton::Jali_State_Wrapper& source_state_wrapper,
           Wonton::Jali_State_Wrapper& target_state_wrapper,
           Wonton::Executor_type* executor,
           std::vector<std::vector<int>> const& source_cells,
           std::vector<std::vector<int>> const& target_cells);

/**
 * @brief Parse and store app parameters.
 *
 * @param path: the JSON parameter file path.
 * @return parsing status flag.
 */
int parse(int argc, char* argv[]);

/* -------------------------------------------------------------------------- */

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
    "    \e[32m'limiter'\e[0m: <boolean>,                                  \n"
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
    "          \e[32m'source'\e[0m: <math. predicate>,                     \n"
    "          \e[32m'target'\e[0m: <math. predicate>                      \n"
    "        },                                                            \n"
    "        {                                                             \n"
    "          \e[32m'uid'\e[0m: 2,                                        \n"
    "          \e[32m'source'\e[0m: <math. predicate>,                     \n"
    "          \e[32m'target'\e[0m: <math. predicate>                      \n"
    "        }                                                             \n"
    "      ]                                                               \n"
    "    }                                                                 \n"
    "  ]                                                                   \n"
    " }                                                                    \n"
  );
}

/**
 * @brief Handle runtime errors.
 *
 * @param message: a message to be displayed if any.
 * @param show_usage: hint for usage instructions print.
 * @return status
 */
int abort(std::string message, bool show_usage) {

  if (my_rank == 0) {
    if (show_usage)
      print_usage();

    std::fprintf(stderr,
      " ---------------------------------------------------------- \n"
      " \e[31m Error: %s. \e[0m                                    \n"
      " ---------------------------------------------------------- \n",
      message.data()
    );
  }

  MPI_Finalize();
  return EXIT_FAILURE;
}


/**
 * @brief Parse and store app parameters.
 *
 * @param path: the JSON parameter file path.
 * @return parsing status flag.
 */
int parse(int argc, char* argv[]) {

  if (argc != 2 or argv == nullptr)
    return abort("wrong arguments");

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

    // then store them
    params.dimension = json["mesh"]["dimension"];
    params.conformal = json["mesh"]["conformal"];
    params.dump      = json["mesh"]["export"];
    params.source    = json["mesh"]["source"];
    params.target    = json["mesh"]["target"];
    params.order     = json["remap"]["order"];
    params.kind      = json["remap"]["kind"] == "cell" ? entity::cell : entity::node;
    params.limiter   = json["remap"]["limiter"] ? limiter::none : limiter::barth;

    /* fixup */
    params.max_fix_iter     = json["remap"]["fixup"]["max-iter"];
    std::string parts_fixup = json["remap"]["fixup"]["partial"];
    std::string empty_fixup = json["remap"]["fixup"]["empty"];

    if (parts_fixup == "locally_conservative")
      params.partial_fixup = Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE;
    else if (parts_fixup == "constant")
      params.partial_fixup = Portage::Partial_fixup_type::CONSTANT;
    else if (parts_fixup == "shifted_conservative")
      params.partial_fixup = Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE;
    else
      return abort("unsupported partially filled cells fixup scheme");

    if (empty_fixup == "leave_empty")
      params.empty_fixup = Portage::Empty_fixup_type::LEAVE_EMPTY;
    else if (empty_fixup == "extrapolate")
      params.empty_fixup = Portage::Empty_fixup_type::EXTRAPOLATE;
    else
      return abort("unsupported empty cells fixup scheme");

    /* parts field */
    for (auto&& scalar : json["remap"]["fields"])
      params.fields[scalar["name"]] = scalar["expr"];

    for (auto&& entry : json["parts"]) {
      std::string field = entry["field"];
      for (auto&& pair : entry["pairs"]) {
        int id = pair["uid"];
        std::string source = pair["source"];
        std::string target = pair["target"];
        params.parts[field].emplace_back(id, source, target);
      }
      assert(not params.parts[field].empty());
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
  return EXIT_SUCCESS;
}

/**
 *
 * @param field
 * @param nb_parts
 * @param source_mesh
 * @param target_mesh
 * @param source_mesh_wrapper
 * @param target_mesh_wrapper
 * @param source_state_wrapper
 * @param target_state_wrapper
 * @param executor
 * @param source_cells
 * @param target_cells
 * @return
 */
template<>
bool remap<2>(std::string field, int nb_parts,
              std::shared_ptr<Jali::Mesh> source_mesh,
              std::shared_ptr<Jali::Mesh> target_mesh,
              Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
              Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
              Wonton::Jali_State_Wrapper& source_state_wrapper,
              Wonton::Jali_State_Wrapper& target_state_wrapper,
              Wonton::Executor_type* executor,
              std::vector<std::vector<int>> const& source_cells,
              std::vector<std::vector<int>> const& target_cells) {

  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;

  using Parts = Portage::Parts<2, Wonton::Entity_kind::CELL,
                                  Wonton::Jali_Mesh_Wrapper,
                                  Wonton::Jali_State_Wrapper>;

  std::vector<Parts> parts_manager;
  parts_manager.reserve(nb_parts);

  // filter cells and populate lists
  for (int i = 0; i < nb_parts; ++i) {
    parts_manager.emplace_back(source_mesh_wrapper, target_mesh_wrapper,
                               source_state_wrapper,target_state_wrapper,
                               source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  for (int i = 0; i < nb_parts; ++i) {
    // test for mismatch and compute volumes
    parts_manager[i].test_mismatch(weights);

    // interpolate field for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      field, field, weights, lower_bound, upper_bound,
      params.limiter, params.partial_fixup,
      params.empty_fixup, params.tolerance,
      params.max_fix_iter, &(parts_manager[i])
    );
  }
}

/**
 *
 * @param field
 * @param nb_parts
 * @param source_mesh
 * @param target_mesh
 * @param source_mesh_wrapper
 * @param target_mesh_wrapper
 * @param source_state_wrapper
 * @param target_state_wrapper
 * @param executor
 * @param source_cells
 * @param target_cells
 * @return
 */
template<>
bool remap<3>(std::string field, int nb_parts,
              std::shared_ptr<Jali::Mesh> source_mesh,
              std::shared_ptr<Jali::Mesh> target_mesh,
              Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
              Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
              Wonton::Jali_State_Wrapper& source_state_wrapper,
              Wonton::Jali_State_Wrapper& target_state_wrapper,
              Wonton::Executor_type* executor,
              std::vector<std::vector<int>> const& source_cells,
              std::vector<std::vector<int>> const& target_cells) {

  using Remapper = Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;

  using Parts = Portage::Parts<3, Wonton::Entity_kind::CELL,
                                  Wonton::Jali_Mesh_Wrapper,
                                  Wonton::Jali_State_Wrapper>;

  std::vector<Parts> parts_manager;
  parts_manager.reserve(nb_parts);

  // filter cells and populate lists
  for (int i = 0; i < nb_parts; ++i) {
    parts_manager.emplace_back(source_mesh_wrapper, target_mesh_wrapper,
                               source_state_wrapper,target_state_wrapper,
                               source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  for (int i = 0; i < nb_parts; ++i) {
    // test for mismatch and compute volumes
    parts_manager[i].test_mismatch(weights);

    // interpolate field for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      field, field, weights, lower_bound, upper_bound,
      params.limiter, params.partial_fixup,
      params.empty_fixup, params.tolerance,
      params.max_fix_iter, &(parts_manager[i])
    );
  }
}

/**
 * @brief Run the application.
 *
 * @param argc: arguments count
 * @param argv: arguments values
 * @return status code
 */
int main(int argc, char* argv[]) {

  auto tic = timer::now();

  // init MPI
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threading);
  MPI_Comm_size(comm, &nb_ranks);
  MPI_Comm_rank(comm, &my_rank);

  if (my_rank == 0) {
    std::printf(" ---------------------------------------------------------- \n");
    std::printf("  Demonstration app for multi-part field interpolation.     \n");
    std::printf("  It handles pure cells remap of multi-material meshes and  \n");
    std::printf("  non-smoothed remap of fields with sharp discontinuities.  \n");
    std::printf(" ---------------------------------------------------------- \n");
  }

  // check and parse parameters
  if (parse(argc, argv) == EXIT_FAILURE)
    return EXIT_FAILURE;

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

  Wonton::MPIExecutor_type mpi_executor(comm);
  Wonton::Executor_type *executor = (nb_ranks > 1 ? &mpi_executor : nullptr);

  // ensure that source and target mesh have the same dimension,
  // and that it corresponds to the one specified by the user.
  assert(source_mesh->space_dimension() == target_mesh->space_dimension());
  assert(params.dimension == source_mesh->space_dimension());

  // retrieve mesh resolutions
  long const nb_source_cells = source_mesh_wrapper.num_owned_cells();
  long const nb_target_cells = target_mesh_wrapper.num_owned_cells();
  int const nb_fields = params.fields.size();

  MPI_Barrier(comm);

  auto init_time = timer::elapsed(tic, true);
  if (my_rank == 0)
    std::printf(" done. \e[32m(%.3f s)\e[0m\n", init_time);

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("Running part-by-part remap ... \n");

  long total_source_cells = 0;
  long total_target_cells = 0;
  MPI_Reduce(&nb_source_cells, &total_source_cells, 1, MPI_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(&nb_target_cells, &total_target_cells, 1, MPI_LONG, MPI_SUM, 0, comm);

  // print some infos for the user
  if (my_rank == 0) {
    std::printf(" - source mesh has %ld cells.\n", total_source_cells);
    std::printf(" - target mesh has %ld cells.\n", total_target_cells);
    std::printf(" - specified numerical fields: \n");
    for (auto&& field : params.fields)
      std::printf("   \u2022 %s: %s\n", field.first.data(), field.second.data());

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
    if(source_field.initialize(params.dimension, field.second)) {
      for (long c = 0; c < nb_source_cells; c++) {
        field_data[c] = source_field(source_mesh->cell_centroid(c));
        #if DEBUG_PART_BY_PART
          std::printf("field_data[%d]: %.3f\n", c, field_data[c]);
        #endif
      }
    } else
      return abort("cannot parse numerical field "+ field.second, false);
  }

  MPI_Barrier(comm);

  std::vector<std::vector<int>> source_cells;
  std::vector<std::vector<int>> target_cells;

  // remap each field
  for (auto&& entry : params.fields) {
    std::string field = entry.first;
    int const nb_parts = params.parts[field].size();
    assert(nb_parts > 1);

    // first determine source and target parts for this field.
    filter_t filter;
    source_cells.resize(nb_parts);
    target_cells.resize(nb_parts);

    if (my_rank == 0)
      std::printf(" - Field '%s' (%d parts):\n", field.data(), nb_parts);

    for (int i = 0; i < nb_parts; ++i) {
      auto const& part = params.parts[field][i];
      source_cells[i].reserve(nb_source_cells);
      target_cells[i].reserve(nb_target_cells);

      // populate source part entities
      if (filter.initialize(params.dimension, part.source)) {
        for (long s = 0; s < nb_source_cells; ++s) {
          if (filter(source_mesh->cell_centroid(s)))
            source_cells[i].push_back(s);
        }
        source_cells[i].shrink_to_fit();
      }
      else
        return abort("cannot filter source part cells for field "+field, false);

      // populate target part entities
      if (filter.initialize(params.dimension, part.target)) {
        for (long t = 0; t < nb_target_cells; ++t) {
          if (filter(target_mesh->cell_centroid(t)))
            target_cells[i].push_back(t);
        }
        target_cells[i].shrink_to_fit();
      }
      else
        return abort("cannot filter target part cells for field "+field, false);

      long local_source_part = source_cells[i].size();
      long local_target_part = target_cells[i].size();
      long total_source_part = 0;
      long total_target_part = 0;
      MPI_Reduce(&local_source_part, &total_source_part, 1, MPI_LONG, MPI_SUM, 0, comm);
      MPI_Reduce(&local_target_part, &total_target_part, 1, MPI_LONG, MPI_SUM, 0, comm);

      if (my_rank == 0) {
        std::printf(
          "   \u2022 part[%d]: source: %ld cells, target: %ld cells\n",
          i, total_source_part, total_target_part
        );
      }
    }

    MPI_Barrier(comm);

    // then process part-by-part remapping.
    // need to explicitly instantiate the driver due to template arguments
    // forwarding failure when instantiating search and intersection kernels.
    // part-by-part remap each field

    switch (params.dimension) {
      case 2: remap<2>(field, nb_parts, source_mesh, target_mesh,
                       source_mesh_wrapper, target_mesh_wrapper,
                       source_state_wrapper, target_state_wrapper,
                       executor, source_cells, target_cells); break;

      case 3: remap<3>(field, nb_parts, source_mesh, target_mesh,
                       source_mesh_wrapper, target_mesh_wrapper,
                       source_state_wrapper, target_state_wrapper,
                       executor, source_cells, target_cells); break;

      default: return abort("invalid dimension", false);
    }

    source_cells.clear();
    target_cells.clear();
    MPI_Barrier(comm);
  }

  auto remap_time = timer::elapsed(tic, true);
  if (my_rank == 0)
    std::printf("Remap done. \e[32m(%.3f s)\e[0m\n", remap_time);

  /* ------------------------------------------------------------------------ */
  if (params.dump) {
    if (my_rank == 0)
      std::printf("Dump data to exodus files ... ");

    // all field data are already attached to meshes states
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source.exo");
    target_mesh->write_to_exodus_file("target.exo");

    MPI_Barrier(comm);

    auto dump_time = timer::elapsed(tic, true);
    if (my_rank == 0)
      std::printf(" done. \e[32m(%.3f s)\e[0m\n", dump_time);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */