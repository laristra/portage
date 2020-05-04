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
#include <portage/interpolate/interpolate_2nd_order.h>
#include "mpi.h" // mandatory

// meshes and states
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "MeshFactory.hh"

// remap kernels and driver
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/driver/coredriver.h"
#include "portage/support/mpi_collate.h"

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
namespace timer {
  // get rid of long namespaces
  using time_t = std::chrono::high_resolution_clock::time_point;

  // get current time point
  inline time_t now() { return std::chrono::high_resolution_clock::now(); }

  // retrieve elapsed time in seconds.
  inline float elapsed(time_t& tic, bool reset = true) {
    auto const toc = now();
    auto secs = static_cast<float>(
      std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()
    ) * 1.E-3;

    if (reset)
      tic = now();

    return secs;
  }
}

/* -------------------------------------------------------------------------- */
struct Params {

  /* mesh */
  int dimension = 2;       
  bool conformal = true;   
  bool dump = false;       
  std::string source {};   
  std::string target {};   

  /* remap */
  int order = 1;
  double tolerance = Portage::DEFAULT_NUMERIC_TOLERANCES<2>.relative_conservation_eps;
  Wonton::Entity_kind kind = Wonton::Entity_kind::CELL;
  std::map<std::string, std::string> fields {};
  std::map<std::string, std::vector<entity::part>> parts {};

  /* fixups */
  int fix_iter = 5;
  Portage::Limiter_type limiter = Portage::Limiter_type::NOLIMITER;
  Portage::Boundary_Limiter_type bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;
  Portage::Empty_fixup_type empty_fixup = Portage::Empty_fixup_type::EXTRAPOLATE;
  Portage::Partial_fixup_type partial_fixup =
    Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE;
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
 * @brief Display input parameters.
 *
 */
void print_params(nlohmann::json const& json);

/**
 * @brief Get number of digits of the given scalar.
 *
 * @tparam type_t: scalar type
 * @param number: scalar value
 * @return its number of digits for print.
 */
template<typename type_t>
int get_number_digit(type_t number);

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
 * @brief Process part-by-part remapping for the given field.
 *
 * @tparam dim: source/target meshes dimension.
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 */
template <int dim>
void remap(std::string field, int nb_parts,
           std::shared_ptr<Jali::Mesh> source_mesh,
           std::shared_ptr<Jali::Mesh> target_mesh,
           Wonton::Jali_Mesh_Wrapper&  source_mesh_wrapper,
           Wonton::Jali_Mesh_Wrapper&  target_mesh_wrapper,
           Wonton::Jali_State_Wrapper& source_state_wrapper,
           Wonton::Jali_State_Wrapper& target_state_wrapper,
           Wonton::Executor_type* executor,
           std::vector<std::vector<int>> const& source_cells,
           std::vector<std::vector<int>> const& target_cells);

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
    "    \e[32m'bnd_limiter'\e[0m: <boolean>,                              \n"    
    "    \e[32m'fixup'\e[0m: {                                             \n"
    "      \e[32m'partial'\e[0m: 'constant|<locally|shifted>_conservative',\n"
    "      \e[32m'empty'\e[0m: '<leave_empty|extrapolate>',                \n"
    "      \e[32m'max-iter'\e[0m: <unsigned integer>,                      \n"
    "    },                                                                \n"
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
 * @brief Display input parameters.
 *
 */
void print_params(nlohmann::json const& json) {

  // extract file path base name
  auto base = [](const std::string& path) {
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

  std::string remap_kind  = json["remap"]["kind"];
  std::string dump_result = json["mesh"]["export"] ? "yes": "no";
  std::string use_limiter = json["remap"]["limiter"] ? "yes": "no";
  std::string use_bnd_limiter = json["remap"]["bnd_limiter"] ? "yes": "no";
  std::string partial_fix = json["remap"]["fixup"]["partial"];
  std::string empty_fix   = json["remap"]["fixup"]["empty"];

  std::printf("Parameters: \n");
  std::printf(" \u2022 MPI ranks: \e[32m%d\e[0m\n", nb_ranks);
  std::printf(" \u2022 dimension: \e[32m%d\e[0m\n", params.dimension);
  std::printf(" \u2022 source mesh: '\e[32m%s\e[0m'\n", base(params.source).data());
  std::printf(" \u2022 target mesh: '\e[32m%s\e[0m'\n", base(params.target).data());
  std::printf(" \u2022 dump results: \e[32m%s\e[0m\n", dump_result.data());
  std::printf(" \u2022 remap order: \e[32m%d\e[0m\n", params.order);
  std::printf(" \u2022 remap kind: \e[32m%s\e[0m\n", remap_kind.data());
  std::printf(" \u2022 use limiter: \e[32m%s\e[0m\n", use_limiter.data());
  std::printf(" \u2022 use boundary limiter: \e[32m%s\e[0m\n", use_bnd_limiter.data());
  std::printf(" \u2022 partial filled fixup: \e[32m%s\e[0m\n", partial_fix.data());
  std::printf(" \u2022 empty cells fixup: \e[32m%s\e[0m\n", empty_fix.data());
  std::printf(" \u2022 fixup iterations: \e[32m%d\e[0m\n", params.fix_iter);
  std::printf(" \u2022 numerical fields: \n");

  for (auto&& field: params.fields) {
    std::printf("   - %s: '\e[32m%s\e[0m'\n", field.first.data(), field.second.data());

    for (auto&& part : params.parts[field.first]) {
      std::printf("   - [part: \e[32m%d\e[0m,", part.id);
      std::printf(" source: '\e[32m%s\e[0m',", part.source.data());
      std::printf(" target: '\e[32m%s\e[0m']\n", part.target.data());
    }
    std::printf("\n");
  }
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
    params.fix_iter  = json["remap"]["fixup"]["max-iter"];

    bool const use_limiter  = json["remap"]["limiter"];
    bool const use_bnd_limiter  = json["remap"]["bnd_limiter"];
    std::string remap_kind  = json["remap"]["kind"];
    std::string parts_fixup = json["remap"]["fixup"]["partial"];
    std::string empty_fixup = json["remap"]["fixup"]["empty"];

    if (remap_kind == "cell")
      params.kind = Wonton::Entity_kind::CELL;
    else
      params.kind = Wonton::Entity_kind::NODE;

    if (use_limiter)
      params.limiter = Portage::Limiter_type::BARTH_JESPERSEN;
    else
      params.limiter = Portage::Limiter_type::NOLIMITER;

    if (use_bnd_limiter)
      params.bnd_limiter = Portage::Boundary_Limiter_type::BND_ZERO_GRADIENT;
    else
      params.bnd_limiter = Portage::Boundary_Limiter_type::BND_NOLIMITER;

    if (parts_fixup == "locally_conservative")
      params.partial_fixup = Portage::Partial_fixup_type::LOCALLY_CONSERVATIVE;
    else if (parts_fixup == "constant")
      params.partial_fixup = Portage::Partial_fixup_type::CONSTANT;
    else
      params.partial_fixup = Portage::Partial_fixup_type::SHIFTED_CONSERVATIVE;

    if (empty_fixup == "leave_empty")
      params.empty_fixup = Portage::Empty_fixup_type::LEAVE_EMPTY;
    else
      params.empty_fixup = Portage::Empty_fixup_type::EXTRAPOLATE;

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

    if (params.kind == Wonton::Entity_kind::NODE)
      return abort("multi-part node remap is not supported", false);

    if (params.order < 1 or params.order > 2)
      return abort("invalid order of accuracy for remap [1|2]", false);

    // check that params.fields and params.parts have same keys
    auto same_keys = [](auto const& a, auto const& b) { return a.first == b.first; };

    bool have_same_size = params.fields.size() == params.parts.size();
    bool have_same_keys = std::equal(params.fields.begin(), params.fields.end(),
                                     params.parts.begin(), same_keys);

    if (not have_same_size or not have_same_keys)
      return abort("numerical fields and per-part fields mismatch");


  } catch(nlohmann::json::parse_error& e) {
    return abort(e.what());
  }

  // everything was ok
  if (my_rank == 0)
    print_params(json);

  return EXIT_SUCCESS;
}


/**
 * @brief Process part-by-part remapping for the given 2D field.
 *
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 */
template<>
void remap<2>(std::string field, int nb_parts,
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

  using PartPair = Portage::PartPair<2, Wonton::Jali_Mesh_Wrapper,
                                        Wonton::Jali_State_Wrapper>;


  std::vector<PartPair> parts_manager;
  parts_manager.reserve(nb_parts);

  for (int i = 0; i < nb_parts; ++i) {
    // create source-target mesh parts manager and
    // populate cell lists for the current part.
    parts_manager.emplace_back(source_mesh_wrapper, source_state_wrapper,
                               target_mesh_wrapper, target_state_wrapper,
                               source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  // use the right interpolator according to the requested order of remap.
  auto interpolate = [&](auto* current_part) {
    Portage::vector<Wonton::Vector<2>> *gradients = nullptr;
    auto const source_part = current_part->source();

    switch (params.order) {
      case 1: 
      
        remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
          field, field, weights, current_part
        );
              
      break;

      case 2: *gradients = remapper.compute_source_gradient(field, params.limiter,
                                                            params.bnd_limiter,0,
                                                            &source_part);

        remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
          field, field, weights, current_part, gradients
        );
        
      break;

      default: throw std::runtime_error("wrong remap order");
    }
    
    if (current_part->has_mismatch())
      current_part->fix_mismatch(field, field, lower_bound, upper_bound, 
                                 params.tolerance, params.fix_iter,
                                 params.partial_fixup, params.empty_fixup);

  };

  for (int i = 0; i < nb_parts; ++i) {
  
    // compute volumes of intersection and test for parts boundaries mismatch.
    parts_manager[i].check_mismatch(weights);

    // interpolate field for each part and fix partially filled or empty cells.
    interpolate(parts_manager.data() + i);
  }
}

/**
 * Process part-by-part remapping for the given 3D field.
 *
 * @param field: a string expression of the numerical field.
 * @param nb_parts: the number of source-target parts couples.
 * @param source_mesh: a pointer to the source mesh.
 * @param target_mesh: a pointer to the target mesh.
 * @param source_mesh_wrapper: a wrapper to access source mesh data.
 * @param target_mesh_wrapper: a wrapper to access target mesh data.
 * @param source_state_wrapper: a wrapper to access source state data.
 * @param target_state_wrapper: a wrapper to access source state data.
 * @param executor: a pointer to the MPI executor.
 * @param source_cells: list of source cells for each part.
 * @param target_cells: list of target cells for each part.
 */
template<>
void remap<3>(std::string field, int nb_parts,
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

  using PartPair = Portage::PartPair<3, Wonton::Jali_Mesh_Wrapper,
                                        Wonton::Jali_State_Wrapper>;

  std::vector<PartPair> parts_manager;
  parts_manager.reserve(nb_parts);

  // filter cells and populate lists
  for (int i = 0; i < nb_parts; ++i) {
    parts_manager.emplace_back(source_mesh_wrapper, source_state_wrapper,
                               target_mesh_wrapper, target_state_wrapper,
                               source_cells[i], target_cells[i], executor);
  }

  // perform remap kernels.
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR3D>(candidates);

  // use the right interpolator according to the requested order of remap.
  auto interpolate = [&](auto* current_part) {
    Portage::vector<Wonton::Vector<3>> gradients;

    switch (params.order) {
      case 1: 
      
        remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
          field, field, weights, current_part
        );
              
      case 2: 
      
        gradients = remapper.compute_source_gradient(field, params.limiter,
                                                     params.bnd_limiter,0,
                                                     &(current_part->source()));

        remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
          field, field, weights, current_part, &gradients
        );
        
      break;

      default: throw std::runtime_error("wrong remap order");
    }
    
    if (current_part->has_mismatch())
      current_part->fix_mismatch(field, field, lower_bound, upper_bound, 
      params.tolerance, params.fix_iter,
      params.partial_fixup, params.empty_fixup);
  };

  for (int i = 0; i < nb_parts; ++i) {
    // compute volumes of intersection and test for parts boundaries mismatch.
    parts_manager[i].check_mismatch(weights);

    // interpolate field for each part and fix partially filled or empty cells.
    interpolate(parts_manager.data() + i);
  }
}

/**
 * @brief Get number of digits of the given scalar.
 *
 * @tparam type_t: scalar type
 * @param number: scalar value
 * @return its number of digits for print.
 */
template<typename type_t>
int get_number_digit(type_t number) {
  return (number > 0 ? static_cast<int>(std::floor(std::log10(number))) + 1 : 0);
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
  mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);
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
  assert(unsigned(params.dimension) == source_mesh->space_dimension());

  // retrieve mesh resolutions
  // nb: only consider owned cells for source mesh to avoid errors.
  long const nb_source_cells = source_mesh_wrapper.num_owned_cells();
  long const nb_target_cells = target_mesh_wrapper.num_owned_cells();

  MPI_Barrier(comm);

  if (my_rank == 0)
    std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));

  /* ------------------------------------------------------------------------ */
  if (my_rank == 0)
    std::printf("\nRunning part-by-part remap ... \n");

  long total_source_cells = 0;
  long total_target_cells = 0;
  MPI_Reduce(&nb_source_cells, &total_source_cells, 1, MPI_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(&nb_target_cells, &total_target_cells, 1, MPI_LONG, MPI_SUM, 0, comm);

  // for formatting
  int format = get_number_digit(std::max(total_source_cells, total_target_cells));

  // print some infos for the user
  if (my_rank == 0) {
    std::printf(" - source mesh has %*ld cells.\n", format, total_source_cells);
    std::printf(" - target mesh has %*ld cells.\n", format, total_target_cells);
    std::printf(" - specified numerical fields: \n");
    for (auto&& field : params.fields)
      std::printf("   \u2022 %s: \e[32m%s\e[0m\n", field.first.data(), field.second.data());

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
      std::printf(" - Remapping \e[32m%s\e[0m field on %d parts:\n", field.data(), nb_parts);

    long max_source_parts = 0;
    long max_target_parts = 0;
    long total_source_part[nb_parts];
    long total_target_part[nb_parts];
    std::fill(total_source_part, total_source_part + nb_parts, 0);
    std::fill(total_target_part, total_target_part + nb_parts, 0);

    // Filter source and target cells for each part
    // with respect to user-defined predicates.
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
      MPI_Reduce(&local_source_part, total_source_part + i, 1, MPI_LONG, MPI_SUM, 0, comm);
      MPI_Reduce(&local_target_part, total_target_part + i, 1, MPI_LONG, MPI_SUM, 0, comm);

      max_source_parts = std::max(total_source_part[i], max_source_parts);
      max_target_parts = std::max(total_target_part[i], max_target_parts);
    }

    MPI_Barrier(comm);

    if (my_rank == 0) {
      int const source_digits = get_number_digit(max_source_parts);
      int const target_digits = get_number_digit(max_target_parts);

      for (int i = 0; i < nb_parts; ++i) {
        std::printf(
          "   \u2022 ["
          " source: %*ld cells \e[32m(%5.2f %%)\e[0m,"
          " target: %*ld cells \e[32m(%5.2f %%)\e[0m ]\n",
          source_digits, total_source_part[i],
          100. * double(total_source_part[i]) / double(total_source_cells),
          target_digits, total_target_part[i],
          100. * double(total_target_part[i]) / double(total_target_cells)
        );
      }
      std::printf("\n");
    }

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

    if (my_rank == 0) {
      std::fflush(stdout);
      std::printf(" %s \n", std::string(58,'-').data());
    }
  }

  if (my_rank == 0)
    std::printf("Remap done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));

  /* ------------------------------------------------------------------------ */
  if (params.kind == Wonton::Entity_kind::CELL) {
    // compute error for each field
    for (auto&& field : params.fields) {

      double error_l1 = 0;
      double error_l2 = 0;
      double relative_error = 0;
      user_field_t exact_value;

      if (my_rank == 0)
        std::printf("\nComputing error for \e[32m%s\e[0m field ... ", field.first.data());

      if (exact_value.initialize(params.dimension, field.second)) {
        // keep track of min and max field values on both meshes
        double min_source_val = std::numeric_limits<double>::max();
        double max_source_val = std::numeric_limits<double>::min();
        double min_target_val = min_source_val;
        double max_target_val = max_source_val;
        double source_extents[] = { min_source_val, max_source_val };
        double target_extents[] = { min_target_val, max_target_val };
        double source_mass = 0;
        double target_mass = 0;
        double total_mass[] = { 0, 0 };
        double global_error[] = { 0, 0, 0 };

        // retrieve field data on both meshes
        double* source_field_data = nullptr;
        double* target_field_data = nullptr;
        source_state_wrapper.mesh_get_data(entity::cell, field.first, &source_field_data);
        target_state_wrapper.mesh_get_data(entity::cell, field.first, &target_field_data);

        // compute total mass on the source mesh to check conservation
        for (int s = 0; s < nb_source_cells; ++s) {
          min_source_val = std::min(source_field_data[s], min_source_val);
          max_source_val = std::max(source_field_data[s], max_source_val);
          source_mass += source_field_data[s] * source_mesh_wrapper.cell_volume(s);
        }

        // compute cell error
        for (int t = 0; t < nb_target_cells; ++t) {
          min_target_val = std::min(target_field_data[t], min_target_val);
          max_target_val = std::max(target_field_data[t], max_target_val);
          // compute difference between exact and remapped value
          auto const& centroid = target_mesh->cell_centroid(t);
          auto const error = exact_value(centroid) - target_field_data[t];
          auto const cell_volume = target_mesh_wrapper.cell_volume(t);
          // update L^p norm error and target mass
          error_l1 += std::abs(error) * cell_volume;
          error_l2 += error * error * cell_volume;
          relative_error += std::abs(target_field_data[t]) * cell_volume;
          target_mass += target_field_data[t] * cell_volume;
        }

        // accumulate all local values on rank 0
        MPI_Reduce(&error_l1, global_error, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&error_l2, global_error+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&relative_error, global_error+2, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&min_source_val, source_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&max_source_val, source_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&min_target_val, target_extents, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&max_target_val, target_extents+1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&source_mass, total_mass, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&target_mass, total_mass+1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

        // update local values then
        error_l1 = global_error[0];
        error_l2 = sqrt(global_error[1]);
        relative_error = error_l1 / global_error[2];
        source_mass = total_mass[0];
        target_mass = total_mass[1];

        MPI_Barrier(comm);

        if (my_rank == 0) {
          std::printf( " done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
          std::printf(" \u2022 L1-norm error     = %lf\n", error_l1);
          std::printf(" \u2022 L2-norm error     = %lf\n", error_l2);
          std::printf(" \u2022 relative L1 error = %.15f\n", relative_error);
          std::printf(" \u2022 source values     = [%.15f, %.15f]\n",
            source_extents[0], source_extents[1]);
          std::printf(" \u2022 target values     = [%.15f, %.15f]\n",
            target_extents[0], target_extents[1]);
          std::printf(" \u2022 source total mass = %.15f\n", source_mass);
          std::printf(" \u2022 target total mass = %.15f\n", target_mass);
          std::printf(" \u2022 mass discrepancy  = %.15f\n",
            std::abs(source_mass - target_mass));
        }
      } else
        return abort("cannot parse numerical field "+ field.second, false);
    }
  } else
    return abort("part-by-part node remap is not supported", false);

  /* ------------------------------------------------------------------------ */
  if (params.dump) {

    MPI_Barrier(comm);

    if (my_rank == 0)
      std::printf("\nDump data to exodus files ... ");

    // dump meshes with attached data
    source_state->export_to_mesh();
    target_state->export_to_mesh();
    source_mesh->write_to_exodus_file("source.exo");
    target_mesh->write_to_exodus_file("target.exo");

    MPI_Barrier(comm);

    std::vector<int> index_helper;
    std::vector<int> local_index;
    std::vector<int> global_index;
    std::vector<double> local_field;
    std::vector<double> global_field;
    double* field_data = nullptr;

    local_index.resize(nb_target_cells);
    local_field.resize(nb_target_cells);

    // dump each field in a separate file
    for (auto&& field : params.fields) {
      // retrieve cell global indices and field values for current rank
      for (auto c = 0; c < nb_target_cells; ++c)
        local_index[c] = target_mesh->GID(c, Jali::Entity_kind::CELL);

      target_state_wrapper.mesh_get_data(entity::cell, field.first, &field_data);
      std::copy(field_data, field_data + nb_target_cells, local_field.begin());

      // append local index and values lists to master rank global lists
      Portage::collate(comm, my_rank, nb_ranks, local_index, global_index);
      Portage::collate(comm, my_rank, nb_ranks, local_field, global_field);

      if (my_rank == 0) {
        // sort field values by global ID
        Portage::argsort(global_index, index_helper);
        Portage::reorder(global_index, index_helper);
        Portage::reorder(global_field, index_helper);

        // dump sorted data eventually
        std::ofstream file("remap_"+ field.first +".dat");
        file << std::scientific;
        file.precision(17);

        for (long c = 0; c < nb_target_cells; ++c)
          file << global_index[c] << "\t" << global_field[c] << std::endl;

        index_helper.clear();
        global_index.clear();
        global_field.clear();
      }

      MPI_Barrier(comm);
    }

    if (my_rank == 0)
      std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
