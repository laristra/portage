/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
*/

#pragma once

#include <fstream>
#include "mpi.h"

// parsers
#include "json.h"

#include "wonton/support/wonton.h"
#include "portage/support/portage.h"
#include "portage/accumulate/accumulate.h"

using Wonton::Point;

/**
 * @enum Coordinate frame.
 *
 * It specifies the coordinate frame to consider when generating the
 * source or target points. Here the global cartesian frame corresponds
 * to absolute coordinates (x,y) whereas the local cartesian frame
 * are relative coordinates (x',y') with respect to a given reference
 * point. The local frame is meant to be used only for generating points
 * from a radial grid.
 */
enum Frame : int { GLOBAL, LOCAL };

/**
 * @enum Point distribution type.
 *
 * It specifies the distribution to consider when generating the points.
 * The latter can be issued from a cartesian grid, a radial grid or
 * randomly using a uniform distribution.
 */
enum Distrib : int { CARTESIAN, RADIAL, RANDOM };

class Params {
public:

  struct {
    std::string file;              /** source points file if any */
    std::string field;             /** source field expression */
    Distrib distrib = CARTESIAN;   /** distribution for random points */
    Point<2> min, max;             /** coordinates extents */
    double span = 0;               /** span angle of radial mesh */
    double width = 0;              /** width of radial mesh */
    double radius = 0;             /** radius of radial mesh */
    Point<2> center;               /** center of radial mesh */
    Frame frame = GLOBAL;          /** coordinates frame */
    int size[2] = { 0, 0 };        /** number of points in each direction */
  } source;

  struct {
    std::string file;              /** target points file if any */
    Distrib distrib = CARTESIAN;   /** distribution for random points */
    Point<2> min, max;             /** coordinates extents */
    double span = 0;               /** span angle of radial mesh */
    double width = 0;              /** width of radial mesh */
    double radius = 0;             /** radius of radial mesh */
    Point<2> center;               /** center of radial mesh */
    Frame frame = GLOBAL;          /** coordinates frame */
    int size[2] = { 0, 0 };        /** number of points in each direction */
    std::string output;            /** output file for target points */
  } target;

  struct {
    using Matrix = std::vector<std::vector<double>>;
    Wonton::vector<Matrix> smooth_length;                           /** smoothing lengths */
    Portage::basis::Type basis = Portage::basis::Unitary;           /** weights basis type */
    Portage::WeightCenter weight = Portage::Gather;                 /** weights center form */
    Portage::Weight::Kernel kernel = Portage::Weight::B4;           /** weights kernel type */
    Portage::Weight::Geometry support = Portage::Weight::ELLIPTIC;  /** support geometry */
    Portage::EstimateType estimator = Portage::LocalRegression;     /** estimator type */
  } remap;

  struct {
    Point<2> coords = { 1, 1 };    /** scaling of coordinates */
    Point<2> smooth = { 1, 1 };    /** scaling of smoothing lengths */
    double field = 1;                    /** scaling of source field */
  } scale;

  struct {
    std::string source = "source.dat";  /** source field output */
    std::string exact = "exact.dat";    /** exact values output */
    std::string remap = "remap.dat";    /** remmapped field output */
    std::string error = "error.dat";    /** error map output */
  } output;

  struct {
    int rank = 0;                                /** current rank */
    int num_ranks = 1;                           /** number of ranks */
    MPI_Comm comm = MPI_COMM_WORLD;              /** communicator */
    Wonton::MPIExecutor_type executor { comm };  /** executor type */
  } mpi;

  /**
   * @brief Parse and initialize parameters.
   *
   * @param argc: number of arguments.
   * @param argv: list of arguments.
   * @return parsing status flag.
   */
  bool parse(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(mpi.comm, &mpi.num_ranks);
    MPI_Comm_rank(mpi.comm, &mpi.rank);
    mpi.executor.mpicomm = mpi.comm;

    if (argc != 2 or argv == nullptr) {
      return abort("wrong arguments");
    } else if (std::string(argv[1]) == "-h") {
      return abort();
    }

    std::ifstream file(argv[1]);
    if (not file.good()) {
      return abort("unable to open input file", false);
    }

    nlohmann::json json;

    try {

      file >> json;
      file.close();

      // verify and store attributes
      // source
      if (json.count("source")) {
        if (json["source"].count("size")) {
          int k = 0;
          for (auto&& count : json["source"]["size"]) { source.size[k++] = count; }

          if (json["source"].count("distrib")) {
            std::string distrib = json["source"]["distrib"];
            if (distrib == "cartesian")   { source.distrib = CARTESIAN; }
            else if (distrib == "radial") { source.distrib = RADIAL; }
            else if (distrib == "random") { source.distrib = RANDOM; }
            else { return abort("invalid source points distribution"); }
          } else { return abort("missing source points distribution"); }

          if (source.distrib == RADIAL) {
            if (json["source"].count("span")) {
              source.span = json["source"]["span"];
            } else { return abort("missing span angle for radial source mesh"); }

            if (json["source"].count("width")) {
              source.width = json["source"]["width"];
            } else { return abort("missing width for radial source mesh"); }

            if (json["source"].count("radius")) {
              source.radius = json["source"]["radius"];
            } else { return abort("missing radius for radial source mesh"); }

            if (json["source"].count("center")) {
              int i = 0;
              for (auto&& coord : json["source"]["center"]) { source.center[i++] = coord; }
            } else { return abort("missing center point for radial source mesh"); }

          } else {
            if (json["source"].count("min") and json["source"].count("max")) {
              int i = 0, j = 0;
              for (auto&& value : json["source"]["min"]) { source.min[i++] = value; }
              for (auto&& value : json["source"]["max"]) { source.max[j++] = value; }
            } else { return abort("missing source points extents"); }
          }

        } else if (json["source"].count("file")) {
          source.file = json["source"]["file"];
        } else { return abort("missing source points"); }

        if (json["source"].count("frame")) {
          source.field = json["source"]["field"];
          std::string const frame_type = json["source"]["frame"];
          if (frame_type == "local")       { source.frame = Frame::LOCAL; }
          else if (frame_type == "global") { source.frame = Frame::GLOBAL; }
          else { return abort("invalid source coords frame"); }
        } else { return abort("missing source coords frame"); }

      } else { return abort("missing source attributes"); }

      // target
      if (json.count("target")) {
        if (json["target"].count("size")) {
          int k = 0;
          for (auto&& count : json["target"]["size"]) { target.size[k++] = count; }

          if (json["target"].count("distrib")) {
            std::string const distrib = json["target"]["distrib"];
            if (distrib == "cartesian")   { target.distrib = CARTESIAN; }
            else if (distrib == "radial") { target.distrib = RADIAL; }
            else if (distrib == "random") { target.distrib = RANDOM; }
            else { return abort("invalid target points distribution"); }
          } else { return abort("missing target points distribution"); }

          if (target.distrib == RADIAL) {
            if (json["target"].count("span")) {
              target.span = json["target"]["span"];
            } else { return abort("missing span angle for radial target mesh"); }

            if (json["target"].count("width")) {
              target.width = json["target"]["width"];
            } else { return abort("missing width for radial target mesh"); }

            if (json["target"].count("radius")) {
              target.radius = json["target"]["radius"];
            } else { return abort("missing radius for radial target mesh"); }

            if (json["target"].count("center")) {
              int i = 0;
              for (auto&& coord : json["target"]["center"]) { target.center[i++] = coord; }
            } else { return abort("missing center point for radial target mesh"); }

          } else {
            if (json["target"].count("min") and json["target"].count("max")) {
              int i = 0, j = 0;
              for (auto&& value : json["target"]["min"]) { target.min[i++] = value; }
              for (auto&& value : json["target"]["max"]) { target.max[j++] = value; }
            } else { return abort("missing target points extents"); }
          }

        } else if (json["target"].count("file")) {
          target.file = json["target"]["file"];
        } else { return abort("missing target points"); }

        if (json["target"].count("frame")) {
          std::string const frame_type = json["target"]["frame"];
          if (frame_type == "local")       { target.frame = Frame::LOCAL; }
          else if (frame_type == "global") { target.frame = Frame::GLOBAL; }
          else { return abort("invalid target coords frame"); }
        } else { return abort("missing target coords frame"); }

        if (json["target"].count("output")) {
          target.output = json["target"]["output"];
        }

      } else { return abort("missing target attributes"); }

      // remap
      if (json.count("remap")) {
        if (json["remap"].count("basis")) {
          std::string const basis = json["remap"]["basis"];
          if (basis == "unitary")        { remap.basis = Portage::basis::Type::Unitary; }
          else if (basis == "linear")    { remap.basis = Portage::basis::Type::Linear; }
          else if (basis == "quadratic") { remap.basis = Portage::basis::Type::Quadratic; }
          else { return abort("invalid basis type"); }
        }

        if (json["remap"].count("weight")) {
          std::string const weight = json["remap"]["weight"];
          if (weight == "gather") { remap.weight = Portage::WeightCenter::Gather; }
          else if (weight == "scatter") { remap.weight = Portage::WeightCenter::Scatter; }
          else { return abort("invalid weights form"); }
        }

        if (json["remap"].count("kernel")) {
          std::string const kernel = json["remap"]["kernel"];
          if (kernel == "b4")                { remap.kernel = Portage::Weight::B4; }
          else if (kernel == "square")       { remap.kernel = Portage::Weight::SQUARE; }
          else if (kernel == "epanechnikov") { remap.kernel = Portage::Weight::EPANECHNIKOV; }
          else if (kernel == "polyramp")     { remap.kernel = Portage::Weight::POLYRAMP; }
          else if (kernel == "invsqrt")      { remap.kernel = Portage::Weight::INVSQRT; }
          else if (kernel == "coulomb")      { remap.kernel = Portage::Weight::COULOMB; }
          else if (kernel == "step")         { remap.kernel = Portage::Weight::STEP; }
          else { return abort("invalid weights kernel"); }
        }

        if (json["remap"].count("support")) {
          std::string const support = json["remap"]["support"];
          if (support == "elliptic")     { remap.support = Portage::Weight::ELLIPTIC; }
          else if (support == "tensor")  { remap.support = Portage::Weight::TENSOR; }
          else if (support == "faceted") { remap.support = Portage::Weight::FACETED; }
          else { return abort("invalid support geometry"); }
        }

        if (json["remap"].count("estimator")) {
          std::string const estimator = json["remap"]["estimator"];
          if (estimator == "kernel-density") { remap.estimator = Portage::KernelDensity; }
          else if (estimator == "lre")       { remap.estimator = Portage::LocalRegression; }
          else if (estimator == "operator")  { remap.estimator = Portage::OperatorRegression; }
          else { return abort("invalid estimator"); }
        }

      } else { return abort("missing remap attributes"); }

      // output
      if (json.count("output")) {

        if (json["output"].count("source")) {
          output.source = json["output"]["source"];
        } else { return abort("missing source field output"); }

        if (json["output"].count("exact")) {
          output.exact = json["output"]["exact"];
        } else { return abort("missing exact values output"); }

        if (json["output"].count("remap")) {
          output.remap = json["output"]["remap"];
        } else { return abort("missing remap values output"); }

        if (json["output"].count("error")) {
          output.error = json["output"]["error"];
        } else { return abort("missing error map output"); }

      } else { return abort("missing output attributes"); }

      // scale
      if (json.count("scale")) {

        if (json["scale"].count("coords")) {
          int i = 0;
          for (auto&& factor : json["scale"]["coords"]) { scale.coords[i++] = factor; }
        }

        if (json["scale"].count("field")) { scale.field = json["scale"]["field"]; }

        if (json["scale"].count("smooth")) {
          int i = 0;
          for (auto&& factor : json["scale"]["smooth"]) { scale.smooth[i++] = factor; }
        }
      }
    } catch(nlohmann::json::parse_error& e) { return abort(e.what()); }

    // everything was fine
    if (mpi.rank == 0) { print(json); }

    return true;
  }

  /**
   * @brief Handle runtime errors.
   *
   * @param message: a message to be displayed if any.
   * @param show_usage: hint for usage instructions print.
   * @return status
   */
  bool abort(std::string const& message = "", bool show_usage = true) const {
    if (mpi.rank == 0) {
      if (not message.empty()) {
        std::fprintf(stderr,
          " ---------------------------------------------------------- \n"
          " \e[31m Error: %s. \e[0m                                    \n"
          " ---------------------------------------------------------- \n",
          message.data()
        );
      }
      if (show_usage) { print_usage(); }
    }

    MPI_Finalize();
    return false;
  }

  /**
   * @brief
   *
   * @return
   */
  Wonton::Executor_type* get_executor() {
    return mpi.rank > 1 ? &mpi.executor : nullptr;
  }

private:
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

    std::printf("Parameters: \n");
    std::printf(" \u2022 MPI ranks: \e[32m%d\e[0m\n", mpi.num_ranks);

    if (json["source"].count("file")) {
      std::printf(" \u2022 source file: '\e[32m%s\e[0m'\n", base(source.file).data());
      std::printf(" \u2022 source field: '\e[32m%s\e[0m'\n", source.field.data());
    } else {
      std::string const frame = json["source"]["frame"];
      std::string const distrib = json["target"]["distrib"];
      std::printf(" \u2022 source size: '\e[32m[%d, %d]\e[0m'\n", source.size[0], source.size[1]);
      std::printf(" \u2022 source distrib: '\e[32m%s\e[0m'\n", distrib.data());

      if (source.distrib == RADIAL) {
        std::printf(" \u2022 span angle: '\e[32m%f\e[0m'\n", source.span);
        std::printf(" \u2022 width: '\e[32m%.2f\e[0m'\n", source.width);
        std::printf(" \u2022 radius: '\e[32m%.2f\e[0m'\n", source.radius);
        std::printf(" \u2022 center: '\e[32m[%.2f, %.2f]\e[0m'\n", source.center[0], source.center[1]);
      } else {
        std::printf(" \u2022 [x_min, y_min]: '\e[32m[%.2f, %.2f]\e[0m'\n", source.min[0], source.min[1]);
        std::printf(" \u2022 [x_max, y_max]: '\e[32m[%.2f, %.2f]\e[0m'\n", source.max[0], source.max[1]);
      }

      std::printf(" \u2022 source frame: '\e[32m%s\e[0m'\n", frame.data());
      std::printf(" \u2022 source field: '\e[32m%s\e[0m'\n", source.field.data());
    }

    if (json["target"].count("file")) {
      std::printf(" \u2022 target file: '\e[32m%s\e[0m'\n", base(target.file).data());
    } else {
      std::string const frame = json["target"]["frame"];
      std::string const distrib = json["target"]["distrib"];
      std::printf(" \u2022 target size: '\e[32m[%d, %d]\e[0m'\n", target.size[0], target.size[1]);
      std::printf(" \u2022 target distrib: '\e[32m%s\e[0m'\n", distrib.data());

      if (target.distrib == RADIAL) {
        std::printf(" \u2022 span angle: '\e[32m%.2f\e[0m'\n", target.span);
        std::printf(" \u2022 width: '\e[32m%.2f\e[0m'\n", target.width);
        std::printf(" \u2022 radius: '\e[32m%.2f\e[0m'\n", target.radius);
        std::printf(" \u2022 center: '\e[32m[%.2f, %.2f]\e[0m'\n", target.center[0], target.center[1]);
      } else {
        std::printf(" \u2022 [x_min, y_min]: '\e[32m[%.2f, %.2f]\e[0m'\n", target.min[0], target.min[1]);
        std::printf(" \u2022 [x_max, y_max]: '\e[32m[%.2f, %.2f]\e[0m'\n", target.max[0], target.max[1]);
      }

      std::printf(" \u2022 target frame: '\e[32m%s\e[0m'\n", frame.data());
    }

    if (json.count("scale")) {
      std::printf(" \u2022 coords scale: '\e[32m[%.e, %.e]\e[0m'\n", scale.coords[0], scale.coords[1]);
      std::printf(" \u2022 field scale: '\e[32m%.e\e[0m'\n", scale.field);
      std::printf(" \u2022 smooth scale: '\e[32m[%.2f, %.2f]\e[0m'\n", scale.smooth[0], scale.smooth[1]);
    }

    std::string const basis = json["remap"]["basis"];
    std::string const weight = json["remap"]["weight"];
    std::string const kernel = json["remap"]["kernel"];
    std::string const support = json["remap"]["support"];
    std::string const estimator = json["remap"]["estimator"];

    std::printf(" \u2022 weight basis: '\e[32m%s\e[0m'\n", basis.data());
    std::printf(" \u2022 weight center: '\e[32m%s\e[0m'\n", weight.data());
    std::printf(" \u2022 weight kernel: '\e[32m%s\e[0m'\n", kernel.data());
    std::printf(" \u2022 weight support: '\e[32m%s\e[0m'\n", support.data());
    std::printf(" \u2022 estimator: '\e[32m%s\e[0m'\n", estimator.data());

    std::printf(" \u2022 exact output: '\e[32m%s\e[0m'\n", output.exact.data());
    std::printf(" \u2022 remap output: '\e[32m%s\e[0m'\n", output.remap.data());
    std::printf(" \u2022 error output: '\e[32m%s\e[0m'\n", output.error.data());
  }

  /**
   * @brief Display run instructions and input file format.
   *
   */
  static void print_usage() {
    std::printf(
      " Usage: ./particle-analysis params.json                             \n\n"
      " params.json:                                                         \n"
      " {                                                                    \n"
      "  \e[32m'source'\e[0m: {                                              \n"
      "    \e[32m'file'\e[0m: '/path/to/wavelets.csv',                       \n"
      "    \e[32m'field'\e[0m: <math>                                        \n"
      "    \e[32m'frame'\e[0m: <global|local>,                               \n"
      "  },                                                                  \n"
      "  \e[32m'target'\e[0m: {                                              \n"
      "    \e[32m'size'\e[0m: [<int>, <int>],                                \n"
      "    \e[32m'distrib'\e[0m: <radial|cartesian|random>,                  \n"
      "    \e[32m'span'\e[0m: <double>,                             # radial \n"
      "    \e[32m'width'\e[0m: <double>,                            # radial \n"
      "    \e[32m'radius'\e[0m: <double>,                           # radial \n"
      "    \e[32m'center'\e[0m: [<double>, <double>],               # radial \n"
      "    \e[32m'min'\e[0m: [<double>, <double>],        # cartesian|random \n"
      "    \e[32m'max'\e[0m: [<double>, <double>],        # cartesian|random \n"
      "    \e[32m'frame'\e[0m: <global|local>                                \n"
      "    \e[32m'output'\e[0m: <string>,                                    \n"
      "  },                                                                  \n"
      "  \e[32m'scale'\e[0m: {                                               \n"
      "    \e[32m'coords'\e[0m: [<double>, <double>],                        \n"
      "    \e[32m'field'\e[0m: <double>,                                     \n"
      "    \e[32m'smooth'\e[0m: [<double>, <double>]                         \n"
      "  },                                                                  \n"
      "  \e[32m'remap'\e[0m: {                                               \n"
      "    \e[32m'basis'\e[0m: <unitary|linear|quadratic>                    \n"
      "    \e[32m'weight'\e[0m: <gather|scatter>                             \n"
      "    \e[32m'kernel'\e[0m: <b4|square|polyramp|invsqrt|coulomb|step>    \n"
      "    \e[32m'support'\e[0m: <elliptic|tensor|faceted>                   \n"
      "    \e[32m'estimator'\e[0m: <lre|kernel-density|operator>             \n"
      "  },                                                                  \n"
      "  \e[32m'output'\e[0m: {                                              \n"
      "    \e[32m'source'\e[0m: <string>,                                    \n"
      "    \e[32m'exact'\e[0m: <string>,                                     \n"
      "    \e[32m'remap'\e[0m: <string>,                                     \n"
      "    \e[32m'error'\e[0m: <string>                                      \n"
      "  }                                                                   \n"
      "}                                                                     \n"
    );
  }
};