



#ifndef MPI_COLLATE_H_
#define MPI_COLLATE_H_

#include <algorithm>
#include <numeric>
#include <vector>

#include <mpi.h>

namespace Portage {

  // Collates local vectors lvec (of various lengths) into a global vector gvec
  // on rank 0. The `gvec` vector is only defined on rank 0, and it is resized
  // there to the proper size. On other ranks it is resized to size 1 and left
  // undefined.
  template <typename T>
  void collate_type(MPI_Comm comm, const int rank, const int numpe,
               const MPI_Datatype mpi_type,
               std::vector<T> &lvec, std::vector<T> &gvec) {
    std::vector<int> lvec_sizes(numpe);
    int lvec_size = lvec.size();
    MPI_Gather(&lvec_size, 1, MPI_INT, &lvec_sizes[0], 1, MPI_INT, 0, comm);
    std::vector<int> displs;
    if (rank == 0) {
      int gvec_size = std::accumulate(lvec_sizes.begin(), lvec_sizes.end(),
                                      0);
      gvec.resize(gvec_size);
      displs.resize(lvec_sizes.size());
      int idx = 0;
      for (int i=0; i < lvec_sizes.size(); i++) {
        displs[i] = idx;
        idx += lvec_sizes[i];
      }
    } else {
      // We resize to size 1, so that the expressions &gvec[0], &displs[0] below
      // are well defined on all ranks.
      gvec.resize(1);
      displs.resize(1);
    }
    MPI_Gatherv(&lvec[0], lvec.size(), mpi_type, &gvec[0], &lvec_sizes[0],
        &displs[0], mpi_type, 0, comm);
  }

  void collate(MPI_Comm comm, const int rank, const int numpe,
               std::vector<int> &lvec, std::vector<int> &gvec) {
    collate_type(comm, rank, numpe, MPI_INT, lvec, gvec);
  }

  void collate(MPI_Comm comm, const int rank, const int numpe,
               std::vector<double> &lvec, std::vector<double> &gvec) {
    collate_type(comm, rank, numpe, MPI_DOUBLE, lvec, gvec);
  }

  // Returns the indices that would sort an array.
  template <typename T>
  void argsort(const std::vector<T> &x, std::vector<int> &idx) {
    idx.resize(x.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&x](int a, int b){ return x[a] < x[b]; });
  }

  // Reorders a vector x using x = x[idx], where idx is a vector of indices
  template <typename T>
  void reorder(std::vector<T> &x, const std::vector<int> &idx) {
    std::vector<T> y(x.size());
    for (int i=0; i < x.size(); i++) y[i] = x[idx[i]];
    x = y;
  }

} // namespace Portage

#endif // MPI_COLLATE_H_
