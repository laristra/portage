/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



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
