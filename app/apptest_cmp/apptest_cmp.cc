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




#include <sys/time.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <cmath>
#include <iostream>
#include <stdexcept>

void print_usage() {
  std::printf("Usage: apptest_cmp file_gold file abs_eps [rel_eps=abs_eps]\n");
}

void load_field(std::iostream &s, std::vector<int> &gid,
                std::vector<double> &values) {
  int g;
  double v;
  while (s >> g >> v) {
    gid.push_back(g);
    values.push_back(v);
  }
}

int main(int argc, char** argv) {
  if (argc < 4 || argc > 5) {
    print_usage();
    return 1;
  }
  double abs_eps = std::stod(argv[3]);
  double rel_eps = (argc == 5) ? std::stod(argv[4]) : abs_eps;

  std::fstream f1(argv[1]), f2(argv[2]);
  if (!f1) throw std::runtime_error("First file cannot be opened.");
  if (!f2) throw std::runtime_error("Second file cannot be opened.");

  std::vector<int> gid1, gid2;
  std::vector<double> values1, values2;
  load_field(f1, gid1, values1);
  load_field(f2, gid2, values2);

  std::cout << std::scientific;
  std::cout.precision(17);
  std::cout << "Comparing files: " << argv[1] << " " << argv[2] << std::endl;
  std::cout << "Absolute Epsilon: " << abs_eps << std::endl;
  std::cout << "Relative Epsilon: " << rel_eps << std::endl;
  std::cout << "Field sizes: " << gid1.size() << " " << gid2.size() <<
    std::endl;
  if (gid1.size() != gid2.size()) {
    throw std::runtime_error("The field sizes do not match.");
  }
  double delta, maxval;
  for (int i=0; i < gid1.size(); i++) {
    if (gid1[i] != gid2[i]) {
      std::cout << i << " " << gid1[i] << " " << gid2[i] << std::endl;
      throw std::runtime_error("The field global IDs do not match.");
    }
    // Test absolute tolerance
    delta = std::abs(values1[i] - values2[i]);
    if (delta <= abs_eps)
      continue;
    // Still might be the same if within the relative tolerance
    maxval = (std::abs(values1[i]) > std::abs(values2[i])) ?
      values1[i] : values2[i];
    if ((delta / std::abs(maxval)) <= rel_eps)
      continue;
    // Fails -- abort
    std::cout << i << " " << values1[i] << " " << values2[i]
              << " Abs Err: " << delta
              << " Rel Err: " << delta/std::abs(maxval) << std::endl;
    throw std::runtime_error("The field values do not match.");
  }

  return 0;
}
