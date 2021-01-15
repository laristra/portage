/*
This file is part of the Ristra Portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

// Dummy  app to  check if  we end  up with  multiply defined  symbols
// because of function definitions included in header files (class
// variables are always inlined)

#include <iostream>
#include "portage_includes_to_check.h"
#include "portage_includes_to_check.h"  // check for proper include guards

int main(int argc, char **argv) {
  std::cout << "Test app to see if we get multiply defined symbols\n";
  return 0;
}
