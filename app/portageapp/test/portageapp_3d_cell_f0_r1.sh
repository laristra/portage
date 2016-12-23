#!/bin/bash
: <<'END'
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
END

# Exit on error
set -e
# Echo each command
set -x

# 3D, 1st order accurate remap of constant cell-centered function

# SERIAL RUN

$APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=y \
--entity_kind=cell --field_order=0 \
--remap_order=1 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f0_r1.txt field_3d_cell_f0_r1.txt 1e-12

# PARALLEL RUN

mpirun -np 4 $APPDIR/portageapp_jali \
--dim=3 --nsourcecells=5 --ntargetcells=7 \
--conformal=n \
--entity_kind=cell --field_order=0 \
--remap_order=1 \
--output_results=y

# Compare the values for the field
$APPDIR/apptest_cmp GOLD-field_3d_cell_f0_r1_nc.txt.0 field_3d_cell_f0_r1_nc.txt.0 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f0_r1_nc.txt.1 field_3d_cell_f0_r1_nc.txt.1 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f0_r1_nc.txt.2 field_3d_cell_f0_r1_nc.txt.2 1e-12
$APPDIR/apptest_cmp GOLD-field_3d_cell_f0_r1_nc.txt.3 field_3d_cell_f0_r1_nc.txt.3 1e-12
