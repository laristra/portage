'''
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
'''

import ingen
from ingen import gwiz, csv, altair
from ingen.gwiz import solid
from ingen.materials import material

innerCaseRadius = 1.75
slugDX = 0.04

# Distribute N points as evenly as possible along a list of segments
def multiSegNDistrib(npts,segList):
    tarc = sum([s.arclen() for s in segList])
    counts = [int(round(npts*s.arclen()/tarc)) for s in segList]
    reqCount = npts+len(segList)-1
    if sum(counts) > reqCount: counts[-1] = counts[-1]-1
    elif sum(counts) < reqCount: counts[-1] = counts[-1]+1
    for s,c in zip(segList,counts): s.equalArcDistrib(npts=c)


ingen.startModel('mesh/shotshell',convenience=globals())
gwiz.loadDirectory('contours',[('*',gwiz.noRule(),gwiz.dxDistribRule(0.005))],cntr)
mat.pb = material(1)
mat.plastic = material(2)
mat.void = material(3)

seg.slug_fin_jMin = altair.segment(cntr.slug,cntr.slug.pnt.p14,cntr.slug.pnt.p13).equalArcDistrib(slugDX)
seg.slug_fin_jMax = altair.segment(cntr.slug,cntr.slug.pnt.p12,cntr.slug.pnt.p11).equalArcDistrib(slugDX)
seg.slug_fin_iMin = altair.segment(cntr.slug,cntr.slug.pnt.p14,cntr.slug.pnt.p12).equalArcDistrib(slugDX)
blk.slug_fin = altair.block3(seg.slug_fin_iMin,seg.slug_fin_jMin,seg.slug_fin_jMax,material=mat.pb)

seg.slug_hp = altair.segment(cntr.slug,cntr.slug.pnt.p7,cntr.slug.pnt.p8).equalArcDistrib(slugDX)
seg.slug_nose = altair.segment(cntr.slug,cntr.slug.pnt.p8,cntr.slug.pnt.p9).equalArcDistrib(slugDX)
seg.slug_base = altair.segment(cntr.slug,cntr.slug.pnt.p6,cntr.slug.pnt.p13).equalArcDistrib(slugDX)

seg.slug_body_fore = altair.segment(cntr.slug,cntr.slug.pnt.p9,cntr.slug.pnt.p10).equalArcDistrib(slugDX)
seg.slug_body_aft = altair.segment(cntr.slug,cntr.slug.pnt.p10,cntr.slug.pnt.p11).equalArcDistrib(slugDX)
blk.slug_body = altair.block3(seg.slug_body_aft+seg.slug_body_fore,
                              seg.slug_base+blk.slug_fin.iMax(),
                              seg.slug_hp+seg.slug_nose,
                              material=mat.pb,
                              rule=altair.copyDistrib(),
                              feather=altair.fthr())

cntr.sabot_outer = gwiz.rLine(innerCaseRadius,cntr.slug.pnt.p12.z+0.25,cntr.slug.pnt.p9.z+0.25)
seg.sabot_outer = altair.segment(cntr.sabot_outer).equalArcDistrib(blk.slug_body.iMin()+blk.slug_fin.jMax())
blk.sabot = altair.block2(blk.slug_body.iMin()+blk.slug_fin.jMax(),
                          seg.sabot_outer,rule=altair.squareDistrib(),
                          material=mat.plastic,
                          feather=altair.fthr())
blk.sabot.jMin().slide('ss_slide')

cntr.wad_aft = gwiz.zLine(3.5,innerCaseRadius)
cntr.wad_fore = gwiz.zLine(cntr.slug.pnt.p14.z,cntr.slug.pnt.p14.r)
seg.wad_fore = altair.segment(cntr.wad_fore).equalArcDistrib(slugDX)
seg.wad_aft = altair.segment(cntr.wad_aft).equalArcDistrib(seg.wad_fore+blk.slug_fin.iMin()+blk.sabot.iMax())
blk.wad = altair.block2(seg.wad_aft,
                        seg.wad_fore+blk.slug_fin.iMin()+blk.sabot.iMax(),
                        material=mat.plastic,
                        rule=altair.squareDistrib(),
                        feather=altair.fthr())
blk.wad.jMax().slide('wf_slide')
seg.slug_fin_jMin.slide('wf_slide')
seg.slug_base.slide('wf_slide')

cntr.cush_aft = gwiz.zLine(2.75,innerCaseRadius)
seg.cush_aft = altair.segment(cntr.cush_aft).equalArcDistrib(seg.wad_aft)
blk.cushion = altair.block2(seg.cush_aft,seg.wad_aft,material=mat.void,rule=altair.squareDistrib())
seg.wad_aft.slide('cw-slide')


ingen.endModel(gmv=True,x3d=True)
altair.blocks2Regions(blk,reg)
solid.write(reg,'shotshell')

