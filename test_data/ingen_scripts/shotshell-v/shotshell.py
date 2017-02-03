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
from ingen.materials import material

innerCaseRadius = 1.75
slugDX = 0.05

def suBlock(f,mat):
  # HACK: names is internal
  uf=altair.uFrame([getattr(f,a) for a in f.names],mat)
  # FIXME: without factor, this crashed with close generators
  uf.dxStipple(altair.voronoiMesh,f.dx()*1.1,100)
  return uf.block()

def block3(iMin,jMin,jMax,rule=altair.copyDistrib(),material=None,feather=None):
  # HACK: only copyDistrib() supported
  f=altair.frame3(iMin,jMin,jMax,material)
  f.copy()
  return suBlock(f,material),f.iMax

def block2(jMin,jMax,rule,material,feather=None):
  # HACK: only squareDistrib() supported
  f=altair.frame2(jMin,jMax,material)
  f.equal(f.dx())
  return suBlock(f,material)

# For the original structured mesh:
# from ingen.altair import block2,block3

ingen.init(globals(),None)
gwiz.loadDirectory('contours',[('*',gwiz.noRule(),gwiz.dxDistribRule(0.005))],cntr)
mat.pb = material(1)
mat.plastic = material(2)
mat.void = material(3)

seg.slug_fin_jMin = altair.segment(cntr.slug,cntr.slug.pnt.p14,cntr.slug.pnt.p13).equalArcDistrib(slugDX)
seg.slug_fin_jMax = altair.segment(cntr.slug,cntr.slug.pnt.p12,cntr.slug.pnt.p11).equalArcDistrib(slugDX)
seg.slug_fin_iMin = altair.segment(cntr.slug,cntr.slug.pnt.p14,cntr.slug.pnt.p12).equalArcDistrib(slugDX)
blk.slug_fin,sfiMax = block3(seg.slug_fin_iMin,seg.slug_fin_jMin,seg.slug_fin_jMax,material=mat.pb)

seg.slug_hp = altair.segment(cntr.slug,cntr.slug.pnt.p7,cntr.slug.pnt.p8).equalArcDistrib(slugDX)
seg.slug_nose = altair.segment(cntr.slug,cntr.slug.pnt.p8,cntr.slug.pnt.p9).equalArcDistrib(slugDX)
seg.slug_base = altair.segment(cntr.slug,cntr.slug.pnt.p6,cntr.slug.pnt.p13).equalArcDistrib(slugDX)

seg.slug_body_fore = altair.segment(cntr.slug,cntr.slug.pnt.p9,cntr.slug.pnt.p10).equalArcDistrib(slugDX)
seg.slug_body_aft = altair.segment(cntr.slug,cntr.slug.pnt.p10,cntr.slug.pnt.p11).equalArcDistrib(slugDX)
sbiMin=seg.slug_body_aft+seg.slug_body_fore
blk.slug_body,_ = block3(sbiMin,
                              seg.slug_base+sfiMax,
                              seg.slug_hp+seg.slug_nose,
                              material=mat.pb,
                              rule=altair.copyDistrib(),
                              feather=altair.fthr())

def avgDX(*ss): return sum(s.dx() for s in ss)/len(ss)
cntr.sabot_outer = gwiz.rLine(innerCaseRadius,cntr.slug.pnt.p12.z+0.25,cntr.slug.pnt.p9.z+0.25)
seg.sabotJMin=sbiMin+seg.slug_fin_jMax
seg.sabotJMin.slide('ss_slide')
seg.sabot_outer = altair.segment(cntr.sabot_outer)\
    .equalArcDistrib(len(seg.sabotJMin))
seg.sabotIMax=altair.segment(None,cntr.slug.pnt.p9,cntr.sabot_outer[-1])\
    .equalArcDistrib(avgDX(seg.sabot_outer,seg.sabotJMin))

cntr.wad_aft = gwiz.zLine(3.5,innerCaseRadius)
cntr.wad_fore = gwiz.zLine(cntr.slug.pnt.p14.z,cntr.slug.pnt.p14.r)
seg.wad_fore = altair.segment(cntr.wad_fore).equalArcDistrib(slugDX)
seg.wadJMax=seg.wad_fore+seg.slug_fin_iMin
seg.wadJMax.slide('wf_slide')
seg.wad_aft = altair.segment(cntr.wad_aft).equalArcDistrib(seg.wadJMax.dx())
seg.slug_base.slide('wf_slide')
wadDX=avgDX(seg.wad_aft,seg.wadJMax)
seg.wadIMax=altair.segment(None,cntr.wad_aft[-1],
                           cntr.sabot_outer[0]).equalArcDistrib(wadDX)
seg.wadIMin=altair.segment(None,cntr.wad_aft[0],
                           cntr.wad_fore[0]).equalArcDistrib(wadDX)

cntr.cush_aft = gwiz.zLine(2.75,innerCaseRadius)
seg.cush_aft = altair.segment(cntr.cush_aft).equalArcDistrib(seg.wad_aft)
blk.cushion = block2(seg.cush_aft,seg.wad_aft,material=mat.void,rule=altair.squareDistrib())
seg.wad_aft.slide('cw-slide')

def interp(a,b,t): return a*(1-t)+b*t
genWadR,genWadZ=int(20/slugDX),int(80/slugDX)
corner=cntr.slug.pnt.p12
gens=[gwiz.rz((i+.5)*innerCaseRadius/genWadR,(3.5+corner.z)/2)
      for i in xrange(genWadR)]+\
    [gwiz.rz((corner.r+innerCaseRadius)/2,
             interp(corner.z,cntr.slug.pnt.p9.z,(i+.5)/genWadZ))
     for i in xrange(genWadZ)]
#for p in gens: pnt().add(p,"gen")
blk.wad_sabot=altair.uFrame([seg.wad_fore+seg.slug_fin_iMin,seg.sabotJMin,
                             seg.sabotIMax,seg.sabot_outer+seg.wadIMax,
                             seg.wad_aft,seg.wadIMin],mat.plastic)\
                             .rule(altair.voronoiMesh(gens,1000)).block()

altair.finalize()
altair.writeGMV('mesh/shotshell-v.gmv')
altair.writeX3D('mesh/shotshell-v.x3d')
