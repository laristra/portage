

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
