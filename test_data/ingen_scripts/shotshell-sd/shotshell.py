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

