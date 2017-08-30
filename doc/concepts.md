# portage Concepts      {#concepts}

The remapping algorithm within portage is divided into three phases,
which are roughly labelled as

1. _search_ - find candidate cells/particles that will contribute to
   remap of a given target cell/particle
2. _intersect_ - calculate the weight of each candidate's contribution
   to the remap of a given target cell/particle; this may include
   higher-order moments if requested.
3. _interpolate_ - using the weights and moments, along with
   appropriate limiters, reconstruct the field data for a given target
   cell/particle.

All of these operations work with your underlying mesh/particles and
state manager through wrappers that provide an interface the queries
needed to perform any particular step.  For an example of the
requirements of the wrappers, see the [Example Use](@ref example)
paeg.  In the below, when we refer to _mesh_ or _particles_ in terms
of the operations, we really mean _mesh wrappers_ and _particle swarm
wrappers_.

All operations consist of a _source_ mesh/particle swarm and a
_target_ mesh/particle swarm.  The _source_ entity is the one where we
have existing field data, and the _target_ entity is object to which
we would like to remap the data.

## Search

Given a source and target entity, this step simply _identifies_ which
parts of the source contribute to which parts of the target.
Concretely, for mesh-mesh remap this step would identify which source
cells overlap each target cell. **ADD IMAGE**

portage has several search algorithms with varying degrees of
sophistication/speed.  **ADD LINKS TO API DOCS**

- SearchSimple - 2d, bounding box search
- SearchKDTree - 2d or 3d parallel k-d tree search
- SearchSimplePoints - any-d quadratic time search over particle swarms
- SearchPointsByCells - any-d linear time search over particle swarms

## Intersect

Given the list of source candidates for intersection for a target
entity, this step calculates the actual weights going into the
intersection.

For meshes, this step uses exact intersection methods to calculate
various _moments_ of the polygon/polyhedron of intersection; moments
of higher order than the 0th (i.e. the area or volume of the
intersection) are needed for higher order remap.  It is possible that
some candidates are determined to have zero intersection, or that some
candidates have multiple intersections in the case of non-convex
cells.  The available intersectors for meshes are: **LINKS**

- IntersectClipper - 2d, exact intersection method based on
  the [Clipper](www.angusj.com/delphi/clipper.php) library for polygon
  intersection and clipping
- IntersectR2D - 2d, fast, exact polygonal intersection method based
  on the [r3d](https://github.com/laristra/r3d) library.
- IntersectR3D - 3d, fast, exact polyhedral intersection method base
  don the [r3d](https://github.com/laristra/r3d) library.

For particles, this step is refferred to as _accumulation_.  The
distinction in terminology stems from the fact that for particles,
local regression estimators (LRE) are used to do the remap.  In this
stage, the LRE weights from particle contributions are accumulated by
computing the weight functions and local regression corrections to
those weights.

## Interpolate

