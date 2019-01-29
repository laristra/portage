# Mismatch Fixup  {#mismatch}

The Portage library includes a repair algorithm to handle remap errors
arising from the mismatch of the boundaries of two meshes.
happen because curved boundaries are discretized with different
resolutions in the source and target meshes or because two physics
packages exchanging data through Portage view the geometry of the
domain a bit differently. If all of the source mesh is not covered by
the target mesh or vice versa, the result may violate conservation or
introduce artifacts in the fields.

The four possible cases of mismatch are:
    * Target cell is partially covered by source cells (partially filled).
    * Target cell does not overlap any source cells (empty).
    * Source cell is partially covered by target cells (partially contributing).
    * Source cell does not overlap any target cells (non-contributing)

If a source cell is partially covered or not covered
at all by the result is a conservation error (missing material in the
target mesh). If a target cell is partially covered by source cells, the amount of
material in the cell or its resulting field value may be incorrect
depending on the weighting used in the interpolation step (more on
that below). If a target cell is empty, then it's field value
is undefined or some value that all cells are initialized to, which is
clearly an error. 

Consider a target cell, \\( \\Omega^t_i \\) that
intersects a set of source cells \\( \\{ \\Omega^s_j \\}_{j=1,n} \\)
and is fully covered by the source mesh. Assume that the field value at each
of the source cells is \\( u^s_j \\). The sum of the intersection
volumes is \\( \\sum_j |\\Omega^t_i \\cap
\\Omega^s_j| \\) and assuming a first order accurate
interpolation, the accumulated contribution of the overlapping source
cells to the target cell is \\( U^t_j = \\sum_j u^s_j|\\Omega^t_i \\cap
\\Omega^s_j| \\). We derive a cell-averaged value \\( u^t_j
\\) from the integral quantity \\( U^t_j \\) as
\\[ u^t_i = \\frac{U^t_j}{|\\Omega^t_i|} = \\frac{U^t_j}{\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|} \\]

It is guaranteed that this form of interpolation is **locally conservative**,
i.e., any material from overlapping
parts of the source mesh are accounted for in the target cell. If, in
addition, \\[ u^s_j = c \\quad \\forall j \\quad \\mbox{then} \\quad
u^t_i = \\frac{U^t_j}{\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|} = \\frac{\\sum_j u^s_j|\\Omega^t_i \\cap
\\Omega^s_j|}{\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|} = \\frac{c\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|}{\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|} = c \\]
i.e. it is **constant preserving**.

## The issue with partially covered target cells

When \\( \\Omega^t_i \\) is partially covered i.e. \\( |\\Omega^t_i|
\\neq \\sum_j |\\Omega^t_i \\cap \\Omega^s_j| \\), how we
convert the integral value \\( U^t_j \\) to a cell-averaged quantity
\\( u^t_j \\) has significant implications. If we say, \\[ u^t_j =
\\frac{U^t_j}{|\\Omega^t_i|} \\] then we ensure that the remap onto the
target cell is **locally conservative**. However, since \\(
\\sum_j |\\Omega^t_i \\cap \\Omega^s_j|  \\lt  |\\Omega^t_i|
\\), this will cause the cell averaged
value to be lower than it would be if it were fully covered. If, as
before, \\[ u^s_j = c \\quad \\forall j \\qquad \\mbox{then} \\qquad
u^t_i = \\frac{U^t_j}{|\\Omega^t_i|} = \\frac{\\sum_j u^s_j|\\Omega^t_i \\cap
\\Omega^s_j|}{|\\Omega^t_i|} = \\frac{c\\sum_j |\\Omega^t_i \\cap
\\Omega^s_j|}{|\\Omega^t_i|} = c_1 < c
\\]
i.e. it is **NOT constant-preserving**.

On the the other hand, if we compute the cell-averaged quantity as
\\[
u^t_i = \\frac{U^t_i}{\\sum_j |\\Omega^t_i \\cap \\Omega^s_j|}
\\]
the result is **constant-preserving** but **NOT locally conservative**
since \\[ u^t_i|\\Omega^t_i|
\\neq \\sum_j u^s_j|\\Omega^t_i \\cap
\\Omega^s_j| \\]

The supplied interpolation routines in Portage adopt the second
method, i.e. they are constant-preserving but not
conservative. Conservation is restored, if requested, in the repair
algorithm described below.

## Discrepancy computation

Given a target mesh \\( \\{\\Omega^s_j\\} \\) and a source mesh \\(
\\{\\Omega^s_j\\} \\), we compute the following (over all partitions
in a distributed mesh):

\\[ \mbox{Target} \; \mbox{Volume} = V^t = \\sum_j|\\Omega^t_j| \\]
\\[ \mbox{Source} \; \mbox{Volume} = V^s = \\sum_i|\\Omega^s_i| \\]
\\[ \mbox{Intersection} \; \mbox{Volume} = V^i = \\sum_j|\\Omega^t_j \\cap \\{\\Omega^s_i\\}| \\]

If \\( V^t \\neq V^i \\), we can conclude that one or more cells of
the target mesh are not covered by source mesh. If \\( V^s \\neq V^i
\\), one or more cells of the source mesh are not
covered by target mesh. It is possible that both conditions exist in a
pair of meshes.

The conservation error for a field whose values on the target mesh are
\\( \\{u^s_i\\} \\) and on the source mesh are \\( \\{u^t_i\\} \\), is
given by
\\[
\\Delta{U} = \\sum_i{u^s_i|\\Omega^s_i|} - \\sum_j{u^t_j|\\Omega^t_j|}
\\]

## Repair

The repair function in Portage performs a global repair. While its
exact details depends on user choices for what happens in empty cells
and partially filled cells, the algorithm has these main steps:

### Empty layers detection

The algorithms labels empty cells starting from
completely filled or partially filled cells and fanning out. Completely or
partially filled are labeled as belonging to layer 1; subsequent
layers are built up of empty cells going outward. Empty cells that
border at least one cell that belongs to layer *N*  are labeled as
belonging to layer *N+1*. The process continues until no more empty cells
are available for placing into layers.

### Filling empty cells

The calling application may choose to leave the empty cells
untouched (Portage::LEAVE_EMPTY) or to fill it using
extrapolation (Portage::EXTRAPOLATE).

If the application asks for empty cells to be filled by extrapolation,
the algorithm populates them layer by layer in ascending order. The
cell-centered value of any empty cell is considered to be the average
of all populated neighbors in the previous layer. If the field is a
constant, this will preserve the constant.

### Restoring conservation

The methods used for populating partially filled cells in the
interpolation step and for populating empty cells in the previous step
ensure that perturbations in the field are minimized (constant fields
are preserved) but they may cause a violation of conservation. In this
step, the repair algorithm offers the opportunity of violating
conservation but leaving the fields as computed in the interpolate
step (Portage::CONSTANT), restoring conservation
in partially filled cells (Portage::CONSERVATIVE,
by scaling the source contribution by the reciprocal of \\(
\\sum_j|\\Omega^t_j \\cap \\{\\Omega^s_i\\}| \\) instead of the
reciprocal of \\( |\\Omega^t_j| \\) ) or by restoring conservation
through a global distribution of the excess
(Portage::SHIFTED_CONSERVATIVE, described below).

Given a conservation error of \\( \\Delta{U} \\), the algorithm for
the Portage::SHIFTED_CONSERVATIVE option attempts to add (or subtract) a
fraction of the error from each cell. The fraction of the error is
proportional to fractional volume of the cell with respect to the
total mesh volume.  \\[ \\delta{U}_i =
\\Delta{U}\\frac{|\\Omega^t_i|}{\\sum_j|\\Omega^t_j|} \\]

If the calling application has specified global bounds for the field,
then the amount that is added (or subtracted) is capped so as to not
violate these bounds; then the global deficit and the deficit for
subsequent cells is recalculated. The algorithm may require multiple
iterations for complete redistribution due to bounds
preservation. In the example driver (Portage::MMDriver), the bounds
are computed from the source field, with a little slack for constant
fields, if the global bounds are not specified. 

The effect of restoring conservation on the original mesh using the
above algorithm is that field maintains its character to the extent
possible but its actual values are shifted. Thus, a constant field
might become a different constant field on the target mesh.


