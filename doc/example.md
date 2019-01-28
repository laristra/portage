# Simple Mesh Example    {#example}

Portage provides very crude mesh and state manager frameworks aptly
called `Simple_Mesh` and `Simple_State` through its support package, 
[Wonton](https://github.com/laristra/wonton).  The goal of these frameworks
is to show how one can wrap their favorite mesh and state manager for
use with portage - _they should not be used in any production sense._
The details on the `Wonton::Simple_Mesh` and `Wonton::Simple_State` 
frameworks, and its wrappers `Wonton::Simple_Mesh_Wrapper` and 
`Wonton::Simple_State_Wrapper`, can be found in the
[Wonton](https://github.com/laristra/wonton) documentation.

# Wrappers

As mentioned on the [Concepts](@ref concepts) page, the search,
intersect, and interpolate actions in portage all operate on
mesh/particle and state _wrappers_.  The reason for this is that these
three steps may need to ask for some information that may not be
readily available within the original mesh/particle and state
frameworks, but could be constructed within a wrapper.

Portage's support package, [Wonton](https://github.com/laristra/wonton), 
provides a helper class, `Wonton::AuxMeshTopology`, that assists in
extending a basic mesh's topological entities needed for some remap
capabilities.  One does not _need_ to use the AuxMeshTopology class,
especially if one's mesh already efficiently supports the advanced
mesh topologies. 
In particular, the advanced mesh topologies and entities, which we
term _sides_, _corners_, and _wedges_ are utilized in node-centered
remapping, but are also _required_ if there are cells with non-planar
faces in the remap.  In 2d, a side is triangle composed of the cell
center and the two nodes of an edge; a wedge is a triangle composed of
half of a side, by the cell center and the edge's midpoint and one of
its nodes; a corner is a quadrilateral formed by the two wedges in a
cell attached to a node.  In 3d, the idea is extended such that a side
is a tetrahedron composed of the cell center, the two nodes of an
edge, and the face center of a face attached to that edge; wedges are
again half of a side, but this time a tetrahedron composed of the cell
center, an edge's midpoint and one of its nodes, and the face center;
a corner is the collection of all wedges in a cell attached to a node.

# Applications and Tests

Nearly all of the generic unit tests have been designed to work with
the `Wonton::Simple_Mesh_Wrapper` and `Wonton::Simple_State_Wrapper`
data structures when a full up mesh is required.  These should serve
as examples of the types of things one may want to do with a mesh and
state framework.

## Simple_mesh_app

Better examples of how to use these wrappers to actually do a remap of
field data are in the application tests.  In particular, the
`app/simple_mesh_app/simple_mesh_app.cc` program shows how to wrap
mesh and state objects, adds some field data to the source state, and
utilize `Portage::Driver` with various search, intersect, and
interpolate algorithms to perform the remap.  The `Portage::Driver` is
templated on mesh wrapper type and state wrapper type, and can be used
with other frameworks.  It need not be used at all, but is a nice
starting point for writing one's own remap application.

## Portageapp_jali

This application demonstrates how to remap (scalar) field data using Jali mesh 
infrastructure. The application `app/portageapp/portageapp_jali.cc` is capable
of remapping cell- and nodal fields defined on square/cubic meshes or on  
unstructured meshes in 2D and 3D, respectively. 
![Visualization of a source mesh (red) and a part of a target mesh (black) with remapped data using Paraview.](doxygen/images/jaliapp_example.png)

The app allow users to specify an analytic field on a source mesh. This is done 
by an algebraic expression in the standard input, see options and examples bellow.   
An error or remap is calculated by comparing remapped values to exact values in the 
cell-centroids/nodes on a target mesh and the error is printed to standard output. 
The remapped field on a target mesh is saved into a file `output.exo` 
if it is not specified otherwise. The app can be run using e.g.     
~~~sh
portageapp_jali \
    --dim=3 \
    --nsourcecells=12 \
    --ntargetcells=8 \
    --field='x+2*y-3*z' \
    --remap_order=2
~~~

More options can by explored by executing 
~~~sh
portageapp_jali --help

Usage: portageapp --dim=2|3 --nsourcecells=N --ntargetcells=M --conformal=y|n 
--entity_kind=cell|node --field="your_math_expression" --remap_order=1|2 
--limiter=barth_jespersen --mesh_min=0. --mesh_max=1. 
--output_meshes=y|n --results_file=filename --convergence_study=NREF 
~~~

## Portageapp_multimat_jali

Similarly to `portageapp_jali` for scalar fields, `app/portageapp/portageapp_multimat_jali` 
demonstrates Portage's ability to remap multi-material fields. In order to preserve sharp
interfaces between materials, Portage relies on [Tangram](https://github.com/laristra/tangram) 
to perform an interface reconstruction
on a source mesh. Reconstructed material polygons and further intersected with cells on a target 
mesh. Therefore, an input file specifying volume fractions (and optionally material centroids 
for MOF) is required in addition to specification of the meshes and multiple material fields. 
The file can be generated e.g. by a `vfgen` application in [Tangram](https://github.com/laristra/tangram).
~~~sh
portageapp_multimat_jali \
    --dim=2 \
    --source_file=reg10x10.exo \
    --target_file=polymesh_10x10.exo \
    --material_fields=5*x-y,2*x+4+3*y,-2*x-y \
    --material_file=reg10x10.bvf 
~~~ 
