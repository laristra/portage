# Simple Mesh Example    {#simple_mesh}

portage provides very crude mesh and state manager frameworks, aptly called
`Simple_Mesh` and `Simple_State`.  The goal of these frameworks is to show
how you can wrap your favorite mesh and state manager for use with portage.

## High-level portage Workflow

portage provides library functions for remapping data from one mesh (the
_source_ mesh) to a a different mesh (the _target mesh_).  In order to do this,
portage needs to _search_ through the source mesh for possible intersections
with each mesh entity in the target mesh, actually calculate the amount of
_intersection_ with the candidates, and then _interpolate_ the weighted results
onto the target mesh.

## Why Wrappers?

In order to perform the _search_, _intersect_, and _interpolate_ steps of a
remap, portage needs to ask the mesh and statemanager for information about
data, connectivity/neighbors, etc.  In order to do this in a general way,
we have adopted wrappers around mesh frameworks that should implement the
queries required by portage.  It is then up to the user of the mesh and state
manager of choice to make sure those queries are satisfied, in whatever
(read: optimal) choice they choose.

## `AuxMeshTopology`

To facilitate wrapping, we provide a class called Portage::AuxMeshTopology,
which implements
the
[Curiously Recurring Template Pattern](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to
implement static polymorphism.  In short, this means that
AuxMeshTopology is templated on a mesh wrapper class, which should be
able to answer basic questions about cells and nodes (and faces and
edges, depending on the dimension).  AuxMeshTopology then can build
more advanced mesh connectivity information, such as sides, wedges,
and corners (see the *Simple Mesh* section below).

A concrete example is the Portage::AuxMeshTopology::num_entities()
function.  This function takes a Portage::Entity_kind (
e.g. Portage::CELL or Portage::NODE) and a Portage::Entity_type
(e.g. Portage::PARALLEL_OWNED or Portage::PARALLEL_GHOST) and returns
the appropriate number of entities of that type and kind.  Under the
hood, it calls the base mesh wrapper class' functions
(e.g. Portage::Simple_Mesh::num_owned_faces()) where applicable, or
can call its own functions for information about more advanced
entities (e.g. Portage::AuxMeshTopology::num_ghost_sides()).

**NOTE:** When wrapping a new mesh framework, one does not _need_ to use
AuxMeshTopology, especially if your mesh framework already supports
the advanced mesh entities; it is simply a convenience tool.

## Methods Wrapped

This is a list of methods needed by a mesh wrapper during a typical
flow through portage.  This is the list of things needed with and in addition to
using Portage::AuxMeshTopology class.  The examples here are concrete examples
using the Portage::Simple_Mesh class as a reference.

<dl>
	<dt>General mesh info</dt>
		<dd>Portage::Simple_Mesh::space_dimension()</dd>
		<dd>Portage::Simple_Mesh::cell_get_type()</dd>
		<dd>Portage::Simple_Mesh::cell_get_element_type()</dd>
		<dd>Portage::Simple_Mesh::get_global_id()</dd>
	<dt>`num_X_Y` for various entities</dt>
		<dd>Portage::Simple_Mesh::num_owned_cells()</dd>
		<dd>Portage::Simple_Mesh::num_owned_faces()</dd>
		<dd>Portage::Simple_Mesh::num_owned_nodes()</dd>
		<dd>Portage::Simple_Mesh::num_ghost_cells()</dd>
		<dd>Portage::Simple_Mesh::num_ghost_faces()</dd>
		<dd>Portage::Simple_Mesh::num_ghost_nodes()</dd>
	<dt>Connectivity information</dt>
		<dd>Portage::Simple_Mesh::node_get_cell_adj_nodes()</dd>
		<dd>Portage::Simple_Mesh::cell_get_node_adj_cells()</dd>
		<dd>Portage::Simple_Mesh::cell_get_faces_and_dirs()</dd>
		<dd>Portage::Simple_Mesh::cell_get_nodes()</dd>
		<dd>Portage::Simple_Mesh::face_get_nodes()</dd>
	<dt>Spatial information</dt>
		<dd>Portage::Simple_Mesh::node_get_coordinates()</dd>
</dl>

## Simple Mesh

## Simple State
