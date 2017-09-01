# Simple Mesh Example    {#example}

portage provides very crude mesh and state manager frameworks aptly
called `Simple_Mesh` and `Simple_State`.  The goal of these frameworks
is to show how one can wrap their favorite mesh and state manager for
use with portage - _they should not be used in any production sense._

----

# The Mesh and State

## Portage::Simple_Mesh

This mesh framework is a non-distributed (i.e. has no ghost
information), 3D, regular Cartesian mesh framework.  A Simple_Mesh is
constructed by specifying the extents of the box and the number of
cells in each direction.  The constructor then builds connectivity
information between cell, node, and face indices.  The ordering of
things like the nodes making up a cell, or the faces making up a cell
are specified consistently, but the choice of ordering does not
matter.  There are a few helper functions like
Portage::Simple_Mesh::cell_get_nodes() that will retrieve the indices
of some connected mesh entities given another mesh entity's index.

## Portage::Simple_State

The state manager for Portage::Simple_Mesh is essentially a collection
of field data specified by some name (e.g. "density") and location on
the Portage::Simple_Mesh where they live (e.g. Portage::CELL).  The
constructor for a Portage::Simple_State takes a pointer to a
Portage::Simple_Mesh so that it has access to things like the number
of nodes in the mesh.  Data are added to and retrieved from the state
manager via `add` and `get` methods.  Iterators are provided for the
map between `(name, location)` pairs and the data.

# Wrappers

As mentioned on the [Concepts](@ref concepts) page, the search,
intersect, and interpolate actions in portage all operate on
mesh/particle and state _wrappers_.  The reason for this is that these
three steps may need to ask for some information that may not be
readily available within the original mesh/particle and state
frameworks, but could be constructed within a wrapper.

We provide a helper class, Portage::AuxMeshTopology, that assists in
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

## Portage::AuxMeshTopology

This helper class will build the additional mesh entities and
connectivities from a basic mesh framework wrapper.  The basic mesh
framework wrapper must support cells, nodes, and faces as well as
connectivity and number queries about these entities.

This class is not a complete class design.  In particular, it is
designed to be used within
the
[Curiously Recurring Template Pattern](https://en.m.mwikipedia.org/wiki/Curiously_recurring_template_pattern)(CRTP)
design pattern to achieve static polymorphism.  Under the CRTP in this
case, the basic mesh framework wrapper looks something like

~~~{.cc}
class Basic_Mesh_Wrapper : public AuxMeshTopology<Basic_Mesh_Wrapper> {...};
~~~

In this way, the `Basic_Mesh_Wrapper` can use its own methods, or
defer to AuxMeshTopology to perform more advanced queries.

In addition to the advanced mesh entities (sides, wedges, and
corners), AuxMeshTopology also resolves some advanced connectivity
information.  An example is
Portage::AuxMeshTopology::node_get_cell_adj_nodes(), which, given a node index
in the mesh, returns a vector of all the nodes that are attached to
all cells attached to the given node.  AuxMeshTopology additionally
creates iterators over all the various types of mesh entities,
provides geometric information (e.g. volumes, centroids, etc.) of the
mesh entities, and a method to determine if an entitiy is on the
domain boundary (needed for limiting in higher-order remaps).

_If one does **not** use AuxMeshTopology to help extend their mesh
wrapper's functionality, one should ensure that their mesh wrapper at
least has the same public functions as AuxMeshTopology._ For this
reason, it is advised to utilize AuxMeshTopology where possible, but
to defer to one's own mesh framework wrapper when more efficient
methods are available.

## Portage::Simple_Mesh_Wrapper

This wraps the Portage::Simple_Mesh framework.  Simple_Mesh, as its
name suggests, is quite simple and does not know about advanced mesh
entities nor connectivities.  It lets AuxMeshTopology do the vast
majority of the heavy lifting by automatically creating the advanced
mesh features.

Where possible, Simple_Mesh_Wrapper provides quick and efficient
answers to queries that AuxMeshTopology would otherwise solve in a
general sense.  Two trivial examples are:

1. Portage::Simple_Mesh_Wrapper::cell_get_type(), which determines the
   Portage::Entity_type (e.g. PARALLEL_OWNED, PARALLEL_GHOST, etc.).
   We know Portage::Simple_Mesh does not know anything about ghost
   information, so we simple return
   Portage::Entity_type::PARALLEL_OWNED.
2. Portage::Simple_Mesh_Wrapper::cell_get_element_type(), which
   determines the geometric shape of a given cell from one of the
   Portage::Element_type's.  Simple mesh is only a 3d, structured
   Cartesian mesh, so we always return Portage::Element_type::HEX.

There are a few other examples within Portage::Simple_Mesh_Wrapper
where the AuxMeshTopology methods are overwritten to take advantage of
information that is cached within the Portage::Simple_Mesh.  This is a
prime example of how the CRTP and AuxMeshTopology are intended to be
used.  In fact, even our wrapper (wonton::flecsi_mesh_t) to the
sophisticated Burton
mesh [specialization](https://github.com/laristra/flecsi-sp) of
the [FleCSI](https://github.com/laristra/flecsi) mesh framework uses
AuxMeshTopology to answer some queries.

## Portage::Simple_State_Wrapper

There is no equivalent of AuxMeshTopology for state wrappers.  This is
simply because the requirements of a state manager are much less
intense than those of the mesh framework.  In paticular, the state
wrappers only need to know how to add data, get data, and query things
like the size of the data and where the data lives (Portage::CELL or
Portage::NODE).  Portage::Simple_State_Wrapper exposes this
functionality and provides some error checking for missing or
duplicate data in the underlying Portage::Simple_State object.

# Applications and Tests

Nearly all of the generic unit tests have been designed to work with
the Portage::Simple_Mesh_Wrapper and Portage::Simple_State_Wrapper
data structures when a full up mesh is required.  These should serve
as examples of the types of things one may want to do with a mesh and
state framework.

Better examples of how to use these wrappers to actually do a remap of
field data are in the application tests.  In particular, the
`app/simple_mesh_app/simple_mesh_app.cc` program shows how to wrap
mesh and state objects, adds some field data to the source state, and
utilize Portage::Driver with various search, intersect, and
interpolate algorithms to perform the remap.  The Portage::Driver is
templated on mesh wrapper type and state wrapper type, and can be used
with other frameworks.  It need not be used at all, but is a nice
starting point for writing one's own remap application.
