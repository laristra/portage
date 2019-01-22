# High Level Concepts of Data Distribution During Parallel Remap {#distributed_concepts}

When running on multiple processors, Portage uses a very simple data
distribution scheme which will be described in this document. 

At the beginning of the distributed remap step in Portage, which for the purposes
of this document we will consider in the default remap case to be the call to 
```
remap_distributed(...)
```
within the basic driver call of
```
driver.run(distributed)
```
, each processor node has a partition of both the source and target meshes. The 
spatial extent of the source and target meshes
on a given partition need not have any special relationship. Portage uses 
polygonal intersections of the overlap between source and target cells (or dual
cells in the case of nodal remap) to compute interpolated field values. From this
point forward we will use the word "cell" to refer to either a cell in cell-based
remap or a dual cell in node-based remap.

In order for each target cell to compute interpolated field values it is necessary
to know the geometries of all source cells that intersect it. This presents a 
problem when as mentioned, the source and target meshes on a given processor 
node are independent. The essence of Portage's data distribution is to bring all
source cells that could intersect the target cells on a particular partition onto
the target partition. Portage's paradigm of "search, intersect, interpolate" 
requires that any source cells that could potentially intersect a target cell
be available for computation on that target cell's partition. In 
```
remap_distributed(...)
```
data distribution is the first nontrivial thing done so that the "search" step
has all the data necessary to work correctly. After data distribution, each
processor node operates exactly as if it were a serial problem. This is the 
source of Portage's embarrasingly parallel characterization.


Portage uses a simple algorithm for data distribution which is described 
presently. A bounding box is computed for the source and target meshes on each
partition. The bounding boxes are computed by finding the minimum and maximum
coordinate in each direction for each mesh on each partition. An MPI_Bcast is 
used to distribute 
the bounding boxes from all target meshes to each source partition. Using the 
bounding boxes, it is determined which source partitions intersect 
each target partition. Bounding box intersections are used to create a vector of booleans, 
different on each partition, that determine to which target partitions each source 
partition should send its data. The same communication topology is used for all
subsequent data communication. Clearly this sends too much data as only a few
source cells may intersect target cells, but the algorithm is easy to easy to
compute.

Many pieces of data are sent from souce partition to target partition. For
each of these pieces of information, two MPI calls need to be made. The first
call establishes the counts of the data that will be sent. The
second MPI call does the actual data distribution. 

It is worthwhile to note that most data is sent making a distinction between 
ghost cells and owned cells. Both ghost and owned data are sent in the same
request, but the counts are kept separately and the sent data
is ordered by owned data first followed by ghost data.
 
The data distributed for remap includes global ids for all entities, node 
coordinates, adjacency information within the mesh and field values, which is
potentially multi-material. There are
many subtleties in distributing the data. The four most important of which are
removing duplicates, converting local references, handling vector data 
 and multi-material data. 

The first subtlety is that because of the way that distribution is handled, 
namely whole patches are sent, the same cell can be received from multiple ranks.
A cell can be only be owned by one rank, but may be a ghost on multiple other
ranks. A cell may also only exist as a ghost on any rank. The global id is part
of the data sent from each rank, so duplicated entities can be removed based on
only keeping the data for each global identifier once. 

The second subtlety is that each reference in the adjacency data is represented
inherently by local references (local ids). The same global id for each mesh
entity (cells, faces and nodes) can be used to update the references to the
new indexing scheme on the partition.


The third subtlety is that MPI data transfer only works with sending vectors of 
uniform type. So any data
transfer needs to conform to this description. In other words, if we need to 
transfer data that is inherently itself vector data, such as a cell centroid, 
it must first be "flattened" to a vector of doubles where the individual 
components of the cell's data are listed in order. Such data must be serialized
to the form (e.g. in 2D)
```
(x1, y1, x2, y2, ..., xn, yn).
```
Upon arrival on the target mesh, the data must then be deserialized to the
original vector form, still with duplicates removed.

The final subtlety is that multimaterial fields greatly complicate the data 
distribution. Multimaterial data in Portage is stored in material dominant 
form, meaning a vector of cells and cell data is kept for each material. Material
dominant data looks like:
```
material 1: cell_id_1, cell_id_3, ...
material 2: cell_id_1, cell_id_4, ...
...
```
And the data
```
material 1: x11, x13, ...
material 2: x21, x24, ...
...
```
Where in the various ```xij's```, the ```x``` data corresponds to material ```i```
in cell ```j```. Multimaterial data greatly complicates the bookkeeping of doing
data distribution because each material acts like it's own mini-mesh. All
multimaterial data is assumed to exist on the same set of cells, so the 
multimaterial cell data is shared between fields. This implies 
that if a material exists in a cell, all multi-material field values must be 
known for that 
material in that cell. Distributing multimaterial field data requires sending 
for each partition the number of materials, the material ids, the number of cells
having each material, the cell ids, then finally the data itself. Just as with 
all other data, since MPI can handle only vectors of uniform type, each of these
"ragged right" data structures must be serialized to vectors of common type on
the source partition before sending and deserialized from vectors to the ragged
right form required by Portage. Finally, multimaterial data is inherently more
complicated because of its shape but we still need to remove duplicated data
and update local ids to the new partition.

