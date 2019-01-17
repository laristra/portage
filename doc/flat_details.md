# Flat Mesh and State Details {#flat_details}

Portage can perform distributed remap even when the source and target partitioning is mismatched i.e. the source cells overlapping the cells of a target partition are on a different partition. It does this by communicating the bounding box of a target partition to every other partition (All-to-All). If a source partition on, say rank 0, overlaps a target partition bounding box from a different rank, say rank 1, then the source partition from rank 0 is sent to rank 1 and "concatenated" to the already existing source partition on that rank. Clearly, there are inefficiencies in this process and a more selective filtering of which cells to send to a different rank would be warranted. However, for now, we will not focus on future improvements but focus instead on the nitty-gritty's of the algorithm as it is implemented. The code is somewhat involved and often takes several days of staring at it to divine what's going on. 

To understand what is being done in the code, we will do a thought experiment with a remap between two 1D meshes distributed over two ranks as shown in the figure below (*Portage does not yet have a 1D capability*). The global mesh domain goes from 0.0 to 1.0 and it is assumed that all the nodes are equally spaced. The meshes on the two ranks are shown with the correct spatial overlap along the X-axis. The solid lines are owned cells and the dashed lines are ghost cells. These meshes are shown below with the following meaning for labels: <code>n<sub>i</sub></code> denotes the local number of a node on a rank, <code>N<sub>i</sub></code>  denotes the **global** number of a node in the mesh, <code>c<sub>i</sub></code> denotes the local number of a cell on a rank, <code>C<sub>i</sub></code> denotes the **global** number of a cell in the mesh.

<br/><br/>

<p align="center">
  <img src="parallel-mesh.png" width="500" title="Parallel Meshes">
</p>
<br/><br/>
Then on rank 1, the list of global cell IDs and list of global node IDs for the source will look like this:<br/><br/>

Array Index	| 0	| 1	| 2
--- | --- | --- | ---
<b>CellGlobalIDs</b>	| <span style="color:red">C<sub>2</sub> | C<sub>3</sub> | C<sub>1</sub></span>

	

	


Array Index
	0	1	2	3
NodeGlobalIDs	

	

	

	


The pink boxes indicate ghost cells which always come after owned cells of a partition; similarly, ghost nodes.


Looking further, the cellToNodeCounts array may look like this
CellNodeCounts	2	2	2


and the cellToNodeList array may look like this, right after the moveField call
CellToNodeList	

	

	

	


	

	


The CellToNodeList is referring to local node indices



On rank 0, for the source we have
Array Index	0	1	2
CellGlobalIDs	

	

	


Array Index
	0	1	2	3
NodeGlobalIDs	

	

	

	


The pink boxes indicate ghost cells which always come after owned cells of a partition; similarly, ghost nodes.


Looking further, the cellToNodeCounts array may look like this
CellNodeCounts	2	2	2


and the cellToNodeList array may look like this, right after the moveField call
CellToNodeList	

	

	

	

	

	




The first step is to create a flat_mesh_wrapper and flat_state_wrapper object from the regular mesh_wrapper and state wrapper classes. The flat_mesh_wrapper and flat_state_wrapper are special mesh/state wrappers templated on the regular mesh/state wrappers that facilitate the redistribution of data across processor (see flat_mesh_wrapper.h and flat_state_wrapper.h or flat_state_mm_wrapper.h from the Wonton repository). The flat_mesh_wrapper in particular copies mesh data from an arbitrary mesh wrapper into its own flat data structures that can be sent over MPI and updated easily.

We should check for bounding box overlap (see below) BEFORE we create flat_mesh_wrapper to avoid unnecessary copy if there is no overlap and nothing needs to be sent over. The bounding box check needs nothing special from the flat mesh wrappers


Then the code creates an MPI_distributor (See mpi_bounding_boxes.h) that does the job of sending and receiving data across and reconciling the data. The finish_init method of the flat_mesh_wrapper is then called to finalize the process of reconciling the mesh from different source ranks
MPI_Bounding_Boxes::distribute(...)


MPI_Bounding_Boxes should be given a more general name


For each target partition, we compute its bounding box (that includes only the owned target cells) and communicate it to all the other processors. These bounding boxes are then compared with the bounding box of each source partition (including owned and ghost cells).

If the source bounding box on a rank overlaps the target bounding box from another rank, a flag is set to indicate that the source partition from this rank has to be sent over to the target rank. Thus each rank has a set of sendFlags and recvFlags indicating which ranks it is sending data to and which ranks it is receiving data from.


 Looking at the example, since the target bounding box from rank 1 [0.33,1.0] overlaps the source bounding box on rank 0 [0.0,0.75], the source partition from rank 0 must be sent to rank 1. Similarly since  target bounding box from rank 0 [0.0,0.33] overlaps the source bounding box from rank 1 [0.25,1.0], the source partition from rank 1 has to be sent to rank 0. This example illustrates that our bounding box check scoops up too many partitions. In reality, the target partitions on both ranks are fully contained within the bounding boxes of their own source partitions.



We really should use something more fine grained than bounding boxes to check for overlap.


We really should use a finer-grained filter to send only relevant cells from a source partition to a target partition instead of sending the entire source partition over.


The actual sending and receiving of mesh and field information is accomplished by the functions setInfo and moveField. For each type of data that needs to be communicated (e.g. global IDs of cells, or IDs of nodes of cells), setInfo communicates the counts of owned and (owned+ghost) entities that will be sent and received. It helps to populate a data structure of type comm_info_t which carries these counts. Then moveField actually sends the data around (actually moveField calls moveData twice, once for owned entities and once for ghost entities). The data received from each processor is concatenated to the existing data so that the first set is always data that is local to this partition




Referring back to our example and focusing our attention on rank 1, our source cell global ID and source node global ID lists look like
Array Index	0	1	2	3	4	5
CellGlobalIDs	


	


	


	


	


	



Array Index
	0	1	2	3	4	5	6	7
NodeGlobalIDs	

	

	

	

	

	

	

	


The pink boxes indicate ghost cells which always come after owned cells of a partition; similarly, ghost nodes.


Looking further, the cellToNodeCounts array may look like this
CellNodeCounts	2	2	2	2	2	2


and the cellToNodeList array may look like this, right after the moveField call
CellToNodeList	

	

	

	

	

	

	

	

	

	

	

	



Clearly, the local node indices are wrong in the rank 1 part of the cellToNodeList of the concatenated list as they refer to indices on rank 1. So, at this stage, the routine fixListIndices is called to find the offsets for cells and nodes for each rank in the concatenated list.



Referring back to our example, we can see that the offset for where the imported nodes from rank 0 is stored on rank 1 is 0 and the offset for the "native" nodes is 4. Therefore all the local node indices in the rank 1 part of CellToNodeList are bumped up by 4 giving us
CellToNodeList	

	

	

	

	

	

	

	

	

	

	

	


The code now calls the finish_init method of the flat_mesh_wrapper which does some de-duplication of topological references so that the flat mesh wrapper can answer mesh questions meaningfully.
Flat_Mesh_Wrapper::finish_init()

The first thing that finish_init does is to make a unique global to local map. To do this it actually uses the map data structure from standard C++ so that when we can easily check if a global ID is within the map before trying to add it. At the same time, it also creates a helper array for de-duplication, called cellUniqueRep and nodeUniqueRep. These point from local IDs of duplicate cells/nodes to the first occurrence of that cell/node.



In our example, this results in a globalCellMap data structure on rank 1 that looks like this
globalCellMap	

	

	

	


	

0
	

1
	

2
	

4

and a globalNodeMap that looks like this
nodalCellMap	

	

	

	

	


	0	1	2	

3
	

6


The cellUniqRep data array looks like this
Array Index	0	1	2	3	4	5
cellUniqRep	0	1	2	2	4	1

This array tells us that local cell 3 is the duplicate of 2 and cell 5 is the duplicate of 1. We can refer back to our list of global IDs to see that this is true.


Similarly, the nodeUniqRep array looks like this
Array Index	0	1	2	3	4	5	6	7
nodeUniqRep	0	1	2	3	2	3	6	1




Using the nodeUniqRep array, one can now simplify the CellToNodeList array to refer to the first occurrence of the nodes. This will be useful later in building the correct nodeToCellList arrays. Also, the node counts of duplicate cells are set to 0 so that they are not double counted during the actual remap operations.



So now our CellToNodeList array looks like this
CellToNodeList	

	

	

	

	

	

	

	

	

	

	

	

Now if we try to build a NodeToCellList from the CellToNodeList we will get the right upward adjacencies (these are used in the AuxMeshTopology class to get the cells surrounding a cell which is useful for gradient calculations).

At the same time the CellNodeCounts list is updated to be
CellNodeCounts	2	2	2	2	0	0

so that if we scoop up all 6 source cells in a search for intersection candidates for a target cell, we will not account for the last two which are duplicates


Now the flat_mesh_wrapper will answer the basic mesh wrapper questions correctly for the concatenated meshes. Finally, AuxMeshTopology::build_aux_entities() method is called to build up additional adjacencies.

finish_init with all of its baggage is also called when we call initialize for the flat_mesh_wrapper and it ends up doing a lot of redundant work related to deduplication just so that it can build a few connectivities - not sure which ones. It should be broken up into parts that deduplicate and parts that build these connectivities.



When we do a remap, the "run" routine in the driver code takes a distributed flag as an argument. If this flag is true the remap is done using remap_distributed instead of the plain remap. The only place flat mesh and flat state are used is in remap_distributed. In short, a serial remap will never use flat mesh or state. Flat mesh and state have a single purpose and that is to enable distributed remap.


The life cycle of the flat mesh is as follows. The flat mesh is initialized from the underlying mesh wrapper on node. At this point the two meshes are equivalent conceptually. Bounding boxes are computed for all partitions on their respective nodes. The bounding box information is distributed from every node to every node and flags are determined for whether every source and target partition intersect. These flags are reused for the remainder of the data transfer and if a source and target pair intersect, then data is transferred from source to target.  This step of transferring data both about bounding boxes and then ultimately the data is done by the mpi_bounding_boxes.h routine named "distribute". After the distribute step, all the source data required by a target partition is located on node. The issue is that there may be duplicates because of ghost cells. A source cell on a target node may be owned by one source partition, a ghost on any number of other partitions, or only exist as a ghost cell. Cells can be identified by a unique global id so it is possible to tell when a cell is duplicated on the target partition.


The flat mesh is used in two ways that are actually quite distinct. Prior to distribution, each source partition needs to know its own bounding box in order to determine overlaps with target partitions. That's it, nothing more is required. After distribution, then a full mesh API is required in order to do the remap.


Flat meshup to a point works quite well. If source and target partitions overlap then a flag is set and all data from the source is sent to the target. There is a long discussion that can be had as to whether this is a good thing to do. The existing method is very simple but will send too much data. If source and target only intersect in a single cell, why send the entire partition? However, we should defer this discussion to a later date. The current goal is to get MM distributed remap working.


Another discussion we should defer is about ghosts. When computing what data should be sent, we use the bounding box of owned cells, that means ghosts are excluded. The ghosts are used to compute gradients and nothing more. Flat mesh and distribute devote a fair amount of code to handling ghosts. Essentially all data transfer is doubled, once for owned cells and once for ghosts. There is separate bookkeeping for owned versus ghost cells and the counts are kept separately. On the target, in the flat mesh the data is striped [owned partition1, ghost partition1, owned partition 2, ghost partition2, etc]. Once cells are moved to the target partition, their "ghostness" becomes irrelevant. The interpolation does not care whether the cell is/was a ghost. Again, any source cell may appear as owned, a ghost or both. The discussion to be had is whether we could simplify the distribution by keeping the bounding box computation over owned cells but sending a single list of cells both owned and ghost, losing the distinction between them and still being able to do the necessary remap steps. I (DWS) believe this is possible since the flat mesh only needs to support the "remap API" and not a full mesh topological API which would distinguish between owned and ghost cells. Flat mesh and state have a single purpose which is to enable distributed remap and the interface could be economized to support this.


Finally, we get to the primary issue with the flat mesh at the moment and why it needs to fixed. Flat mesh is redundant. It has multiple copies of identical cells which may be owned or ghost in any combination. Therefore a flat mesh is too big in that it has redundant cells. The problem comes to when we try to do interface reconstruction. Tangram requires the length of data (volume fractions and centroids) to be the same as the number of cells in the mesh. One option is to redundantly include the data for duplicated cells, but interface reconstruction is expensive and why do it multiple times when it will never be used? (RVG: Does the rest of the code work if you do this?


Distributing the data across partitions has a post processing step which finds unique cells based on global cell id. Any cell in the flat mesh is referenced to the first occurrence that is equivalent. With Portage's material dominant data representation each material has a list of cells. Tangram requires a cell dominant data representation. Each cells needs a list of materials in it. The problem from flat state is that with multiple equivalent cells being added to flat state we end up with repeated cells in the cell lists. When inverted to the cell centric for Tangram there are now cells (the non primary duplicates) that have no material data... and Tangram breaks.


The proposed solution is to simply disambiguate the cells in flat mesh. After the distribution, which again works well, we simply need to merge unique cells based on global id, so each cell appears only once. The mesh and state will be correctly constructed with no missing cell material or duplicated material cell ids. Removing the duplication is straightforward. Finish_init already constructs a map based on global id that determines the correct cell to reference for any global id. The idea is when constructing private data members in flat_mesh:


private:

  std::vector<T>    nodeCoords_;
  std::vector<int>  cellToNodeList_;
  std::vector<int>  cellNodeCounts_;
  std::vector<int>  cellNodeOffsets_;
  std::vector<int>  cellToFaceList_;
  std::vector<bool> cellToFaceDirs_;  // unused in 2D (identical 
                                      // to cellNodeCounts_)
  std::vector<int>  cellFaceCounts_;  // unused in 2D (identical 
                                      // to cellNodeOffsets_)
  std::vector<int>  cellFaceOffsets_;
  std::vector<int>  faceToNodeList_;
  std::vector<int>  faceNodeCounts_;  // unused in 2D (always 2)
  std::vector<int>  faceNodeOffsets_; // unused in 2D (can be computed)
  std::vector<int>  nodeToCellList_;
  std::vector<int>  nodeCellCounts_;
  std::vector<int>  nodeCellOffsets_;
  std::vector<int>  cellGlobalIds_;
  std::vector<int>  nodeGlobalIds_;


We don't loop over all cells in the post distribute flat mesh, we only loop over unique global ids (I would suggest in sort order to facilitate debugging). Using this strategy everything can be consistently book kept.


This strategy does require a data copy. The various private data members cannot be kept intact. However, with redundancy in the cells, extra work was already being done. The only vectors that need to be modified in place are a subset of those listed. The others are constructed in finish_init and AuxMeshTopology.h.


After sleeping on this, I have a few more observations. Flat mesh is really a generic mesh, like simple mesh in that it supplies the functionality required for remap while keeping a local copy of the data. The internal representation of the mesh being composed entirely of standard vectors is inherently suited to distribution. The distribution is handled by the mpi_bounding_boxes code, not the mesh itself. There is nothing wrong with this, it is just the way we implement. However, once on the target partition, no effort is made to turn the mesh data into a topologically correct mesh.


A different organization of the code might be (analogous to the flat state code) a mesh which inherits from a generic mesh but provides serialization and deserialization capability. That is really what we are doing although the routines are not named as such. Also it is worth noting that only a minimal amount of information is distributed to construct a topologically correct mesh, namely:

    node coordinates
    cell node counts (2D only)
    cell node lists (2D only)
    cell face counts (3D only)
    cell face lists (3D only, data is mangled to include binary direction as well)
    cell face node counts (3D only)
    cell face node lists (3D only)
    cell global ids
    node global ids

The private data members which are offsets are computed as a postprocessing step. The inverse maps of nodeToCellList_, nodeCellCounts_, nodeCellOffsets_ are also computed. In 2D, cellToFaceList_ and cellToFaceDirs_ are computed as well.


A note of curiosity is that internally the cell to face index data transfer is hijacked to transfer both the face id and direction into a single number. The original face id is bit shifted left and the boolean direction is inserted as the first bit.


As my thinking progresses on this, my greatest concern is how to handle faces in 3D. In 2D faces are constructed on the fly. In 3D they are given by the underlying meshes. At present, faces do not have a global id (in the sense that we don't get it, the AuxMeshTopology spec requires it to be available) and therefore cannot be merged across partitions without doing some computationally expensive or hacked step (e.g. set equality). Face global ids are also not distributed. If the underlying mesh API can give a global id for faces as well, then we can distribute and merge just like the other topological entities. This should be the case, but since I have not seen a face global id, I  cannot be sure.



