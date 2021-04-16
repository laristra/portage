/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MPI_BOUNDING_BOXES_H_
#define MPI_BOUNDING_BOXES_H_


#ifdef WONTON_ENABLE_MPI

#include <cassert>
#include <algorithm>
#include <numeric>
#include <memory>
#include <unordered_map>
#include <map>
#include <vector>
#include <set>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/state/state_vector_uni.h"
#include "portage/support/portage.h"
#include "mpi.h"

/*!
  @file mpi_bounding_boxes.h
  @brief Distributes source data using MPI based on bounding boxes
 */

namespace Portage {

using Wonton::Point;
using Wonton::GID_t;
using Wonton::to_MPI_Datatype;

/*!
  @class MPI_Bounding_Boxes
  @brief Distributes source data using MPI based on bounding boxes

         Currently assumes coordinates and all fields are doubles.
*/
class MPI_Bounding_Boxes {
 public:

  /*!
    @brief Constructor of MPI_Bounding_Boxes
   */
  MPI_Bounding_Boxes(Wonton::MPIExecutor_type const *mpiexecutor) {
    assert(mpiexecutor);
    comm_ = mpiexecutor->mpicomm;
  }


  /*!
    @brief Destructor of MPI_Bounding_Boxes
   */
  ~MPI_Bounding_Boxes() = default;


  /*!
    @brief Helper structure containg comms info for a given entity type
   */
  struct comm_info_t {
    //!< Number of total/owned entities in source field
    int sourceNum = 0, sourceNumOwned = 0;
    //!< Number of total/owned entities in new field
    int newNum = 0, newNumOwned = 0;
    //! Array of total/owned send sizes from me to all PEs
    std::vector<int> sendCounts {}, sendOwnedCounts {};
    //! Array of total/owned recv sizes to me from all PEs
    std::vector<int> recvCounts {}, recvOwnedCounts {};
  };


  /*!
    @brief Compute whether this partition (Bob) needs data from other partitions 
          (hungry) or whether all the data is already on the partition
    @param[in] source_mesh  Input mesh 
    @param[in] target_mesh  Target mesh

   */
  template <class Source_Mesh, class Target_Mesh>
  bool is_bob_hungry(const Source_Mesh &source_mesh,
                  const Target_Mesh &target_mesh)
  {
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(comm_, &commSize);
    MPI_Comm_rank(comm_, &commRank);

    dim_ = source_mesh.space_dimension();
    assert(dim_ == static_cast<int>(target_mesh.space_dimension()));

    // sendFlags, which partitions to send data
    // this is computed via intersection of whole partition bounding boxes
    std::vector<bool> sendFlags(commSize);
    compute_sendflags(source_mesh, target_mesh, sendFlags);
    
    std::vector<int> sendFlagsInt(commSize);
    for (int i=0; i<commSize; ++i) sendFlagsInt[i]=sendFlags[i]?1:0;
    
    // Each rank will tell each other rank if it is going to send it
    std::vector<int> recvFlags(commSize);
    MPI_Alltoall(&(sendFlagsInt[0]), 1, MPI_INT,
                 &(recvFlags[0]), 1, MPI_INT, comm_);
                 
    // loop over the partitions
    for (int i=0; i<commSize; ++i) 
      // if any other other partition sends me data, we need to redistribute
      if (i!=commRank && recvFlags[i])
        return true;  
  
    // the fall through case means no redistribution
    return false;
  }


  /*!
    @brief Compute whether any partition needs data from other partitions
    @param[in] source_mesh  Input mesh 
    @param[in] target_mesh  Target mesh

   */
  template <class Source_Mesh, class Target_Mesh>
  bool is_redistribution_needed(const Source_Mesh &source_mesh,
                  const Target_Mesh &target_mesh)
  {
    // does this partition need data from an other partition
    bool r = is_bob_hungry(source_mesh, target_mesh);
    
    // convert to an int 
    int ir = r ? 1 : 0;
    
    // allocate the result
    int result;
    
    // do the MPI_ALLReduce to aggregate the results
    MPI_Allreduce(&ir, &result, 1, MPI_INT, MPI_LOR, comm_);
                    
    return result;
  }



  /*!
    @brief Compute bounding boxes for all partitions, and send source mesh and state
           information to all target partitions with an overlapping bounding box using MPI
    @param[in] source_mesh_flat  Input mesh (must be flat representation)
    @param[in] source_state_flat Input state (must be flat representation)
    @param[in] target_mesh       Target mesh
    @param[in] target_state      Target state (not actually used for now)

   */
  template <class Source_Mesh, class Source_State, class Target_Mesh, class Target_State>
  void distribute(Source_Mesh &source_mesh_flat, Source_State &source_state_flat,
                  Target_Mesh &target_mesh, Target_State &target_state)
  {
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(comm_, &commSize);
    MPI_Comm_rank(comm_, &commRank);

    int dim = dim_ = source_mesh_flat.space_dimension();
    assert(dim == static_cast<int>(target_mesh.space_dimension()));

    // sendFlags, which partitions to send data
    // this is computed via intersection of whole partition bounding boxes
    std::vector<bool> sendFlags(commSize);
    compute_sendflags(source_mesh_flat, target_mesh, sendFlags);

    // set counts for cells
    comm_info_t cellInfo;
    int sourceNumOwnedCells = source_mesh_flat.num_owned_cells();
    int sourceNumCells = sourceNumOwnedCells + source_mesh_flat.num_ghost_cells();
    setSendRecvCounts(&cellInfo, commSize, sendFlags,sourceNumCells, sourceNumOwnedCells);

    // set counts for nodes
    comm_info_t nodeInfo;
    int sourceNumOwnedNodes = source_mesh_flat.num_owned_nodes();
    int sourceNumNodes = sourceNumOwnedNodes + source_mesh_flat.num_ghost_nodes();
    setSendRecvCounts(&nodeInfo, commSize, sendFlags,sourceNumNodes, sourceNumOwnedNodes);

    ///////////////////////////////////////////////////////
    // always distributed
    ///////////////////////////////////////////////////////

    // SEND GLOBAL CELL IDS

    std::vector<GID_t>& sourceCellGlobalIds = source_mesh_flat.get_global_cell_ids();
    std::vector<GID_t> distributedCellGlobalIds(cellInfo.newNum);
    sendField(cellInfo, commRank, commSize, to_MPI_Datatype<GID_t>(), 1,
              sourceCellGlobalIds, &distributedCellGlobalIds);

    // SEND GLOBAL NODE IDS
    std::vector<GID_t>& sourceNodeGlobalIds = source_mesh_flat.get_global_node_ids();
    std::vector<GID_t> distributedNodeGlobalIds(nodeInfo.newNum);
    sendField(nodeInfo, commRank, commSize, to_MPI_Datatype<GID_t>(), 1,
              sourceNodeGlobalIds, &distributedNodeGlobalIds);

    // Using the post distribution global id's, compress the data so that each
    // global id appears only once. The trick here is to get the ghosts correct.
    // In the flat mesh, after distribution, an entity that was owned by any
    // partition is considered owned in the flat mesh. Likewise, any entity that
    // appears only as a ghost will be a ghost in the flat mesh
    compress_with_ghosts(distributedCellGlobalIds, cellInfo.newNumOwned,
      distributedCellIds_, flatCellGlobalIds_, flatCellNumOwned_);
    compress_with_ghosts(distributedNodeGlobalIds, nodeInfo.newNumOwned,
      distributedNodeIds_, flatNodeGlobalIds_, flatNodeNumOwned_);

    // create the map from cell global id to flat cell index
    create_gid_to_flat_map(flatCellGlobalIds_, gidToFlatCellId_);
    create_gid_to_flat_map(flatNodeGlobalIds_, gidToFlatNodeId_);

    // SEND NODE COORDINATES
    std::vector<double>& sourceCoords = source_mesh_flat.get_coords();
    std::vector<double> distributedCoords(dim*nodeInfo.newNum);
    sendField(nodeInfo, commRank, commSize, MPI_DOUBLE, dim,
              sourceCoords, &distributedCoords);

    // merge and set coordinates in the flat mesh
    merge_duplicate_data(distributedCoords, distributedNodeIds_, sourceCoords, dim_);



    ///////////////////////////////////////////////////////
    // 2D distributed
    ///////////////////////////////////////////////////////

    if (dim == 2)
    {

      // send cell node counts
      std::vector<int>& sourceCellNodeCounts = source_mesh_flat.get_cell_node_counts();
      std::vector<int> distributedCellNodeCounts(cellInfo.newNum);
      sendField(cellInfo, commRank, commSize, MPI_INT, 1,
                sourceCellNodeCounts, &distributedCellNodeCounts);

      // merge and set cell node counts
      merge_duplicate_data(distributedCellNodeCounts, distributedCellIds_, sourceCellNodeCounts);

      // mesh data references
      std::vector<int>& sourceCellNodeOffsets = source_mesh_flat.get_cell_node_offsets();
      std::vector<int>& sourceCellToNodeList = source_mesh_flat.get_cell_to_node_list();

      int sizeCellToNodeList = sourceCellToNodeList.size();
      int sizeOwnedCellToNodeList = (
          sourceNumCells == sourceNumOwnedCells ? sizeCellToNodeList :
          sourceCellNodeOffsets[sourceNumOwnedCells]);

      comm_info_t cellToNodeInfo;
      setSendRecvCounts(&cellToNodeInfo, commSize, sendFlags,
              sizeCellToNodeList, sizeOwnedCellToNodeList);

      // send cell to node lists
      std::vector<GID_t> distributedCellToNodeList(cellToNodeInfo.newNum);
      sendField(cellToNodeInfo, commRank, commSize, to_MPI_Datatype<GID_t>(), 1,
                to_gid(sourceCellToNodeList, sourceNodeGlobalIds), &distributedCellToNodeList);


      // merge and map cell node lists
      merge_duplicate_lists(distributedCellToNodeList, distributedCellNodeCounts,
        distributedCellIds_, gidToFlatNodeId_, sourceCellToNodeList);

    }


    ///////////////////////////////////////////////////////
    // 3D distributed
    ///////////////////////////////////////////////////////

    if (dim == 3)
    {

      int sourceNumOwnedFaces = source_mesh_flat.num_owned_faces();
      int sourceNumFaces = sourceNumOwnedFaces + source_mesh_flat.num_ghost_faces();

      comm_info_t faceInfo;
      setSendRecvCounts(&faceInfo, commSize, sendFlags,sourceNumFaces, sourceNumOwnedFaces);

      // SEND GLOBAL FACE IDS
      std::vector<GID_t>& sourceFaceGlobalIds = source_mesh_flat.get_global_face_ids();
      std::vector<GID_t> distributedFaceGlobalIds(faceInfo.newNum);
      sendField(faceInfo, commRank, commSize, MPI_LONG_LONG, 1,
                sourceFaceGlobalIds, &distributedFaceGlobalIds);

      // Create map from distributed gid's to distributed index and flat indices
      compress_with_ghosts(distributedFaceGlobalIds, faceInfo.newNumOwned,
        distributedFaceIds_, flatFaceGlobalIds_, flatFaceNumOwned_);

      // create the map from face global id to flat cell index
      create_gid_to_flat_map(flatFaceGlobalIds_, gidToFlatFaceId_);

      // mesh data references
      std::vector<int>& sourceCellFaceOffsets = source_mesh_flat.get_cell_face_offsets();
      std::vector<int>& sourceCellToFaceList = source_mesh_flat.get_cell_to_face_list();

      int sizeCellToFaceList = sourceCellToFaceList.size();
      int sizeOwnedCellToFaceList = (
          sourceNumCells == sourceNumOwnedCells ? sizeCellToFaceList :
          sourceCellFaceOffsets[sourceNumOwnedCells]);

      comm_info_t cellToFaceInfo;
      setSendRecvCounts(&cellToFaceInfo, commSize, sendFlags,
              sizeCellToFaceList, sizeOwnedCellToFaceList);

      // SEND NUMBER OF FACES FOR EACH CELL
      std::vector<int>& sourceCellFaceCounts = source_mesh_flat.get_cell_face_counts();
      std::vector<int> distributedCellFaceCounts(cellInfo.newNum);
      sendField(cellInfo, commRank, commSize, MPI_INT, 1,
                sourceCellFaceCounts, &distributedCellFaceCounts);

      // merge and set cell face counts
      merge_duplicate_data( distributedCellFaceCounts, distributedCellIds_, sourceCellFaceCounts);

      // SEND CELL-TO-FACE MAP
      // map the cell face list vector to gid's
      std::vector<GID_t> sourceCellToFaceList_ = to_gid(sourceCellToFaceList, sourceFaceGlobalIds);

      // For this array only, pack up face IDs + dirs and send together
      std::vector<bool>& sourceCellToFaceDirs = source_mesh_flat.get_cell_to_face_dirs();
      int const sourceCellToFaceListSize = sourceCellToFaceList.size();

      for (int j = 0; j < sourceCellToFaceListSize; ++j) {
        int f = sourceCellToFaceList_[j];
        int dir = static_cast<int>(sourceCellToFaceDirs[j]);
        sourceCellToFaceList_[j] = (f << 1) | dir;
      }

      std::vector<GID_t> distributedCellToFaceList(cellToFaceInfo.newNum);
      sendField(cellToFaceInfo, commRank, commSize, MPI_LONG_LONG, 1,
                sourceCellToFaceList_, &distributedCellToFaceList);

      // Unpack face IDs and dirs
      std::vector<bool> distributedCellToFaceDirs(cellToFaceInfo.newNum);
      int const distributedCellToFaceListSize = distributedCellToFaceList.size();

      for (int j = 0; j < distributedCellToFaceListSize; ++j) {
        int fd = distributedCellToFaceList[j];
        distributedCellToFaceList[j] = fd >> 1;
        distributedCellToFaceDirs[j] = fd & 1;
      }


      // merge and map cell face lists
      merge_duplicate_lists(distributedCellToFaceList, distributedCellFaceCounts,
        distributedCellIds_, gidToFlatFaceId_, sourceCellToFaceList);

      // merge cell face directions
      merge_duplicate_lists(distributedCellToFaceDirs, distributedCellFaceCounts,
        distributedCellIds_, sourceCellToFaceDirs);

      // mesh data references
      std::vector<int>& sourceFaceNodeOffsets = source_mesh_flat.get_face_node_offsets();
      std::vector<int>& sourceFaceToNodeList = source_mesh_flat.get_face_to_node_list();

      int sizeFaceToNodeList = sourceFaceToNodeList.size();
      int sizeOwnedFaceToNodeList = (
          sourceNumFaces == sourceNumOwnedFaces ? sizeFaceToNodeList :
          sourceFaceNodeOffsets[sourceNumOwnedFaces]);

      comm_info_t faceToNodeInfo;
      setSendRecvCounts(&faceToNodeInfo, commSize, sendFlags,
              sizeFaceToNodeList, sizeOwnedFaceToNodeList);

      // SEND NUMBER OF NODES FOR EACH FACE
      std::vector<int>& sourceFaceNodeCounts = source_mesh_flat.get_face_node_counts();
      std::vector<int> distributedFaceNodeCounts(faceInfo.newNum);
      sendField(faceInfo, commRank, commSize, MPI_INT, 1,
                sourceFaceNodeCounts, &distributedFaceNodeCounts);

      // SEND FACE-TO-NODE MAP
      std::vector<GID_t> distributedFaceToNodeList(faceToNodeInfo.newNum);
      sendField(faceToNodeInfo, commRank, commSize, MPI_LONG_LONG, 1,
                to_gid(sourceFaceToNodeList, sourceNodeGlobalIds), &distributedFaceToNodeList);

      // merge and set face node counts
      merge_duplicate_data( distributedFaceNodeCounts, distributedFaceIds_, sourceFaceNodeCounts);

      // merge and map face node lists
      merge_duplicate_lists(distributedFaceToNodeList, distributedFaceNodeCounts,
        distributedFaceIds_, gidToFlatNodeId_, sourceFaceToNodeList);

      // merge face global ids
      merge_duplicate_data(distributedFaceGlobalIds, distributedFaceIds_,sourceFaceGlobalIds);

      // set counts for faces in the flat mesh
      source_mesh_flat.set_num_owned_faces(flatFaceNumOwned_);

    }

    // SEND FIELD VALUES

    // multimaterial state info
    int nmats = source_state_flat.num_materials();
    comm_info_t num_mat_cells_info {};

    // Is the a multimaterial problem? If so we need to pass the cell indices
    // in addition to the field values
    if (nmats>0){
#if defined(PORTAGE_DEBUG)
      std::cout << "in distribute, this a multimaterial problem with " << nmats << " materials\n";
#endif

      /////////////////////////////////////////////////////////
      // get the material ids across all nodes
      /////////////////////////////////////////////////////////

      // set the info for the number of materials on each node
      comm_info_t num_mats_info;
      setSendRecvCounts(&num_mats_info, commSize, sendFlags, nmats, nmats);

      // get the sorted material ids on this node
      std::vector<int> material_ids=source_state_flat.get_material_ids();

      // get all material ids across
      distributedMaterialIds_.resize(num_mats_info.newNum);

      // send all materials to all nodes, num_mats_info.recvCounts is the shape
      sendData(commRank, commSize, MPI_INT, 1, 0, num_mats_info.sourceNum, 0,
        num_mats_info.sendCounts, num_mats_info.recvCounts,
        material_ids, &distributedMaterialIds_
      );

      /////////////////////////////////////////////////////////
      // get the material cell shapes across all nodes
      /////////////////////////////////////////////////////////

      // get the sorted material shapes on this node
      std::vector<int> material_shapes=source_state_flat.get_material_shapes();

      // resize the post distribute material id vector
      distributedMaterialShapes_.resize(num_mats_info.newNum);

      // send all material shapes to all nodes
      sendData(commRank, commSize, MPI_INT, 1, 0, num_mats_info.sourceNum, 0,
        num_mats_info.sendCounts, num_mats_info.recvCounts,
        material_shapes, &distributedMaterialShapes_
      );

      /////////////////////////////////////////////////////////
      // get the lists of material cell ids across all nodes
      /////////////////////////////////////////////////////////

      // get the total number of material cell id's on this node
      int nmatcells = source_state_flat.num_material_cells();

      // set the info for the number of materials on each node
      setSendRecvCounts(&num_mat_cells_info, commSize, sendFlags, nmatcells, nmatcells);

      // get the sorted material ids on this node
      std::vector<int> material_cells=source_state_flat.get_material_cells();

      // resize the post distribute material id vector
      distributedMaterialCells_.resize(num_mat_cells_info.newNum);

      // send material cells to all nodes, but first translate to gid
      sendData(commRank, commSize, to_MPI_Datatype<GID_t>(), 1, 0,
        num_mat_cells_info.sourceNum, 0, num_mat_cells_info.sendCounts,
        num_mat_cells_info.recvCounts, to_gid(material_cells, sourceCellGlobalIds),
        &distributedMaterialCells_
      );

      /////////////////////////////////////////////////////////
      // We need to turn the flattened material cells into a correctly shaped
      // ragged right structure for use as the material cells in the flat
      // state wrapper. We aren't removing duplicates yet, just concatenating by material
      /////////////////////////////////////////////////////////

      /*!
        @todo It might be possible to make the next two code blocks more
        performant by using some builtin such as insert, std::copy, or transform.
        The first block is likely, although I am not sure of the improvement since it
        is a reshape and would require as many rounds of resizing as partitions.
        The current way uses the intrinsic vector management and is probably more
        efficient. The second block needs to convert from gid to flat local id,
        so that is less promising... Maybe a transform. Don't know if it will be
        any faster though.
      */

      // allocate the material indices
      std::unordered_map<int,std::vector<GID_t>> material_indices;

      // reset the running counter
      int running_counter=0;

      // loop over material ids on different nodes
      int const distributedMaterialIdsSize = distributedMaterialIds_.size();

      for (int i = 0; i < distributedMaterialIdsSize; ++i) {

        // get the current working material
        int mat_id = distributedMaterialIds_[i];

        // get the current number of material cells for this material
        int nmat_cells = distributedMaterialShapes_[i];

        // get or create a reference to the correct material cell vector
        std::vector<GID_t>& these_material_cells = material_indices[mat_id];

        // loop over the correct number of material cells
        for (int j=0; j<nmat_cells; ++j){
            these_material_cells.push_back(distributedMaterialCells_[running_counter++]);
        }
      }

      // cell material indices are added not replaced by vector so we need to
      // start clean
      source_state_flat.clear_material_cells();

      // merge the material cells and convert to local id (gid sort order)
      for (auto &kv: material_indices){

        // get the material id
        int m = kv.first;

        // compress the material ids so each gid appears only once per material
        compress(kv.second, distributedMaterialCellIds_[m]);

        // allocate the flat material cell ids for this material
        std::vector<int> flatMaterialCellIds;
        flatMaterialCellIds.reserve(distributedMaterialCellIds_[m].size());

        // loop of material cell indices, converting to gid, then flat cell
        for (auto id: distributedMaterialCellIds_[m])
          flatMaterialCellIds.push_back(gidToFlatCellId_[kv.second[id]]);

        // add the material cells to the state manager
        source_state_flat.mat_add_cells(kv.first, flatMaterialCellIds);

      }
    }

    // Send and receive each field to be remapped
    for (std::string field_name : source_state_flat.names())
    {

      // this is a packed version of the field with copied field values and is
      // not a pointer to the original field
      std::vector<double> sourceField = source_state_flat.pack(field_name);

      // get the field stride
      int sourceFieldStride = source_state_flat.get_field_stride(field_name);

      comm_info_t info;
      if (source_state_flat.get_entity(field_name) == Entity_kind::NODE){
          // node mesh field
          info = nodeInfo;
      } else if (source_state_flat.field_type(Entity_kind::CELL, field_name) == Wonton::Field_type::MESH_FIELD){
          // mesh cell field
          info = cellInfo;
      } else {
         // multi material field
         info = num_mat_cells_info;
      }

      // allocate storage for the new distributed data, note that this data
      // has raw doubles and will need to be merged and type converted

      std::vector<double> distributedField(sourceFieldStride*info.newNum);

      // send the field
      sendField(info, commRank, commSize, MPI_DOUBLE, sourceFieldStride,
        sourceField, &distributedField);

      std::vector<double> tempDistributedField;

      if (source_state_flat.get_entity(field_name) == Entity_kind::NODE){

        // node mesh field
        // merge duplicates, but data still is raw doubles
        merge_duplicate_data(distributedField, distributedNodeIds_, tempDistributedField, sourceFieldStride);

        // unpack the field, has the correct data type
        source_state_flat.unpack(field_name, tempDistributedField);

      } else if (source_state_flat.field_type(Entity_kind::CELL, field_name) == Wonton::Field_type::MESH_FIELD){

        // cell mesh field
        // merge duplicates, but data still is raw doubles
        merge_duplicate_data(distributedField, distributedCellIds_, tempDistributedField, sourceFieldStride);

        // unpack the field, has the correct data types
        source_state_flat.unpack(field_name, tempDistributedField);

      } else {

        // multi material field
        // as opposed to the preceeding two cases, the merging and type conversion
        // are both done in the unpack routine because there is more to do
        // getting the shapes correct.
        source_state_flat.unpack(field_name, distributedField,
          distributedMaterialIds_, distributedMaterialShapes_,
          distributedMaterialCellIds_);
      }
    }

    // need to do this at the end of the mesh stuff, because converting to
    // gid uses the global id's and we don't want to modify them before we are
    // done converting the old relationships
    // merge global ids and set in the flat mesh
    merge_duplicate_data(distributedCellGlobalIds, distributedCellIds_, sourceCellGlobalIds);
    merge_duplicate_data(distributedNodeGlobalIds, distributedNodeIds_, sourceNodeGlobalIds);

    // set counts for cells and nodes in the flat mesh
    source_mesh_flat.set_num_owned_cells(flatCellNumOwned_);
    source_mesh_flat.set_num_owned_nodes(flatNodeNumOwned_);

    // Finish initialization using redistributed data
    source_mesh_flat.finish_init();


  } // distribute

  private:

  // The communicator we are using
  MPI_Comm comm_ = MPI_COMM_NULL;

  int dim_ = 1;

  // the number of nodes "owned" by the flat mesh. "Owned" is in quotes because
  // a node may be "owned" by multiple partitions in the flat mesh. A node is
  // owned by the flat mesh if it was owned by any partition.
  int flatNodeNumOwned_ = 0;

  // the global id's of the kept nodes in the flat mesh and their indices in the
  // distributed node global id's
  std::vector<GID_t> flatNodeGlobalIds_ {};
  std::vector<int> distributedNodeIds_ {};

  // maps from gid to distributed node index and flat node index
  std::map<GID_t,int> gidToFlatNodeId_ {};

  // the number of faces "owned" by the flat mesh. "Owned" is in quotes because
  // a face may be "owned" by multiple partitions in the flat mesh. A face is
  // owned by the flat mesh if it was owned by any partition.
  int flatFaceNumOwned_ = 0;

  // the global id's of the kept faces in the flat mesh and their indices in the
  // distributed face global id's
  std::vector<GID_t> flatFaceGlobalIds_ {};
  std::vector<int> distributedFaceIds_ {};

  // maps from gid to distributed face index and flat face index
  std::map<GID_t,int> gidToFlatFaceId_ {};

  // the number of cells "owned" by the flat mesh. "Owned" is in quotes because
  // a cell may be "owned" by multiple partitions in the flat mesh. A cell is
  // owned by the flat mesh if it was owned by any partition.
  int flatCellNumOwned_ = 0;

  // the global id's of the kept cells in the flat mesh and their indices in the
  // distributed cell global id's
  std::vector<GID_t> flatCellGlobalIds_ {};
  std::vector<int> distributedCellIds_ {};

  // maps from gid to distributed cell index and flat cell index
  std::map<GID_t,int> gidToFlatCellId_ {};

  // vectors for distributed multimaterial data
  std::vector<int> distributedMaterialIds_ {};
  std::vector<int> distributedMaterialShapes_ {};
  std::vector<GID_t> distributedMaterialCells_ {};

  // map for the distributed material cell indices
  // for each material there is a vector of unique distributed indices
  std::map<int, std::vector<int>> distributedMaterialCellIds_ {};

  /*!
    @brief Compute fields needed to do comms for a given entity type
    @param[in] info              Info data structure to be filled
    @param[in] commSize          Total number of MPI ranks
    @param[in] sendFlags         Array of flags:  do I send to PE n?
    @param[in] sourceNum         Number of entities (total) on this rank
    @param[in] sourceNumOwned    Number of owned entities on this rank
   */
  void setSendRecvCounts(comm_info_t* info,
               const int commSize,
               const std::vector<bool>& sendFlags,
               const int sourceNum,
               const int sourceNumOwned)
  {
    // Set my counts of all entities and owned entities
    info->sourceNum = sourceNum;
    info->sourceNumOwned = sourceNumOwned;

    // Each rank will tell each other rank how many indexes it is going to send it
    info->sendCounts.resize(commSize);
    info->recvCounts.resize(commSize);
    for (int i=0; i<commSize; i++)
      info->sendCounts[i] = sendFlags[i] ? info->sourceNum : 0;
    MPI_Alltoall(&(info->sendCounts[0]), 1, MPI_INT,
                 &(info->recvCounts[0]), 1, MPI_INT, comm_);

    // Each rank will tell each other rank how many owned indexes it is going to send it
    info->sendOwnedCounts.resize(commSize);
    info->recvOwnedCounts.resize(commSize);
    for (int i=0; i<commSize; i++)
      info->sendOwnedCounts[i] = sendFlags[i] ? info->sourceNumOwned : 0;
    MPI_Alltoall(&(info->sendOwnedCounts[0]), 1, MPI_INT,
                 &(info->recvOwnedCounts[0]), 1, MPI_INT, comm_);

    // Compute the total number of indexes this rank will receive from all ranks
    for (int i=0; i<commSize; i++)
      info->newNum += info->recvCounts[i];
    for (int i=0; i<commSize; i++)
      info->newNumOwned += info->recvOwnedCounts[i];
  } // setSendRecvCounts


  /*!
    @brief Send values for a single data field to all ranks as needed
    @tparam[in] T                C++ type of data to be sent
    @param[in] info              Info struct for entity type of field
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be sent
    @param[in] stride            Stride of data field
    @param[in] sourceData        Array of (old) source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void sendField(const comm_info_t& info,
                 int commRank, int commSize,
                 MPI_Datatype mpiType, int stride,
                 const std::vector<T>& sourceData,
                 std::vector<T>* newData)
  {
    // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
    std::vector<int> recvGhostCounts(commSize);
    std::vector<int> sendGhostCounts(commSize);

    for (int i=0; i<commSize; i++)
    {
      recvGhostCounts[i] = info.recvCounts[i] - info.recvOwnedCounts[i];
      sendGhostCounts[i] = info.sendCounts[i] - info.sendOwnedCounts[i];
    }

    sendData(commRank, commSize, mpiType, stride,
             0, info.sourceNumOwned,
             0,
             info.sendOwnedCounts, info.recvOwnedCounts,
             sourceData, newData);
    sendData(commRank, commSize, mpiType, stride,
             info.sourceNumOwned, info.sourceNum,
             info.newNumOwned,
             sendGhostCounts, recvGhostCounts,
             sourceData, newData);

#ifdef DEBUG_MPI
    std::cout << "Number of values on rank " << commRank << ": " << (*newData).size() << std::endl;
    for (int f=0; f<(*newData).size(); f++)
      std::cout << (*newData)[f] << " ";
    std::cout << std::endl << std::endl;
#endif
  } // sendField


  /*!
    @brief Send values for a single range of data to all ranks as needed
    @tparam[in] T                C++ type of data to be sent
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be sent
    @param[in] stride            Stride of data field
    @param[in] sourceStart       Start location in source (send) data
    @param[in] sourceEnd         End location in source (send) data
    @param[in] newStart          Start location in new (recv) data
    @param[in] curSendCounts     Array of send sizes from me to all PEs
    @param[in] curRecvCounts     Array of recv sizes to me from all PEs
    @param[in] sourceData        Array of (old) source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void sendData(int commRank, int commSize,
                MPI_Datatype mpiType, int stride,
                int sourceStart, int sourceEnd,
                int newStart,
                const std::vector<int>& curSendCounts,
                const std::vector<int>& curRecvCounts,
                const std::vector<T>& sourceData,
                std::vector<T>* newData)
  {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = newStart;
    int myOffset = 0;
    std::vector<MPI_Request> requests;
    for (int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curRecvCounts[i] > 0))
      {
        MPI_Request request;
        MPI_Irecv((void *)&((*newData)[stride*writeOffset]),
                  stride*curRecvCounts[i], mpiType, i,
                  MPI_ANY_TAG, comm_, &request);
        requests.push_back(request);
      }
      else if (i == commRank)
      {
        myOffset = writeOffset;
      }
      writeOffset += curRecvCounts[i];
    }

    // Copy data values that will stay on this rank into the
    // proper place in the new vector
    if (curRecvCounts[commRank] > 0)
      std::copy(sourceData.begin() + stride*sourceStart,
                sourceData.begin() + stride*sourceEnd,
                newData->begin() + stride*myOffset);

    // Each rank will send its data values to appropriate ranks
    for (int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curSendCounts[i] > 0))
      {
        MPI_Send((void *)&(sourceData[stride*sourceStart]),
                 stride*curSendCounts[i], mpiType, i, 0, comm_);
      }
    }

    // Wait for all receives to complete
    if (!requests.empty())
    {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }
  } // sendData


  template <class Source_Mesh, class Target_Mesh>
  void compute_sendflags(Source_Mesh & source_mesh, Target_Mesh &target_mesh,
              std::vector<bool> &sendFlags){

    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(comm_, &commSize);
    MPI_Comm_rank(comm_, &commRank);

    int dim = source_mesh.space_dimension();

    // Compute the bounding box for the target mesh on this rank, and put it in an
    // array that will later store the target bounding box for each rank
    int targetNumOwnedCells = target_mesh.num_owned_cells();
    std::vector<double> targetBoundingBoxes(2*dim*commSize);
    for (int i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }
    for (int c=0; c<targetNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      target_mesh.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        // ugly hack, since dim is not known at compile time
        if (dim == 3)
        {
          Point<3> nodeCoord;
          target_mesh.node_get_coordinates(n, &nodeCoord);
          for (int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = nodeCoord[k];
            if (nodeCoord[k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = nodeCoord[k];
          }
        }
        else if (dim == 2)
        {
          Point<2> nodeCoord;
          target_mesh.node_get_coordinates(n, &nodeCoord);
          for (int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = nodeCoord[k];
            if (nodeCoord[k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = nodeCoord[k];
          }
        }
      } // for j
    } // for c

    // Get the source mesh properties and cordinates
    int sourceNumOwnedCells = source_mesh.num_owned_cells();

    // Compute the bounding box for the source mesh on this rank
    std::vector<double> sourceBoundingBox(2*dim);
    for (int i=0; i<2*dim; i+=2)
    {
      sourceBoundingBox[i+0] = std::numeric_limits<double>::max();
      sourceBoundingBox[i+1] = -std::numeric_limits<double>::max();
    }
    for (int c=0; c<sourceNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      source_mesh.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        if (dim ==3 )
        {
          Point<3> nodeCoord;
          source_mesh.node_get_coordinates(n, &nodeCoord);
          for (int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < sourceBoundingBox[2*k])
              sourceBoundingBox[2*k] = nodeCoord[k];
            if (nodeCoord[k] > sourceBoundingBox[2*k+1])
              sourceBoundingBox[2*k+1] = nodeCoord[k];
          }
        }
        else if (dim ==2)
        {
          Point<2> nodeCoord;
          source_mesh.node_get_coordinates(n, &nodeCoord);
          for (int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < sourceBoundingBox[2*k])
              sourceBoundingBox[2*k] = nodeCoord[k];
            if (nodeCoord[k] > sourceBoundingBox[2*k+1])
              sourceBoundingBox[2*k+1] = nodeCoord[k];
          }
        }
      } // for j
    } // for c

    // Broadcast the target bounding boxes so that each rank knows the bounding boxes for all ranks
    for (int i=0; i<commSize; i++)
      MPI_Bcast(&(targetBoundingBoxes[0])+i*2*dim, 2*dim, MPI_DOUBLE, i, comm_);


    // Offset the source bounding boxes by a fudge factor so that we don't send source data to a rank
    // with a target bounding box that is incident to but not overlapping the source bounding box
    const double boxOffset = 2.0*std::numeric_limits<double>::epsilon();
    double min2[dim], max2[dim];
    for (int k=0; k<dim; ++k)
    {
      min2[k] = sourceBoundingBox[2*k]+boxOffset;
      max2[k] = sourceBoundingBox[2*k+1]-boxOffset;
    }

    // For each target rank with a bounding box that overlaps the bounding box for this rank's partition
    // of the source mesh, we will send it all our source cells; otherwise, we will send it nothing
    for (int i=0; i<commSize; ++i)
    {
      double min1[dim], max1[dim];
      bool sendThis = true;
      for (int k=0; k<dim; ++k)
      {
        min1[k] = targetBoundingBoxes[2*dim*i+2*k]+boxOffset;
        max1[k] = targetBoundingBoxes[2*dim*i+2*k+1]-boxOffset;
        sendThis = sendThis &&
            ((min1[k] <= min2[k] && min2[k] <= max1[k]) ||
             (min2[k] <= min1[k] && min1[k] <= max2[k]));
      }

      sendFlags[i] = sendThis;
    }
  }


  /*!
    @brief Convert a vector of integer references to their global id's.

    @param[in] in  Integer references to convert to gid
    @param[in] gids  The vector of gid's for each entity pre distribution
    @return The new vector of gid's after mapping the references

    This is always used prior to distribution. The idea is that in a topological
    map such as sourceCellToNodeList, the second type of entity, node in this
    case, needs to get converted to gid. We always convert local ids to gids
    before distributing.
   */
  std::vector<GID_t> to_gid(std::vector<int> const& in, std::vector<GID_t>const& gids) const {
    std::vector<GID_t> result;
    result.reserve(in.size());
    for (auto x:in) result.push_back(gids[x]);
    return result;
  }


  /*!
    @brief Compress distributed global id's into the vector of distributed
    indices(first occurrence for each global id).
    The global id's are in ascending order.

    @param[in] distributedGlobalIds  The vector of gid's for each entity post distribution
    @param[out] distributedIds  The vector of indices of first occurrence in
      distributedGlobalIds that produce the global id's in the flat mesh

    This is called after distribution of gids for each set of material cell gids.
    The vector of material cell gids (with potentially duplicate gid's) is used
    to create a vector of he first position of occurrence in post distribution
    material cell gids vector.

    The following example will hopefully illustrate what this function is doing.
    Consider the following post distribution set of material cell gids:
    `3 | 2 7 3 | 9 10 7 3`
    where a single vertical bar separates partitions. The original distribution
    returns all entities in rank order. So the example above consists of three
    partitions. The first partition has entity `3`, the second partition has
    entities `2, 7, 3`, and the third partition has entities `9,10, 7, 3`. After
    merging, in the flat mesh, the  entities will `2, 3, 7, 9, 10`, but we don't
    need to return this. The result of the function is
    `distributedIds=[1, 0, 2, 4, 5]` . The  returned argument, `distributedIds`
    is the indices of the unique global ids in the distributed
    material cell data. The vector `distributedIds` is used to merge all
    subsequent distributed material celldata.
  */
  void compress(std::vector<GID_t> const& distributedGlobalIds,
                std::vector<int>& distributedIds) const {

    // a map for keeping track of what we have already seen
    std::map<GID_t,int> uniqueGid;

    int const num_distributed_global_ids = distributedGlobalIds.size();

    // loop over owned cells in the distributed global id's
    for (int i = 0 ; i < num_distributed_global_ids; ++i){

      // get the current gid
      GID_t gid = distributedGlobalIds[i];

      // is this gid new
      if (uniqueGid.find(gid) == uniqueGid.end()){

        // make the gid seen
        uniqueGid[gid]=i;

      }

    }

    // clear and reserve the result
    distributedIds.clear();
    distributedIds.reserve(uniqueGid.size());

    // push the  cells in gid order
    for (auto const& kv : uniqueGid){

       // push to the distributed cells id
      distributedIds.push_back(kv.second);

    }
  }
  /*!
    @brief Compress distributed global id's into the vector of distributed
    indices(first occurrence for each global id) and the vector of global id's.
    The global id's are in ascending order of "owned" global id, followed by "ghost"
    global id. In the flat mesh, an entitiy is considered "owned" if it is
    owned by any partition, and a "ghost" if it is not owned by any partition.

    @param[in] distributedGlobalIds  The vector of gid's for each entity post distribution
    @param[in] distributedNumOwned  The number of entities that are owned after distribution
    @param[out] distributedIds  The vector of indices of first occurrence in
      distributedGlobalIds that produce the global id's in the flat mesh
    @param[out] flatGlobalIds  The global ids in the flat mesh. They are in order
      of ascending gid for owned entities followed by ascending gid for ghosts
    @param[out] flatNumOwned  The number of owned entities in the flat mesh

    This is called after distribution of gids for each entity kind. The vector of
    gids (with potentially duplicate gid's) is used to create two vectors. The first
    vector is the first position of occurrence in post distribution gids vector.
    The second vector is the gid in the flat mesh after compression.

    We handle"ghosts" in the post distributed data. An entity is considered
    "owned" by the flat mesh if it was owned on any partition and a "ghost" if it
    was a ghost on all partitions. This implies that an entity will be owned by as many
    partitions as the original source partition was sent to. So while we preserve
    the concept that a ghost is a ghost in the flat mesh only if it was originally
    only ghosts, we unfortunately have that an entity in the flat mesh will typically
    be owned by many different paritions.

    The following example will hopefully illustrate what this function is doing.
    Consider the following post distribution set of gids:
    `3 | 2 7 | 9 10 || 4* 2* | 5* 4* | 4*`
    where a single vertical bar separates partitions, the double vertical bar
    separates the owned entities from ghosts and a number followed by an asterisk
    means it was a ghost on it's original partition. The original distribution
    returns all owned entities in rank order followed by all ghost entities in
    rank order. So the example above consists of three partitions. The first
    partition has owned entity `3` and ghost entities `4,2`. The second partition
    has owned entities `2,7` and ghost entities `5,4`. The third partition has
    owned entities `9,10` and ghost entity `4`. After merging, in the flat mesh,
    the owned entities will `2, 3, 7, 9, 10` and the ghost entities will be
    `4, 5`. The interesting entity is `2` which is both owned and a ghost on
    different partitions and is therefore considered owned. The result of the
    function is `distributedIds=[1, 0, 2, 3, 4, 5,7]` and
    `flatGlobalIds=[2, 3, 7, 9, 10, 4, 5]`. The first returned argument,
    `distributedIds` is the indices of the unique global ids in the distributed
    data. The second returned argument `flatGlobalIds` is the gids themselves.
    In the code, `distributedIds` is used much more frequently than `flatGlobalIds`.
    The vector `distributedIds` is used to merge all subsequent distributed
    data and topological reference vectors. The number of owned entities would
    be `flatNumOwned=5`. Notice that in `flatGlobalIds` the vector consists
    of two ascending sequences `2, 3, 7, 9, 10` followed by `4, 5`. The first
    sequence is owned entities. The second is ghosts.
  */
  void compress_with_ghosts(std::vector<GID_t> const& distributedGlobalIds,
    int const distributedNumOwned, std::vector<int>& distributedIds,
    std::vector<GID_t>& flatGlobalIds, int &flatNumOwned) const {

    // a map for keeping track of what we have already seen
    std::map<GID_t,int> uniqueGid, uniqueGhostGid;

    // loop over owned entitites in the distributed global id's
    for (int i=0 ; i<distributedNumOwned ; ++i){

      // get the current gid
      GID_t gid = distributedGlobalIds[i];

      // is this gid new
      if (uniqueGid.find(gid)==uniqueGid.end()){

        // make the gid seen
        uniqueGid[gid]=i;

      }

    }

    // We have processed owned entitites in the distributed mesh, so everything
    // we have collected to this point is considered owned
    flatNumOwned = uniqueGid.size();

    // push the owned entitites first and in gid order
    for (auto const& kv : uniqueGid){

      // push to the flat entitites gid
      flatGlobalIds.push_back(kv.first);

      // push to the distributed entitites id
      distributedIds.push_back(kv.second);

    }

    int const num_distributed_global_ids = distributedGlobalIds.size();

    // loop over ghost entitites in the distributed global id's
    for (int i = distributedNumOwned ; i < num_distributed_global_ids; ++i){

      // get the current gid
      GID_t gid = distributedGlobalIds[i];

      // is this gid new
      if (uniqueGid.find(gid)==uniqueGid.end() &&
        uniqueGhostGid.find(gid)==uniqueGhostGid.end()){

        // make the gid seen
        uniqueGhostGid[gid]=i;

      }

    }

    // push the ghost entities in gid order
    for (auto const& kv : uniqueGhostGid){

      // push to the flat entities gid
      flatGlobalIds.push_back(kv.first);

      // push to the distributed entities id
      distributedIds.push_back(kv.second);

    }

  }


  /*!
    @brief Create a map from gid to flat mesh index

    @param[in] flatGlobalIds  The vector of gid's in the flat mesh
    @param[out] gidToFlat  The map from global id to flat cell index

    This function creates a trivial map from global id to flat cell index.
    The index of the global id is the value of the map. We use this later when
    converting global id's into the new flat mesh local index.

    In the above example, create_gid_flat_map would take the vector
    `flatGlobalIds=[2, 3, 7, 9, 10, 4, 5]` and return the map
    `gidToFlat={(2:0),(3:1),(7:2),(9:3),(10:4),(4:5),(5:6)}`.
  */
  void create_gid_to_flat_map(std::vector<GID_t> const& flatGlobalIds,
                              std::map<GID_t,int>& gidToFlat) const {

    int const num_flat_global_ids = flatGlobalIds.size();
    for (int i = 0; i < num_flat_global_ids; ++i)
      gidToFlat[flatGlobalIds[i]]=i;

  }


  /*!
    @brief Merge post distribution data so that each datum appears
           only once for each unique gid.

    @param[in] in   The post distribution data to be merged
    @param[in] distributedIds  The vector of distributed indices to keep
    @param[out] result The new vector of flattened and merged data
    @param[in] stride=1  The number of values associated with the each element of data.
                       Scalar data has stride 1, centroid data has stride D.

    This is called after distribution. It takes post distribution data and
    selects a single occurrence of each gid defined by distributedIds
   */
  template<class T>
  void merge_duplicate_data(std::vector<T>const& in, std::vector<int>const& distributedIds,
    std::vector<T>& result, int stride = 1){

    // clear result and reserve to stride * the number kept
    result.clear();
    result.reserve(distributedIds.size()*stride);

    // since the vector is the correct size
    for (auto id: distributedIds)
      for (int d=0; d<stride; ++d)
        result.push_back(in[stride*id + d]);
  }

  /*!
    @brief Merge post distribution lists of topological references.

    @param[in] in   The post distribution references to be merged
    @param[in] counts  The vector of the number of references associated with
                       this entity, e.g. number of faces for this cell
    @param[in] distributedIds  The map from gid to first occurrence pre distribution of
                         the "from" part of the mapping, e.g. cell in cellToFace
    @param[in] gidToFlatId  The map from gid to position in the post distribution
                         vector of the "to" part of the mapping, e.g. face in
                         cellToFace
    @param[out] result The new vector of merged references

    This is called after distribution. It takes post distribution topological
    reference data and corrects for two things. The first thing is that the input has
    duplicate offset lists that need to be merged. The second thing is that the
    topological references need to get converted from gid to their new flat index id.
  */
  void merge_duplicate_lists(std::vector<GID_t>const& in, std::vector<int> const& counts,
    std::vector<int>const& distributedIds, std::map<GID_t,int>const& gidToFlatId,
    std::vector<int>& result){

    // allocate offsets
    std::vector<int> offsets(counts.size());

    // compute offsets (note the first element is zero and correct)
    std::partial_sum(counts.begin(), counts.end()-1, offsets.begin()+1);

    // make sure the result is clear and approximately sized
    result.clear();
    result.reserve(distributedIds.size()*(dim_+1)); // estimate, lower bound

    // loop over the compressed distributed entities
    for (int distributedId : distributedIds){

      // temp for the offset of this id
      int const thisOffset = offsets[distributedId];

      // loop over the references and map them
      for (int j=0; j<counts[distributedId]; ++j)

        // push the mapped reference
        result.push_back(gidToFlatId.at(in[thisOffset+j]));

    }
  }




  /*!
    @brief Merge post distribution data where the data is variable length lists
           for each entity.

    @param[in] in   The post distribution data to be merged
    @param[in] counts  The vector of the number of references associated with
                       this entity, e.g. number of faces for this cell
    @param[in] distributedIds  The map from gid to first occurrence pre distribution of
                         the "from" part of the mapping, e.g. cell in cellToFace
    @return The new vector of merged data

    This signature is called after distribution. It takes post distribution lists of data
    and corrects the "from" part of the map to remove duplicates. It leaves the
    data values untouched unlike the signature above which updates topological references.
   */
  template<class T>
  void merge_duplicate_lists(std::vector<T>const& in, std::vector<int> const & counts,
    std::vector<int>const& distributedIds, std::vector<T>& result){

    // allocate offses
    std::vector<int> offsets(counts.size());

    // compute offsets (note the first element is zero and correct)
    std::partial_sum(counts.begin(), counts.end()-1, offsets.begin()+1);

    // make sure the result is clear and approximately sized
    result.clear();
    result.reserve(distributedIds.size()*(dim_+1)); // estimate, lower bound

    // loop over the compressed distributed entities
    for (int distributedId : distributedIds){

      // temp for the offset of this id
      int const thisOffset = offsets[distributedId];

      // loop over the references and map them
      for (int j=0; j<counts[distributedId]; ++j)

        // push the mapped reference
        result.push_back(in[thisOffset+j]);

    }
  }



}; // MPI_Bounding_Boxes

} // namespace Portage

#endif  // WONTON_ENABLE_MPI

#endif // MPI_Bounding_Boxes_H_
