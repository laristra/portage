/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MPI_BOUNDING_BOXES_H_
#define MPI_BOUNDING_BOXES_H_

//#define DEBUG_MPI

#include <cassert>
#include <algorithm>
#include <numeric>
#include <memory>
#include <unordered_map>
#include <map>
#include <vector>

#include "portage/support/portage.h"
#include "wonton/support/Point.h"
#include "wonton/state/state_vector_uni.h"
#include "mpi.h"
using Wonton::Point;
/*!
  @file mpi_bounding_boxes.h
  @brief Distributes source data using MPI based on bounding boxes
 */

namespace Portage {


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
  MPI_Bounding_Boxes() {}


  /*!
    @brief Destructor of MPI_Bounding_Boxes
   */
  ~MPI_Bounding_Boxes() {}


  /*!
    @brief Helper structure containg comms info for a given entity type
   */
  struct comm_info_t {
    //!< Number of total/owned entities in source field
    int sourceNum = 0, sourceNumOwned = 0;
    //!< Number of total/owned entities in new field
    int newNum = 0, newNumOwned = 0;
    //! Array of total/owned send sizes from me to all PEs
    std::vector<int> sendCounts, sendOwnedCounts;
    //! Array of total/owned recv sizes to me from all PEs
    std::vector<int> recvCounts, recvOwnedCounts;
  };


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
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    int dim = dim_ = source_mesh_flat.space_dimension();
    assert(dim == target_mesh.space_dimension());

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
    std::vector<int>& sourceCellGlobalIds = source_mesh_flat.get_global_cell_ids();
    std::vector<int> distributedCellGlobalIds(cellInfo.newNum);
    sendField(cellInfo, commRank, commSize, MPI_INT, 1,
              sourceCellGlobalIds, &distributedCellGlobalIds);
              
    // SEND GLOBAL NODE IDS
    std::vector<int>& sourceNodeGlobalIds = source_mesh_flat.get_global_node_ids();
    std::vector<int> distributedNodeGlobalIds(nodeInfo.newNum);
    sendField(nodeInfo, commRank, commSize, MPI_INT, 1,
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
    merge_data(distributedCoords, distributedNodeIds_, sourceCoords, dim_);
    


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
      merge_data(distributedCellNodeCounts, distributedCellIds_, sourceCellNodeCounts);

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
      std::vector<int> distributedCellToNodeList(cellToNodeInfo.newNum);
      sendField(cellToNodeInfo, commRank, commSize, MPI_INT, 1,
                to_gid(sourceCellToNodeList, sourceNodeGlobalIds), &distributedCellToNodeList);


      // merge and map cell node lists
      merge_lists(distributedCellToNodeList, distributedCellNodeCounts, 
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
      std::vector<int>& sourceFaceGlobalIds = source_mesh_flat.get_global_face_ids();
      std::vector<int> distributedFaceGlobalIds(faceInfo.newNum);
      sendField(faceInfo, commRank, commSize, MPI_INT, 1,
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
      merge_data( distributedCellFaceCounts, distributedCellIds_, sourceCellFaceCounts);
      
      // SEND CELL-TO-FACE MAP
      // map the cell face list vector to gid's
      std::vector<int> sourceCellToFaceList__= to_gid(sourceCellToFaceList, sourceFaceGlobalIds);
      
      // For this array only, pack up face IDs + dirs and send together
      std::vector<bool>& sourceCellToFaceDirs = source_mesh_flat.get_cell_to_face_dirs();
      for (unsigned int j=0; j<sourceCellToFaceList.size(); ++j)
      {
        int f = sourceCellToFaceList__[j];
        int dir = static_cast<int>(sourceCellToFaceDirs[j]);
        sourceCellToFaceList[j] = (f << 1) | dir;
      }

      std::vector<int> distributedCellToFaceList(cellToFaceInfo.newNum);
      sendField(cellToFaceInfo, commRank, commSize, MPI_INT, 1,
                sourceCellToFaceList, &distributedCellToFaceList);

      // Unpack face IDs and dirs
      std::vector<bool> distributedCellToFaceDirs(cellToFaceInfo.newNum);
      for (unsigned int j=0; j<distributedCellToFaceList.size(); ++j)
      {
        int fd = distributedCellToFaceList[j];
        distributedCellToFaceList[j] = fd >> 1;
        distributedCellToFaceDirs[j] = fd & 1;
      }

      
      // merge and map cell face lists
      merge_lists(distributedCellToFaceList, distributedCellFaceCounts, 
        distributedCellIds_, gidToFlatFaceId_, sourceCellToFaceList);
        
      // merge cell face directions
      merge_lists(distributedCellToFaceDirs, distributedCellFaceCounts, 
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
      std::vector<int> distributedFaceToNodeList(faceToNodeInfo.newNum);
      sendField(faceToNodeInfo, commRank, commSize, MPI_INT, 1,
                to_gid(sourceFaceToNodeList, sourceNodeGlobalIds), &distributedFaceToNodeList);      
                
      // merge and set face node counts
      merge_data( distributedFaceNodeCounts, distributedFaceIds_, sourceFaceNodeCounts);
      
      // merge and map face node lists
      merge_lists(distributedFaceToNodeList, distributedFaceNodeCounts, 
        distributedFaceIds_, gidToFlatNodeId_, sourceFaceToNodeList);

      // merge face global ids
      merge_data(distributedFaceGlobalIds, distributedFaceIds_,sourceFaceGlobalIds);
      
      // set counts for faces in the flat mesh
      source_mesh_flat.set_num_owned_faces(flatFaceNumOwned_);
                                 
    }
    
    // need to do this at the end of the mesh stuff, because converting to 
    // gid uses the global id's and we don't want to modify them before we are
    // done converting the old relationships
    // merge global ids and set in the flat mesh    
    merge_data(distributedCellGlobalIds, distributedCellIds_, sourceCellGlobalIds);
    merge_data(distributedNodeGlobalIds, distributedNodeIds_, sourceNodeGlobalIds);

    // set counts for cells and nodes in the flat mesh
    source_mesh_flat.set_num_owned_cells(flatCellNumOwned_);
    source_mesh_flat.set_num_owned_nodes(flatNodeNumOwned_);

    // Finish initialization using redistributed data
    source_mesh_flat.finish_init();


    // SEND FIELD VALUES
    
    // multimaterial state info
    int nmats = source_state_flat.num_materials();
    comm_info_t num_mats_info, num_mat_cells_info;
    std::vector<int> all_material_ids, all_material_shapes, all_material_cells;

/*

    // Is the a multimaterial problem? If so we need to pass the cell indices
    // in addition to the field values
    if (nmats>0){
    
      std::cout << "in distribute, this a multimaterial problem with " << nmats << " materials\n";
      
      /////////////////////////////////////////////////////////
      // get the material ids across all nodes
      /////////////////////////////////////////////////////////
      
      // set the info for the number of materials on each node
      setSendRecvCounts(&num_mats_info, commSize, sendFlags, nmats, nmats);
      
      // get the sorted material ids on this node
      std::vector<int> material_ids=source_state_flat.get_material_ids();
      
      // resize the post distribute material id vector
      all_material_ids.resize(num_mats_info.newNum);

      // sendData
      sendData(commRank, commSize, MPI_INT, 1, 0, num_mats_info.sourceNum, 0,
        num_mats_info.sendCounts, num_mats_info.recvCounts,
        material_ids, &all_material_ids
      );
      
      /////////////////////////////////////////////////////////
      // get the material cell shapes across all nodes
      /////////////////////////////////////////////////////////      

      // get the sorted material shapes on this node
      std::vector<int> material_shapes=source_state_flat.get_material_shapes();
      
      // resize the post distribute material id vector
      all_material_shapes.resize(num_mats_info.newNum);

      // sendData
      sendData(commRank, commSize, MPI_INT, 1, 0, num_mats_info.sourceNum, 0,
        num_mats_info.sendCounts, num_mats_info.recvCounts,
        material_shapes, &all_material_shapes
      );
     
      /////////////////////////////////////////////////////////
      // get the material cell ids across all nodes
      /////////////////////////////////////////////////////////
      
      // get the total number of material cell id's on this node
      int nmatcells = source_state_flat.num_material_cells();
      
      // set the info for the number of materials on each node
      setSendRecvCounts(&num_mat_cells_info, commSize, sendFlags, nmatcells, nmatcells);
      
      // get the sorted material ids on this node
      std::vector<int> material_cells=source_state_flat.get_material_cells();
      
      // resize the post distribute material id vector
      all_material_cells.resize(num_mat_cells_info.newNum);

      // sendData
      sendData(commRank, commSize, MPI_INT, 1, 0, num_mat_cells_info.sourceNum, 0,
        num_mat_cells_info.sendCounts, num_mat_cells_info.recvCounts,
        material_cells, &all_material_cells
      );
      
      /////////////////////////////////////////////////////////
      // compute the cell map from local id to global id back to
      // first appearance of the local id. This code was copied
      // verbatim from flat_mesh_wrapper.h
      /////////////////////////////////////////////////////////

      // Global to local maps for cells and nodes
      std::map<int, int> globalCellMap;
      std::vector<int> cellUniqueRep(distributedCellGlobalIds.size());
      for (unsigned int i=0; i<distributedCellGlobalIds.size(); ++i) {
        auto itr = globalCellMap.find(distributedCellGlobalIds[i]);
        if (itr == globalCellMap.end()) {
          globalCellMap[distributedCellGlobalIds[i]] = i;
          cellUniqueRep[i] = i;
        }
        else {
          cellUniqueRep[i] = itr->second;
        }
      }

      /////////////////////////////////////////////////////////
      // Fix the material cell indices to account for the fact
      // that they are local indices on the different nodes and 
      // need to be consistent within this flat mesh. FixListIndices
      // below is hard coded to work with ghost cells which we don't have
      // so I'm going to do this inline.
      /////////////////////////////////////////////////////////
      
      // get the local cell index offsets
      std::vector<int> offsets(commSize);
      offsets[0]=0;

      std::partial_sum(cellInfo.recvCounts.begin(),cellInfo.recvCounts.end()-1,
          offsets.begin()+1);     
          
      // make life easy by keeping a running counter
      int running_counter=0;

      // loop over the ranks
      for (int rank=0; rank<commSize; ++rank){
          
          // get the cell offset for this rank
          int this_offset = offsets[rank];
          
          // get the number of material cells for this rank
          int nmat_cells = num_mat_cells_info.recvCounts[rank];
          
          // loop over material cells on this rank
          for (int i=0; i<nmat_cells; ++i){
              
              // offset the indices in all_material_cells
              all_material_cells[running_counter] += this_offset;
              
              // use cellUniqueRep to map each cell down to its first appearance
              all_material_cells[running_counter++] = cellUniqueRep[all_material_cells[running_counter]];
          }    
      }
      
      /////////////////////////////////////////////////////////
      // We need to turn the flattened material cells into a correctly shaped
      // ragged right structure for use as the material cells in the flat
      // state wrapper. Just as in the flat state mesh field (and associated
      // cell ids), we aren't removing duplicates, just concatnating by material
      /////////////////////////////////////////////////////////
      
      // allocate the material indices
      std::unordered_map<int,std::vector<int>> material_indices;
      
      // reset the running counter
      running_counter=0;
      
      // loop over material ids on different nodes
      for (int i=0; i<all_material_ids.size(); ++i){
      
          // get the current working material
          int mat_id = all_material_ids[i];
          
          // get the current number of material cells for this material
          int nmat_cells = all_material_shapes[i];
          
          // get or create a reference to the correct material cell vector
          std::vector<int>& these_material_cells = material_indices[mat_id];
          
          // loop over the correct number of material cells
          for (int j=0; j<nmat_cells; ++j){
              these_material_cells.push_back(all_material_cells[running_counter++]);
          }
          
      }
      
      // We are reusing the material cells and cell materials. Since we are using
      // maps and sets we want to make sure that we are starting clean and not
      // adding to cruft that is already there.
      source_state_flat.clear_material_cells();
      
      // add the material indices by keys
      for ( auto& kv: material_indices){
          source_state_flat.mat_add_cells(kv.first, kv.second);
      }
      
    }
*/
    // Send and receive each field to be remapped
    for (std::string field_name : source_state_flat.names())
    {
      // this is a packed version of the field with copied field values and is
      // not a pointer to the original field
      std::vector<double> sourceField = source_state_flat.pack(field_name);
      
      // get the field stride
      int sourceFieldStride = source_state_flat.get_field_stride(field_name);

      // Currently only cell and node fields are supported
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
      // is still packed and raw doubles and will need to be unpacked
      std::vector<double> distributedField(sourceFieldStride*info.newNum);


      // send the field
      sendField(info, commRank, commSize, MPI_DOUBLE, sourceFieldStride,
        sourceField, &distributedField);
      
      // merge the data before unpacking
      std::vector<double> tempDistributedField;
      if (source_state_flat.get_entity(field_name) == Entity_kind::NODE){
      
        // node mesh field
        merge_data(distributedField, distributedNodeIds_, tempDistributedField, sourceFieldStride);
        
      } else if (source_state_flat.field_type(Entity_kind::CELL, field_name) == Wonton::Field_type::MESH_FIELD){
      
        // mesh cell field
        merge_data(distributedField, distributedCellIds_, tempDistributedField, sourceFieldStride);

      } else {
        // multi material field             
      }
         
      // unpack the field, has the correct data types, but still is not merged
      source_state_flat.unpack(field_name, tempDistributedField, all_material_ids, all_material_shapes);
          
    } 

  } // distribute

  private:
 
  int dim_;
  
  // the number of nodes "owned" by the flat mesh. "Owned" is in quotes because
  // a node may be "owned" by multiple partitions in the flat mesh. A node is
  // owned by the flat mesh if it was owned by any partition.
  int flatNodeNumOwned_;
  
  // the global id's of the kept nodes in the flat mesh and their indices in the
  // distributed node global id's
  std::vector<int> flatNodeGlobalIds_, distributedNodeIds_;
  
  // maps from gid to distributed node index and flat node index  
  std::map<int,int> gidToFlatNodeId_;

  // the number of faces "owned" by the flat mesh. "Owned" is in quotes because
  // a face may be "owned" by multiple partitions in the flat mesh. A face is
  // owned by the flat mesh if it was owned by any partition.
  int flatFaceNumOwned_;
  
  // the global id's of the kept faces in the flat mesh and their indices in the
  // distributed face global id's
  std::vector<int> flatFaceGlobalIds_, distributedFaceIds_;
  
  // maps from gid to distributed face index and flat face index  
  std::map<int,int> gidToFlatFaceId_;  
  
  // the number of cells "owned" by the flat mesh. "Owned" is in quotes because
  // a cell may be "owned" by multiple partitions in the flat mesh. A cell is
  // owned by the flat mesh if it was owned by any partition.
  int flatCellNumOwned_;
  
  // the global id's of the kept cells in the flat mesh and their indices in the
  // distributed cell global id's
  std::vector<int> flatCellGlobalIds_, distributedCellIds_;
  
  // maps from gid to distributed cell index and flat cell index  
  std::map<int,int> gidToFlatCellId_;
  
  
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
    for (unsigned int i=0; i<commSize; i++)
      info->sendCounts[i] = sendFlags[i] ? info->sourceNum : 0;
    MPI_Alltoall(&(info->sendCounts[0]), 1, MPI_INT,
                 &(info->recvCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Each rank will tell each other rank how many owned indexes it is going to send it
    info->sendOwnedCounts.resize(commSize);
    info->recvOwnedCounts.resize(commSize);
    for (unsigned int i=0; i<commSize; i++)
      info->sendOwnedCounts[i] = sendFlags[i] ? info->sourceNumOwned : 0;
    MPI_Alltoall(&(info->sendOwnedCounts[0]), 1, MPI_INT,
                 &(info->recvOwnedCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of indexes this rank will receive from all ranks
    for (unsigned int i=0; i<commSize; i++)
      info->newNum += info->recvCounts[i];
    for (unsigned int i=0; i<commSize; i++)
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
                 const int commRank, const int commSize,
                 const MPI_Datatype mpiType, const int stride,
                 const std::vector<T>& sourceData,
                 std::vector<T>* newData)
  {
    // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
    std::vector<int> recvGhostCounts(commSize);
    std::vector<int> sendGhostCounts(commSize);

    for (unsigned int i=0; i<commSize; i++)
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
    for (unsigned int f=0; f<(*newData).size(); f++)
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
  void sendData(const int commRank, const int commSize,
                const MPI_Datatype mpiType, const int stride,
                const int sourceStart, const int sourceEnd,
                const int newStart,
                const std::vector<int>& curSendCounts,
                const std::vector<int>& curRecvCounts,
                const std::vector<T>& sourceData,
                std::vector<T>* newData)
  {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = newStart;
    int myOffset;
    std::vector<MPI_Request> requests;
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curRecvCounts[i] > 0))
      {
        MPI_Request request;
        MPI_Irecv((void *)&((*newData)[stride*writeOffset]),
                  stride*curRecvCounts[i], mpiType, i,
                  MPI_ANY_TAG, MPI_COMM_WORLD, &request);
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
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curSendCounts[i] > 0))
      {
        MPI_Send((void *)&(sourceData[stride*sourceStart]),
                 stride*curSendCounts[i], mpiType, i, 0, MPI_COMM_WORLD);
      }
    }

    // Wait for all receives to complete
    if (requests.size() > 0)
    {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }
  } // sendData


  template <class Source_Mesh, class Target_Mesh>        
  void compute_sendflags(Source_Mesh & source_mesh_flat, Target_Mesh &target_mesh, 
              std::vector<bool> &sendFlags){
                
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    int dim = source_mesh_flat.space_dimension();

    // Compute the bounding box for the target mesh on this rank, and put it in an
    // array that will later store the target bounding box for each rank
    int targetNumOwnedCells = target_mesh.num_owned_cells();
    std::vector<double> targetBoundingBoxes(2*dim*commSize);
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<targetNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      target_mesh.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (unsigned int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        // ugly hack, since dim is not known at compile time
        if (dim == 3)
        {
          Point<3> nodeCoord;
          target_mesh.node_get_coordinates(n, &nodeCoord);
          for (unsigned int k=0; k<dim; ++k)
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
          for (unsigned int k=0; k<dim; ++k)
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
    int sourceNumOwnedCells = source_mesh_flat.num_owned_cells();
    std::vector<double>& sourceCoords = source_mesh_flat.get_coords();

    // Compute the bounding box for the source mesh on this rank
    std::vector<double> sourceBoundingBox(2*dim);
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      sourceBoundingBox[i+0] = std::numeric_limits<double>::max();
      sourceBoundingBox[i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<sourceNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      source_mesh_flat.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (unsigned int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        for (unsigned int k=0; k<dim; ++k)
        {
          if (sourceCoords[n*dim+k] < sourceBoundingBox[2*k])
            sourceBoundingBox[2*k] = sourceCoords[n*dim+k];
          if (sourceCoords[n*dim+k] > sourceBoundingBox[2*k+1])
            sourceBoundingBox[2*k+1] = sourceCoords[n*dim+k];
        }
      } // for j
    } // for c

    // Broadcast the target bounding boxes so that each rank knows the bounding boxes for all ranks
    for (unsigned int i=0; i<commSize; i++)
      MPI_Bcast(&(targetBoundingBoxes[0])+i*2*dim, 2*dim, MPI_DOUBLE, i, MPI_COMM_WORLD);


    // Offset the source bounding boxes by a fudge factor so that we don't send source data to a rank
    // with a target bounding box that is incident to but not overlapping the source bounding box
    const double boxOffset = 2.0*std::numeric_limits<double>::epsilon();
    double min2[dim], max2[dim];
    for (unsigned int k=0; k<dim; ++k)
    {
      min2[k] = sourceBoundingBox[2*k]+boxOffset;
      max2[k] = sourceBoundingBox[2*k+1]-boxOffset;
    }

    // For each target rank with a bounding box that overlaps the bounding box for this rank's partition
    // of the source mesh, we will send it all our source cells; otherwise, we will send it nothing
    for (unsigned int i=0; i<commSize; ++i)
    {
      double min1[dim], max1[dim];
      bool sendThis = true;
      for (unsigned int k=0; k<dim; ++k)
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
  std::vector<int> to_gid(std::vector<int> const& in, vector<int>const& gids){
    std::vector<int> result;
    result.reserve(in.size());
    for (auto x:in) result.push_back(gids[x]);
    return result;
  }


  /*!
    @brief Compress distributed global id's into the vector of distributed
    indices(first occurrence for each global id) and the vector of global id's.
    The global id's are in ascending order of "owned" global id, followed by "ghost"
    global id. In the flat mesh, an entitiy is considered "owned" if it is 
    owned by any partition, and a "ghost" if it is not owned by any partition
           
    @param[in] distributedGlobalIds  The vector of gid's for each entity post distribution
    @param[in] distributedNumOwned  The number of entities that are owned after distribution
    @param[out] distributedIds  The vector of indices of first occurrence in 
      distributedGlobalIds that produce the global id's in the flat mesh
    @param[out] flatGlobalIds  The global ids in the flat mesh. They are in order
      of ascending gid for owned entities followed by ascending gid for ghosts
    @param[out] flatNumOwned  The number of owned entities in the flat mesh
    
    This is called after distribution of gids for each entity kind. The vector of
    gids (with potentially duplicated gid's) is used to create two vectors. The first
    vector is the first position of occurrence in post distribution gids vector. 
    The second vector is the gid in the flat mesh after compression. 
    
    We handle"ghosts" are in the post distributed data. An entity is considered
    "owned" by the flat mesh if it was owned on any partition and a "ghost" if
    was a ghost on all partitions. This implies that an entity will be owned by as many
    partitions as the original source partition was sent to. So while we preserve
    the concept that a ghost is a ghost in the flat mesh only if it was originally
    ghosts, we unfortunately have that an entity in the flat mesh will typically 
    be owned by many different paritions
  */
  void compress_with_ghosts(std::vector<int> const& distributedGlobalIds, 
    int const distributedNumOwned, std::vector<int>& distributedIds, 
    std::vector<int>& flatGlobalIds, int &flatNumOwned){
    
    // a map for keeping track of what we have already seen
    std::map<int,int> uniqueGid, uniqueGhostGid;
    
    // loop over owned cells in the distributed global id's
    for (int i=0 ; i<distributedNumOwned ; ++i){
    
      // get the current gid
      int gid = distributedGlobalIds[i];
      
      // is this gid new
      if (uniqueGid.find(gid)==uniqueGid.end()){
      
        // make the gid seen
        uniqueGid[gid]=i;      
        
      }
      
    }
    
    // We have processed owned cells in the distributed mesh, so everything
    // we have collected to this point is considered owned
    flatNumOwned = uniqueGid.size();
    
    // push the owned cells first and in gid order
    for (auto const& kv : uniqueGid){
    
      // push to the flat cells gid
      flatGlobalIds.push_back(kv.first);
      
      // push to the distributed cells id
      distributedIds.push_back(kv.second);
    
    }
    
    // loop over owned cells in the distributed global id's
    for (int i=distributedNumOwned ; i<distributedGlobalIds.size() ; ++i){
    
      // get the current gid
      int gid = distributedGlobalIds[i];
      
      // is this gid new
      if (uniqueGid.find(gid)==uniqueGid.end() &&
        uniqueGhostGid.find(gid)==uniqueGhostGid.end()){
      
        // make the gid seen
        uniqueGhostGid[gid]=i;      
        
      }
      
    }
        
    // push the owned cells first and in gid order
    for (auto const& kv : uniqueGhostGid){
    
      // push to the flat cells gid
      flatGlobalIds.push_back(kv.first);
      
      // push to the distributed cells id
      distributedIds.push_back(kv.second);
    
    }
    
  }


  /*!
    @brief Create a map from flat mesh gid to flat mesh index
           
    @param[in] flatGlobalIds  The vector of gid's in the flat mesh
    @param[out] gidToFlat  The map from global id to flat cell index

    This function creates a trivial map from global id to flat cell index.
    The index in the global id is the value of the map. We use this later when
    converting global id's into the new flat mesh local index.
  */
  void create_gid_to_flat_map(std::vector<int> const& flatGlobalIds, 
    std::map<int,int>& gidToFlat){
    
    for (int i = 0; i < flatGlobalIds.size(); ++i) 
      gidToFlat[flatGlobalIds[i]]=i;
    
  }


  /*!
    @brief Merge post distribution flattend data so that each global id appears 
           only once.
           
    @param[in] in   The post distribution data to be merged
    @param[in] distributedIds  The vector of distributed indices to keep
    @param[out] result The new vector of flattened and merged data
    @param[in] stride=1  The number of values associated with the each element of data.
                       Scalar data has stride 1, centroid data has stride D. 
    
    This is called after distribution. It takes post distribution data and 
    selects a single occurrence of each gid defined by distributedIds
   */
  template<class T>
  void merge_data(std::vector<T>const& in, std::vector<int>const& distributedIds, 
    std::vector<T>& result, int const stride=1){

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
    reference data and corrects for two things. The first is the input has
    duplicate offset lists that need to be merged. The second is the references
    need to get converted from gid to new flat index id.   
  */
  void merge_lists(std::vector<int>const& in, std::vector<int> const& counts,
    std::vector<int>const& distributedIds, std::map<int,int>const& gidToFlatId,
    std::vector<int>& result){
    
    // allocate offsets
    std::vector<int> offsets(counts.size());
    
    // compute offsets (note the first element is zero and correct)
    std::partial_sum(counts.begin(), counts.end()-1, offsets.begin()+1);
    
    // make sure the result is clear and approximately sized
    result.clear();
    result.reserve(distributedIds.size()*(dim_+1)); // estimate, lower bound
    
    // loop over the compressed distributed entities
    for (int i=0; i<distributedIds.size(); ++i){
    
      // temp for the offset of this id
      int const thisOffet = offsets[distributedIds[i]];
    
      // loop over the references and map them
      for (int j=0; j<counts[i]; ++j)
      
        // push the mapped reference 
        result.push_back(gidToFlatId.at(in[thisOffet+j]));   
             
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
    
    This is called after distribution. It takes post distribution lists of data 
    and corrects the "from" part of the map to remove duplicates. The only usage
    of this signature is compressing the face direction which is a boolean. There
    are lists for each cell in distributedCellToFaceDirs, but the values are kept
    unchanged unlike the signature above with not template parameter where the
    values are references that need to be remapped as well. Merge_data won't 
    work because there are variable numbers of data for each entity.
   */
  template<class T>
  void merge_lists(std::vector<T>const& in, std::vector<int> const & counts,
    std::vector<int>const& distributedIds, std::vector<T>& result){
    
    // allocate offses
    std::vector<int> offsets(counts.size());
    
    // compute offsets (note the first element is zero and correct)
    std::partial_sum(counts.begin(), counts.end()-1, offsets.begin()+1);
    
    // make sure the result is clear and approximately sized
    result.clear();
    result.reserve(distributedIds.size()*(dim_+1)); // estimate, lower bound

    // loop over the compressed distributed entities
    for (int i=0; i<distributedIds.size(); ++i){
    
      // temp for the offset of this id
      int const thisOffet = offsets[distributedIds[i]];
    
      // loop over the references and map them
      for (int j=0; j<counts[i]; ++j)
      
        // push the mapped reference 
        result.push_back(in[thisOffet+j]);   
             
    }    
  }
  
  

}; // MPI_Bounding_Boxes

} // namespace Portage

#endif // MPI_Bounding_Boxes_H_
