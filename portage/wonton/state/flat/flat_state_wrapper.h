/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/




#ifndef FLAT_STATE_WRAPPER_H_
#define FLAT_STATE_WRAPPER_H_

#include <map>
#include <memory>
#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <set>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

//#define DBG

/*!
  @file flat_state_wrapper.h
  @brief Wrapper for interfacing with the Flat state manager
 */

namespace Wonton {

/*!
  @class Flat_State_Wrapper "flat_state_wrapper.h"
  @brief Stores state data in a flat representation

         Currently all fields must be of the same type
*/
template <class T=double>
class Flat_State_Wrapper {
 public:

  //! pair of name and entity to be used as data key
  using pair_t = std::pair<std::string, Entity_kind>;

  /*!
    @brief Constructor of Flat_State_Wrapper
   */
  Flat_State_Wrapper() { };

  /*!
   * @brief Initialize the state wrapper with another state wrapper and a list of names
   * @param[in] input another state wrapper, which need not be for Flat_State.
   * @param[in] var_names a list of state names to initialize
   *
   * Entities and sizes associated with the given name will be obtained from the input state wrapper.
   *
   * A name can be re-used with a different entity, but a name-entity combination
   * must be unique.
   *
   * A name-entity combination must not introduce a new size for that entity if
   * it has previously been encountered.
   *
   * All existing internal data is forgotten.
   */
  template <class State_Wrapper>
  void initialize(State_Wrapper const & input,
                  std::vector<std::string> var_names)
  {
		clear(); // forget everything
		
		// keep track of whether we have already done multimaterial setup
		bool has_multimaterial_data = false;
		
		// loop over field names
    for (unsigned int i=0; i<var_names.size(); i++)
    {
    
      // get field name
      std::string varname = var_names[i];
      
				  #ifdef DBG
      std::cout << "processing " << varname << std::endl;
					#endif

      // get entity kind
      Entity_kind entity = input.get_entity(varname);

			// define the raw data pointer
			// in a better world we wouldn't need to copy to this intermediate storage
			T const* data;
			
			// store the field types into the private data member
			field_types_[pair_t(varname,entity)]=input.field_type(entity, varname);
			
			// determine whether this field is single field data or multimaterial data
			if (input.field_type(entity, varname)==Portage::Field_type::MULTIMATERIAL_FIELD){
			
				// the field is multimaterial
				  #ifdef DBG
				std::cout << varname << " is a multimaterial field"<<std::endl;
					#endif
								
				// all multimaterial fields are defined on the same cells. If the
				// material exist in the cell, then all multimaterial fields need to be
				// defined. The following block of code which defines the shape of the
				// data only needs to be done once, regardless of the number of multi
				// material fields
				if (!has_multimaterial_data){
				
				  #ifdef DBG
					std::cout << "hopefully single setting of sizes " <<std::endl;
					#endif
				
					// get the number of materials and set the private member
					nmats_ = input.num_materials();
				  #ifdef DBG
					std::cout << "  the state wrapper has " << nmats_ << " materials" <<std::endl;
					#endif
					
					// allocate space so we don't need to dynamically allocate
					mat_names_.reserve(nmats_);
					mat_ncells_.reserve(nmats_);
					mat_cellids_.reserve(nmats_);
					
					// loop over materials
					for (int matid=0; matid < nmats_; ++matid){
					
						// set the material names
						mat_names_.push_back(input.material_name(matid));
					
					  // set the number of cells for this material 
					  mat_ncells_.push_back(input.mat_get_num_cells(matid));
					  
					  // set the cell id's for this material
					  std::vector<int> cellids;
					  input.mat_get_cells(matid, &cellids);
					  mat_cellids_.push_back(cellids);
					  
					  // update the list of materials for each cell id
					  for (int cellid: cellids){
					  	cell_to_mats_[cellid].insert(matid);
					  }
					  
					  // increase the private member, flat_size
						flat_size += mat_ncells_[matid];
						
					}

					// make sure this setup code block is run only once
					has_multimaterial_data = true;
					
				  #ifdef DBG
					std::cout << "  flat_size= " << flat_size<<std::endl;
						#endif
				
				}
							
				// concatenate the data vectors
				std::vector<T> vdata;

				// loop over materials
				for (int matid=0; matid < nmats_; ++matid){
								  
				  // copy this material field data into the temp array
					input.mat_get_celldata(varname, matid, &data);
					
				  // append state data to the flat storage vector
					vdata.insert(vdata.end(), data, data+mat_ncells_[matid]);
					
				}
				
				// finished copying into vdata. Let's check what's there
				  #ifdef DBG
				for (int j=0; j<flat_size; ++j) {std::cout << vdata[j] << " "<<std::endl;}
					#endif
				
				// add to database (by flattening, we are able to shoehorn the multimaterial
				// data into the same format as regular field data)
		    mesh_add_data(entity, varname, std::make_shared<std::vector<T>>(vdata));

			} else {
			
		    // get pointer to data for state from input state wrapper
		    input.mesh_get_data(entity, varname, &data);

		    // copy input state data into new vector storage
		    size_t dataSize = input.get_data_size(entity, varname);
		    auto vdata = std::make_shared<std::vector<T>>(dataSize);
				std::copy(data, data+dataSize, vdata->begin());

				// add to database
		    mesh_add_data(entity, varname, vdata);
      }
    }
  }

  /*!
    @brief Assignment operator (disabled) - don't know how to implement (RVG)
   */
  Flat_State_Wrapper & operator=(Flat_State_Wrapper const &) = delete;

  /*!
    @brief Empty destructor
   */
  ~Flat_State_Wrapper() {};

  /*!
   * @brief Initialize the state wrapper with explicit lists of names, entities and data
   * @param[in] names a list of state names to initialize
   * @param[in] entities entities corresponding to the names
   * @param[in] data the state vectors to be stored
   *
   * A name can be re-used with a different entity, but a name-entity combination
   * must be unique.
   *
   * A name-entity combination must not introduce a new size for that entity if
   * the entity has previously been encountered.
   *
   * All existing internal data is forgotten.
   */
  void initialize(std::vector<std::string> const& names,
                  std::vector<Entity_kind> const& entities,
                  std::vector<std::shared_ptr<std::vector<T>>> const& data)
  {
    if (not (names.size() == entities.size() and names.size() == data.size() and data.size() == entities.size())) {
        throw std::runtime_error("argument sizes do not agree");
    }

    clear();

    size_t index;
    for (size_t i=0; i<names.size(); i++) {
        mesh_add_data(entities[i], names[i], data[i]);
    }
  }

  /*!
   @brief Get names of data fields associated with a given entity
   @param[in] entity The entity type
   @param[out] names The names associated with the entity. Cleared on entry.
  */
  void get_names(Entity_kind on_what, std::vector<std::string>& names) {
    names.clear();
    for (auto iter=entity_map_.begin(); iter!=entity_map_.end(); iter++) {
      if (iter->second == on_what) {
        names.push_back(iter->first);
      }
    }
  }

  /*!
    @brief Number of materials in problem
  */

  int num_materials() const {
    return nmats_;
  }

  /*!
    @brief Name of material
  */

  std::string material_name(int matid) const {
    return mat_names_[matid];
  }

  /*!
    @brief Get number of cells containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @return         Number of cells containing material 'matid'
  */

  int mat_get_num_cells(int matid) const {
    return mat_ncells_[matid];
  }

  /*!
    @brief Get cell indices containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @param matcells Cells containing material 'matid'
  */

  void mat_get_cells(int matid, std::vector<int> *matcells) const {
    matcells->clear();
    *matcells = std::vector<int>(mat_cellids_[matid]);
  }

  /*!
    @brief Get number of materials contained in a cell
    @param cellid  Index of cell in mesh
    @return        Number of materials in cell
  */

  int cell_get_num_mats(int cellid) const {
		return cell_to_mats_.at(cellid).size();
  }

  /*!
    @brief Get the IDs of materials in a cell
    @param cellid    Index of cell in mesh
    @param cellmats  Indices of materials in cell
  */

  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
    // unfortunately, I don't of a way to convert from a set to a vector
    // without an explicit loop or transform
    for (auto x: cell_to_mats_.at(cellid)) cellmats->push_back(x);
  }

  /*!
    @brief Get the local index of mesh cell in material cell list
    @param meshcell    Mesh cell ID
    @param matid       Material ID
    @return             Local cell index in material cell list
  */

  int cell_index_in_material(int meshcell, int matid) const {
  	const std::vector<int> cellids(mat_cellids_[matid]);
    const auto it = std::find(cellids.begin(), cellids.end(), meshcell);
    if (it!=cellids.end()){
    	return std::distance(cellids.begin(), it);
    } else {
    	return -1;
    }
  }

  /*!
    @brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
    @param onwhat    Entity_kind that field is defined on
    @param varname   Name of field
    @return          Field type
  */

  Field_type field_type(Entity_kind on_what, std::string const& var_name) const {
  	pair_t key(var_name,on_what);
    return field_types_.at(key);  
  }

  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which the data is defined
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data. Null on output if data does not exist.
    *
    * Data is associated with the name-entity combination. Both values must be valid.
   */
  void mesh_get_data(Entity_kind on_what, std::string const& var_name,
                     T** data) {
    pair_t pr(var_name, on_what);
    auto iter = name_map_.find(pr);
    if (iter != name_map_.end()) {
      (*data) = (T*)(&((*(state_[iter->second]))[0]));
    } else {
      (*data) = nullptr;
    }
  }

  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which the data is defined
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data. Null on output if data does not exist.
    *
    * Data is associated with the name-entity combination. Both values must be valid.
   */
  void mesh_get_data(Entity_kind on_what, std::string const& var_name,
                     T const **data) const {
    pair_t pr(var_name, on_what);
    auto iter = name_map_.find(pr);
    if (iter != name_map_.end()) {
      (*data) = (T const *)(&((*(state_[iter->second]))[0]));
    } else {
      (*data) = nullptr;
    }
  }

  /*!
    @brief Get pointer to read-only scalar cell data for a particular material
    @param[in] var_name The string name of the data field
    @param[in] matid   Index (not unique identifier) of the material
    @param[out] data   vector containing the values corresponding to cells in the material
   */

  void mat_get_celldata(std::string const& var_name, int matid,
                        T const **data) const {
  }


  /*!
    @brief Get pointer to read-write scalar data for a particular material
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in] matid   Index (not unique identifier) of the material
    @param[out] data   vector containing the values corresponding to cells in the material

    Removing the constness of the template parameter allows us to call
    this function and get const data back (e.g. pointer to double const)
    even if the wrapper object is not const. The alternative is to make
    another overloaded operator that is non-const but returns a pointer
    to const data. Thanks StackOverflow!
   */

  void mat_get_celldata(std::string const& var_name, int matid, T **data) {
  
  		// define the pair_t key
  		pair_t key(var_name, Portage::Entity_kind::CELL);
  		
  		// get the vector index into the state_ vector for this field
  		size_t flat_data_index= name_map_.at(key);  		
  		
  		#ifdef DBG
  		std::cout <<"flat data index for "<<var_name<<": "<<flat_data_index<<std::endl;

  		// get the flat state vector
 			std::vector<T> flat_data = *state_[flat_data_index];
 			
  		std::cout <<"flat data  for "<<var_name<<": ";
 			for (int i =0; i<flat_size; ++i){
 				std::cout<<flat_data[i]<<" ";
 			}
 			std::cout<<std::endl;
 			#endif
 			
 			int start=0;
 			for (int i=0; i<matid; ++i) start += mat_ncells_[i];
 			
  		#ifdef DBG
 			std::cout << "offset for material "<< matid << " is " << start<<std::endl;
 			#endif

			// this is a complicated expression
			// state_[flat_data_index] is the shared pointer to the vector of flattened data
			// *state_[flat_data_index]) is the flattened data
			// (*state_[flat_data_index])[start] is the start of the data for this material
			// &((*state_[flat_data_index])[start]) is a pointer to the data for this material
 			(*data) = (T*)(&((*state_[flat_data_index])[start]));
  }



  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] var_name The name of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   *
   * If the name has previously been associated with more than one entity, then the
   * entity that was most recently associated will be returned. To avoid this
   * ambiguity, please provide entity hints in the field name for yourself.
   *
   * This function is provided to make the class compatible with other state wrappers.
   */
  Entity_kind get_entity(std::string const& var_name) {
    return entity_map_[var_name];
  }

  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] index The index of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   */
  Entity_kind get_entity(int index) const
  {
    return entities_[index];
  }

  /*!
    @brief Get size for entity
  */
  size_t get_entity_size(Entity_kind ent) {
    return entity_size_map_[ent];
  }


  /*!
    @brief Get the data type of the given field
    @param[in] name The string name of the data field
    @return A reference to the type_info struct for the field's data type
   */
  const std::type_info& get_data_type(std::string const& name) const {
    size_t index = -1;
    pair_t name_cell_pair(name, Portage::Entity_kind::CELL);
    auto it = name_map_.find(name_cell_pair);
    if (it != name_map_.end())
      index = get_vector_index(Portage::Entity_kind::CELL, name);
    else {
      pair_t name_node_pair(name, Portage::Entity_kind::NODE, name);
      it = name_map_.find(name_node_pair);
      if (it != name_map_.end())
        index = get_vector_index(Portage::Entity_kind::NODE, name);
    }

    if (index >= 0) {
      return typeid(T);
    } else {
      std::cerr << "Could not find state variable " << name << "\n";
      return typeid(0);
    }
  }


  /*!
   * @brief Get index for entity and name
   */
  size_t get_vector_index(Entity_kind ent, std::string const& name) {
    pair_t pair(name, ent);
    return name_map_[pair];
  }

  /*!
    @brief Get the data vector
  */
  std::shared_ptr<std::vector<T>> get_vector(size_t index)
  {
    return state_[index];
  }

  /*!
    @brief Get gradients
  */
  std::shared_ptr<std::vector<Portage::Point3>> get_gradients(size_t index)
  {
    return gradients_[index];
  }

  /*!
    @brief Get the number of data vectors
  */
  size_t get_num_vectors() { return state_.size(); }

  /*!
    @brief Add a gradient field
  */
  void add_gradients(std::shared_ptr<std::vector<Portage::Point3>> new_grad)
  {
    if (new_grad->size() <= 0) return;
    gradients_.push_back(new_grad);
  }

  /*!
    @brief Get field stride
  */
  size_t get_field_stride(size_t index)
  {
    return 1;
  }

  /*!
    @brief Get the number of gradient vectors
  */
  size_t get_num_gradients() { return gradients_.size(); }

private:
  std::vector<std::shared_ptr<std::vector<T>>> state_;
  int nmats_ = 0;
  size_t flat_size=0;
  std::vector<int> mat_ncells_;
  std::map<pair_t, size_t> name_map_;
  std::vector<Entity_kind> entities_;
  std::map<std::string, Entity_kind> entity_map_;
  std::map<Entity_kind, size_t> entity_size_map_;
  std::vector<std::shared_ptr<std::vector<Portage::Point3>>> gradients_;
  std::vector<std::string> mat_names_;
  std::vector<std::vector<int>> mat_cellids_;
  std::unordered_map<int, std::set<int>> cell_to_mats_;
  std::map<pair_t, Portage::Field_type> field_types_;
  /* 
  In the above declarations, at some point, we should do some diagnostics to 
  determine whether maps or unordered maps are best. My initial philosophy is 
  that if the lookups are small, we might as well use an ordered map so you 
  don't have to hash, but I could be wrong. In general, the number of cells may 
  be large so I made cell_to_mats_ and unordered set while the number of fields 
  is small so best to use an ordered map.
  */
  
  /*!
   * @brief Forget all internal data so initialization can start over
   */
  void clear(){
          state_.clear();
          name_map_.clear();
          entity_map_.clear();
          entity_size_map_.clear();
          gradients_.clear();
          mat_ncells_.clear();
          mat_names_.clear();
          mat_cellids_.clear();
          cell_to_mats_.clear();
  }

  /*!
   @brief Add a scalar data field
   @param[in] on_what The entity type on which the data is defined
   @param[in] var_name The name of the data field
   @param[in] value vector of data to add
   *
   * If data for name-entity combination already exists, then replace data.
   * Size of data must match previous size recorded for requested entity.
   * If entity has not been seen before, make that entity's size equal to that of data.
   */
  void mesh_add_data(Entity_kind on_what, std::string const& var_name,
                     std::shared_ptr<std::vector<T>> data) {
    // if we have seen this entity before - check size match, else
    // store this size
    auto sziter = entity_size_map_.find(on_what);
    if (sziter != entity_size_map_.end()) {
      if (sziter->second != data->size()) {
        throw std::runtime_error(std::string("variable ")+var_name+
                                 " has incompatible size on add");
      }
    } else {
      entity_size_map_[on_what] = data->size();
    }
    
    // store data and update internal book-keeping
    pair_t pair(var_name, on_what);
    auto iter = name_map_.find(pair);
    if (iter == name_map_.end()) {  // have not seen this entity-name combo before, add data
      state_.push_back(data);
      name_map_[pair] = state_.size() - 1;
      entities_.push_back(on_what);
      entity_map_[var_name] = on_what;
      entity_size_map_[on_what] = data->size();
    } else {  // have seen the entity-name combo already, replace data. already checked size.
      std::copy(data->begin(), data->end(), state_[iter->second]->begin());
    }
  }
};  // Flat_State_Wrapper

}  // namespace Wonton

#endif  // FLAT_STATE_WRAPPER_H_
