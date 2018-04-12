/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#ifndef JALI_MMSTATE_WRAPPER_H_
#define JALI_MMSTATE_WRAPPER_H_

#include <string>
#include <vector>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

#include "Mesh.hh"       // Jali mesh declarations
#include "JaliState.h"  // Jali-based state manager declarations

/*!
  @file jali_mmstate_wrapper.h
  @brief Wrapper for interfacing with the Jali state manager
 */
namespace Wonton {

using namespace Portage;

/*!
  @class Jali_MMState_Wrapper "jali_mmstate_wrapper.h"
  @brief Provides access to data stored in Jali_State
*/
class Jali_State_Wrapper {
 public:

  /*!
    @brief Constructor of Jali_MMState_Wrapper
    @param[in] jali_state A reference to a Jali::State instance
   */
  Jali_State_Wrapper(Jali::State & jali_state) : jali_state_(jali_state) {}

  /*!
    @brief Copy constructor of Jali_State_Wrapper - not a deep copy
    @param[in] state A reference to another Jali_State_Wrapper instance
   */
  Jali_State_Wrapper(Jali_State_Wrapper & state) :
      jali_state_(state.jali_state_) {}

  /*!
    @brief Assignment operator (disabled) - don't know how to implement (RVG)
   */
  Jali_State_Wrapper & operator=(Jali_State_Wrapper const &) = delete;

  /*!
    @brief Empty destructor
   */
  ~Jali_State_Wrapper() {}


  /*!
    @brief Initialize fields from mesh file
   */
  void init_from_mesh() { jali_state_.init_from_mesh(); }

  /*!
    @brief Export fields to mesh file
   */
  void export_to_mesh() {jali_state_.export_to_mesh(); }


  /*!
    @brief Number of materials in problem
  */

  int num_materials() const {
    return jali_state_.num_materials();
  }

  /*!
    @brief Name of material
  */

  std::string material_name(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return jali_state_.material_name(matid);
  }

  /*!
    @brief Get number of cells containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @return         Number of cells containing material 'matid'
  */

  int mat_get_num_cells(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return jali_state_.material_cells(matid).size();
  }

  /*!
    @brief Get cell indices containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @param matcells Cells containing material 'matid'
  */

  void mat_get_cells(int matid, std::vector<int> *matcells) const {
    assert(matid >= 0 && matid < num_materials());
    matcells->clear();
    *matcells = jali_state_.material_cells(matid);
  }

  /*!
    @brief Get number of materials contained in a cell
    @param cellid  Index of cell in mesh
    @return        Number of materials in cell
  */

  int cell_get_num_mats(int cellid) const {
    return jali_state_.num_cell_materials(cellid);
  }

  /*!
    @brief Get the IDs of materials in a cell
    @param cellid    Index of cell in mesh
    @param cellmats  Indices of materials in cell
  */

  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
    *cellmats = jali_state_.cell_materials(cellid);
  }

  /*!
    @brief Get the local index of mesh cell in material cell list
    @param meshcell    Mesh cell ID
    @param matid       Material ID
    @return             Local cell index in material cell list
  */

  int cell_index_in_material(int meshcell, int matid) const {
    return jali_state_.cell_index_in_material(meshcell, matid);
  }

  /*!
    @brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
    @param onwhat    Entity_kind that field is defined on
    @param varname   Name of field
    @return          Field type
  */

  Field_type field_type(Entity_kind on_what, std::string const& var_name)
      const {

    Jali::State::const_iterator it = jali_state_.cbegin();
    while (it != jali_state_.cend()) {
      std::shared_ptr<Jali::StateVectorBase> bvec = *it;
      if (bvec->name() == var_name &&
          static_cast<Portage::Entity_kind>(bvec->entity_kind()) == on_what) {
        if (bvec->type() == Jali::StateVector_type::UNIVAL)
          return Field_type::MESH_FIELD;
        else
          return Field_type::MULTIMATERIAL_FIELD;
      }
      it++;
    }
  }

  
  /*!
    @brief Get a pointer to read-only single-valued data on the mesh
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A vector containing the data
   */

  template <class T>
  void mesh_get_data(Entity_kind on_what, std::string const& var_name,
                     T const **data) const {
    Jali::UniStateVector<T, Jali::Mesh> vector;
    if (jali_state_.get<T, Jali::Mesh, Jali::UniStateVector>(var_name,
                                                    jali_state_.mesh(),
                                                    (Jali::Entity_kind) on_what,
                                                    Jali::Entity_type::ALL,
                                                    &vector)) {
      *data = (T const *) vector.get_raw_data();
    } else {
      std::cerr << "Could not find state variable " << var_name << "\n";
      *data = nullptr;
    }
  }


  /*!
    @brief Get a pointer to read-write single-valued data on the mesh
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A vector containing the data

    Removing the constness of the template parameter allows us to call
    this function and get const data back (e.g. pointer to double const)
    even if the wrapper object is not const. The alternative is to make
    another overloaded operator that is non-const but returns a pointer
    to const data. Thanks StackOverflow!
   */
  template <class T>
  void mesh_get_data(Entity_kind on_what, std::string const& var_name,
                     T **data) {
    using T1 = typename std::remove_const<T>::type;
    Jali::UniStateVector<T1, Jali::Mesh> vector;
    if (jali_state_.get<T1, Jali::Mesh, Jali::UniStateVector>(var_name,
                                                    jali_state_.mesh(),
                                                    (Jali::Entity_kind) on_what,
                                                    Jali::Entity_type::ALL,
                                                    &vector)) {
      *data = vector.get_raw_data();
    } else {
      std::cerr << "Could not find state variable " << var_name << "\n";
      *data = nullptr;
    }
  }


  /*!
    @brief Get pointer to read-only scalar cell data for a particular material
    @param[in] var_name The string name of the data field
    @param[in] matid   Index (not unique identifier) of the material
    @param[out] data   vector containing the values corresponding to cells in the material
   */

  template <class T>
  void mat_get_celldata(std::string const& var_name, int matid,
                        T const **data) const {
    Jali::MultiStateVector<T, Jali::Mesh> mmvector;
    if (jali_state_.get<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                      jali_state_.mesh(),
                                                      Jali::Entity_kind::CELL,
                                                      Jali::Entity_type::ALL,
                                                      &mmvector)) {
      // data copy
      *data = (T const *) mmvector.get_raw_data(matid);
    } else {
      std::cerr << "Could not find state variable " << var_name << "\n";
      *data = nullptr;
    }
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

  template <class T>
  void mat_get_celldata(std::string const& var_name, int matid, T **data) {
    using T1 = typename std::remove_const<T>::type;

    Jali::MultiStateVector<T1, Jali::Mesh> mmvector;
    if (jali_state_.get<T1, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                      jali_state_.mesh(),
                                                      Jali::Entity_kind::CELL,
                                                      Jali::Entity_type::ALL,
                                                      &mmvector)) {
      // data copy
      *data = mmvector.get_raw_data(matid);
    } else {
      std::cerr << "Could not find state variable " << var_name << "\n";
      *data = nullptr;
    }
  }


  /*!
   @brief Add a scalar single valued data field
   @param[in] on_what The entity type on which the data is defined
   @param[in] var_name The name of the data field
   @param[in] values Initialize with this array of values
   */
  template <class T>
  void mesh_add_data(Entity_kind on_what, std::string const& var_name,
                     T const * const values) {
    jali_state_.add(var_name, jali_state_.mesh(), (Jali::Entity_kind) on_what,
                    Jali::Entity_type::ALL, values);
  }

  /*!
   @brief Add a scalar single valued data field with uniform values
   @param[in] on_what The entity type on which the data is defined
   @param[in] var_name The name of the data field
   @param[in] value Initialize with this value

   This version of the overloaded operator is being DISABLED for
   pointer and array types (via the line 'typename
   std::enable_if....type') because template deduction rules are
   making the compiler invoke this version, when we call it with a
   const double ** pointer
   
   See, stackoverflow.com Q&A
   
   http://stackoverflow.com/questions/13665574/template-argument-deduction-and-pointers-to-constants

    We could make it work for some cases using

    template <class T, class DomainType,
              template<class, class> class StateVecType>
    auto add(........,
             T const& data) -> StateVecType<decltype(data+data), DomainType>&

    but this does not work if T is a double[3] or std::array<double, 3>
    as there is no + operator defined for these types
   */
  template <class T>
  typename std::enable_if<(!std::is_pointer<T>::value &&
                           !std::is_array<T>::value),
                          void>::type
  mesh_add_data(Entity_kind on_what, std::string const& var_name,
                const T value) {
    // Compiler needs some help deducing template parameters here
    jali_state_.add<T, Jali::Mesh, Jali::UniStateVector>(var_name,
                                                    jali_state_.mesh(),
                                                    (Jali::Entity_kind) on_what,
                                                    Jali::Entity_type::ALL,
                                                    value);
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and initialize its
   material data to a single value
   @param[in] var_name The name of the data field
   @param[in] value Initialize with this value

   The 2D array will be read and values copied according to which materials
   are contained in which cells. If a material+cell combination is not active
   providing a value for the array will have no effect. 

   This version of the overloaded operator is being DISABLED for
   pointer and array types (via the line 'typename
   std::enable_if....type') because template deduction rules are
   making the compiler invoke this version, when we call it with a
   const double ** pointer
   
   See, stackoverflow.com Q&A
   
   http://stackoverflow.com/questions/13665574/template-argument-deduction-and-pointers-to-constants

    We could make it work for some cases using

    template <class T, class DomainType,
              template<class, class> class StateVecType>
    auto add(........,
             T const& data) -> StateVecType<decltype(data+data), DomainType>&

    but this does not work if T is a double[3] or std::array<double, 3>
    as there is no + operator defined for these types
   */

  template <class T>
  typename std::enable_if<(!std::is_pointer<T>::value &&
                           !std::is_array<T>::value),
                          void>::type
  mat_add_celldata(std::string const& var_name, T value) {
    jali_state_.add<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                        jali_state_.mesh(),
                                                        Jali::Entity_kind::CELL,
                                                        Jali::Entity_type::ALL,
                                                        value);
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and initialize its
   material data according to a 2D array
   @param[in] var_name The name of the data field
   @param[in] layout  Whether 2D array is laid out with first index being the cell (CELL_CENRIC) or material (MATERIAL CENTRIC)
   @param[in] value Initialize with this value

   The 2D array will be read and values copied according to which materials
   are contained in which cells. If a material+cell combination is not active
   providing a value for the array will have no effect. 
   */
  template <class T>
  void mat_add_celldata(std::string const& var_name,
                        T const * const *values = nullptr,
                        Data_layout layout = Data_layout::MATERIAL_CENTRIC) {
    jali_state_.add(var_name, jali_state_.mesh(), Jali::Entity_kind::CELL,
                    Jali::Entity_type::ALL, (Jali::Data_layout) layout, values);
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and add data to one of
   its materials
   @param[in] var_name The name of the data field
   @param[in] matid  Index of material in the problem
   @param[in] layout Data layout - 
   @param[in] values Initialize with this array of values

   Subsequent calls to this function with the same name will find the added
   field and just add the data.
   */
  template <class T>
  void mat_add_celldata(std::string const& var_name, int matid,
                        T const * values) {
    auto it = jali_state_.find<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                        jali_state_.mesh(),
                                                        Jali::Entity_kind::CELL,
                                                        Jali::Entity_type::ALL);

    if (it == jali_state_.end()) {
      Jali::MultiStateVector<T, Jali::Mesh>& mmvec =
          jali_state_.add<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                  jali_state_.mesh(),
                                                  Jali::Entity_kind::CELL,
                                                  Jali::Entity_type::ALL);
      std::vector<T>& matdata = mmvec.get_matdata(matid);
      int nmatcells = matdata.size();
      matdata.assign(values, values+nmatcells);
    } else {
      Jali::MultiStateVector<T, Jali::Mesh>& mmvec =
          *(std::dynamic_pointer_cast<Jali::MultiStateVector<T, Jali::Mesh>>(*it));
      std::vector<T>& matdata = mmvec.get_matdata(matid);
      int nmatcells = matdata.size();
      matdata.assign(values, values+nmatcells);
    }
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and initialize one of
   its material data to a uniform value
   @param[in] var_name The name of the data field
   @param[in] matid Index of material in the problem
   @param[in] value Initialize with this value

   Subsequent calls to this function with the same name will find the added
   field and just add the data.

   Need to disable this call for deductions where T is a pointer;
   otherwise when we call mat_add_celldata with a pointer to a
   non-const type, this function is called by interpreting T as 'type
   *' (e.g. this overload is called with T = double *, instead of
   calling the pointer version with T = double)

   See stackoverflow.com Q&A

   http://stackoverflow.com/questions/13665574/template-argument-deduction-and-pointers-to-constants

   */
  template <class T>
  typename
  std::enable_if<(!std::is_pointer<T>::value && !std::is_array<T>::value),
                 void>::type 
  mat_add_celldata(std::string const& var_name, int matid, T value) {
    auto it = jali_state_.find<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                        jali_state_.mesh(),
                                                        Jali::Entity_kind::CELL,
                                                        Jali::Entity_type::ALL);

    if (it == jali_state_.end()) {
      Jali::MultiStateVector<T, Jali::Mesh>& mmvec =
          jali_state_.add<T, Jali::Mesh, Jali::MultiStateVector>(var_name,
                                                  jali_state_.mesh(),
                                                  Jali::Entity_kind::CELL,
                                                  Jali::Entity_type::ALL);
      std::vector<T>& matdata = mmvec.get_matdata(matid);
      int nmatcells = matdata.size();
      matdata.assign(nmatcells, value);
    } else {
      Jali::MultiStateVector<T, Jali::Mesh>& mmvec =
          *(std::dynamic_pointer_cast<Jali::MultiStateVector<T, Jali::Mesh>>(*it));
      std::vector<T>& matdata = mmvec.get_matdata(matid);
      int nmatcells = matdata.size();
      matdata.assign(nmatcells, value);
    }
  }
    


  /*!
    @brief Add cells to material (or add material to cells)
    @param[in] matid  Material ID
    @param[in] newcells Vector of new cells in material
  */

  void mat_add_cells(int matid, std::vector<int> const& newcells) {
    jali_state_.add_cells_to_material(matid, newcells);
  }


  /*!
    @brief Remove cells from material (or remove material from cells)
    @param[in] matid  Material ID
    @param[in] matcells Vector of to be removed cells
  */
  
  void mat_rem_cells(int matid, std::vector<int> const& delcells) {
    jali_state_.rem_cells_from_material(matid, delcells);
  }


  /*!
    @brief Add a material to state
    @param[in] matname  Name of material
    @param[in] matcells Cells containing the material
   */
  void add_material(std::string const& matname,
                    std::vector<int> const& matcells) {
    jali_state_.add_material(matname, matcells);
  }

  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] var_name The string name of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   */
  Entity_kind get_entity(const std::string var_name) const {
    
    Jali::State::const_iterator it = jali_state_.find(var_name);
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::StateVectorBase> vector = *it;
      if (vector)
        return (Portage::Entity_kind) vector->entity_kind();
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    return Portage::UNKNOWN_KIND;
  }


  /*!
    @brief Get the data size for the given field
    @param[in] on_what  The entity type on which the data field is defined
    @param[in] var_name The string name of the data field
    @return The data size for the field with the given name on the given entity type
         
    For multi-material state, this will give the number of materials for now
   */
  int get_data_size(Entity_kind on_what, std::string const& var_name) const {

    Jali::State::const_iterator it = jali_state_.find(var_name);
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::StateVectorBase> vector = *it;
      if (vector && on_what == (Portage::Entity_kind) vector->entity_kind()) {
        if (vector->type() == Jali::StateVector_type::UNIVAL) {
          auto uvec = std::dynamic_pointer_cast<Jali::UniStateVectorBase<Jali::Mesh>>(vector);
          if (uvec && uvec->domain() == jali_state_.mesh())
            return uvec->size();
        } else if (vector->type() == Jali::StateVector_type::MULTIVAL) {
          auto mmvec = std::dynamic_pointer_cast<Jali::MultiStateVectorBase<Jali::Mesh>>(vector);
          if (mmvec && mmvec->domain() == jali_state_.mesh())
              return (vector->size());
        }
      }
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    return 0;
  }


  /*!
    @brief Get the data type of the given field
    @param[in] var_name The string name of the data field
    @return A reference to the type_info struct for the field's data type
   */
  const std::type_info& get_data_type(std::string const& var_name) const {

    Jali::State::const_iterator it =
        jali_state_.find(var_name, jali_state_.mesh());
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::StateVectorBase> vector = *it;
      if (vector)
        return vector->data_type();
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    return typeid(0);
  }

  /*!
    @brief Begin iterator on vector names
    @return Begin iterator on vector of strings
   */
  std::vector<std::string>::iterator names_begin() const {
    return jali_state_.names_begin();
  }

  /*!
    @brief End iterator on vector names
    @return End iterator on vector of strings
   */
  std::vector<std::string>::iterator names_end() const {
    return jali_state_.names_end();
  }

 private:

  Jali::State & jali_state_;

};  // Jali_State_Wrapper

}  // namespace Wonton

#endif  // JALI_STATE_WRAPPER_H_
