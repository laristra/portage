/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#ifndef SRC_SIMPLE_MESH_SIMPLE_STATE_H_
#define SRC_SIMPLE_MESH_SIMPLE_STATE_H_

#include <vector>
#include <memory>
#include <string>
#include <utility>

#include "portage/simple_mesh/simple_mesh.h"

#include "portage/support/portage.h"


namespace Portage {

class Simple_StateVector {
 public:
  Simple_StateVector(std::string const name,
                     Entity_kind const on_what,
                     std::shared_ptr<std::vector<double>> const data) :
  name_(name), on_what_(on_what), data_(data) { }

  ~Simple_StateVector() { }

  void* get_raw_data() { return reinterpret_cast<void*>(&((*data_)[0])); }

  std::shared_ptr<std::vector<double>> get_data() { return data_; }

  typedef std::vector<double>::iterator iterator;
  typedef std::vector<double>::const_iterator const_iterator;

  iterator begin() { return data_->begin(); }
  iterator end() { return data_->end(); }
  const_iterator cbegin() { return data_->cbegin(); }
  const_iterator cend() { return data_->cend(); }

  typedef double& reference;
  typedef double const& const_reference;
  reference operator[](int i) { return (*data_)[i]; }
  const_reference operator[](int i) const { return (*data_)[i]; }

  int size() const { return data_->size(); }
  std::string name() const { return name_; }
  Entity_kind on_what() const { return on_what_; }

  void clear() { data_->clear(); }

 private:
  std::string name_;
  Entity_kind on_what_;
  std::shared_ptr<std::vector<double>> data_;
};  // Simple_StateVector


class Simple_State {
 public:
  explicit Simple_State(const std::shared_ptr<Simple_Mesh> mesh) : mesh_(mesh)
  { }

  ~Simple_State() { }

  // Disabled functions
  Simple_State(const Simple_State & state) = delete;
  Simple_State & operator=(const Simple_State & state) = delete;

  std::shared_ptr<Simple_Mesh> mesh() { return mesh_; }

  typedef std::vector<std::shared_ptr<Simple_StateVector>>::iterator it;
  typedef std::vector<std::shared_ptr<Simple_StateVector>>::const_iterator c_it;
  typedef std::vector<std::string>::iterator string_it;

  it begin() { return state_vectors_.begin(); }
  it end() { return state_vectors_.end(); }
  c_it cbegin() { return state_vectors_.cbegin(); }
  c_it cend() { return state_vectors_.cend(); }
  string_it names_begin() { return names_.begin(); }
  string_it names_end() { return names_.end(); }

  int size() const { return state_vectors_.size(); }

  it find(std::string const name,
          Entity_kind const on_what = Entity_kind::ANY_KIND) {
    it iit = begin();
    while (iit != end()) {
      std::shared_ptr<Simple_StateVector> sv_ptr = *iit;
      if ((sv_ptr->name() == name) &&
          ((on_what == Entity_kind::ANY_KIND) ||
           (sv_ptr->on_what() == on_what)))
        break;
      else
        ++iit;
    }
    return iit;
  }

  int numVars() const { return state_vectors_.size(); }

  /* c_it find(std::string const name, */
  /*           Entity_kind const on_what = Entity_kind::ANY_KIND) const { */
  /*   c_it cit = cbegin(); */
  /*   while (cit != cend()) { */
  /*     std::shared_ptr<Simple_StateVector> const sv_ptr = *cit; */
  /*     if ((sv_ptr->name() == name) && */
  /*         ((on_what == Entity_kind::ANY_KIND) || */
  /*          (sv_ptr->on_what() == on_what))) */
  /*       break; */
  /*     else */
  /*       ++cit; */
  /*   } */
  /*   return cit; */
  /* } */

  Simple_StateVector& add(std::string const name,
                          Entity_kind const on_what,
                          double const * const data = nullptr) {
    auto iit = find(name, on_what);
    if (iit == end()) {
      // this is a new data vector, add it to our known vectors
      int num_entities = mesh_->num_entities(on_what,
                                             Entity_type::PARALLEL_OWNED);

      std::shared_ptr<std::vector<double>> dv;
      if (data == nullptr)
        dv = std::make_shared<std::vector<double>>(num_entities);
      else
        dv = std::make_shared<std::vector<double>>(data, data+num_entities);
      state_vectors_.emplace_back(name, on_what, dv);

      names_.emplace_back(name);
      return (*state_vectors_.back());
    } else {
      // already have this vector!
      std::cerr << "Attempted to add duplicate state vector.  Ignoring."
                << std::endl;
      return **iit;
    }
  }

  bool get(std::string const name,
           Entity_kind const on_what,
           std::shared_ptr<Simple_StateVector> *sv_ptr) {
    auto iit = find(name, on_what);
    if (iit != end()) {
      *sv_ptr = *iit;
      return true;
    } else {
      return false;
    }
  }


 private:
  std::shared_ptr<Simple_Mesh> mesh_;
  std::vector<std::shared_ptr<Simple_StateVector>> state_vectors_;
  std::vector<std::string> names_;
};  // Simple_State



/*   typedef std::pair<std::string, Entity_kind> field_ID; */
/*   typedef std::vector<field_ID> field_ID_list; */
/*   typedef field_ID_list::const_iterator field_ID_iterator; */




/*   std::vector<double>& add(std::string const name, */
/*                            Entity_kind const on_what, */
/*                            double const * const data = nullptr) { */
/*     auto it = find(name, on_what); */
/*     if (it == known_fields_.end()) { */
/*       // This wasn't found in our list; add it */
/*       auto numEntities = mesh_->num_entities(on_what, */
/*                                              Entity_type::PARALLEL_OWNED); */
/*       if (data == nullptr) */
/*         data_.emplace_back(numEntities); */
/*       else */
/*         data_.emplace_back(data, data + numEntities); */
/*       known_fields_.emplace_back(name, on_what); */
/*       return data_[data_.size()-1]; */
/*     } else { */
/*       std::cerr << "Attempted to add duplicate state vector; ignoring..." */
/*                 << std::endl; */
/*       return data_[it - known_fields_.begin()]; */
/*     } */
/*   } */

/*   int numVars() const { return known_fields_.size(); } */


/*  private: */
/*   field_ID_list known_fields_; */
/*   std::vector<std::vector<double>> data_; */
/*   const std::shared_ptr<Simple_Mesh> mesh_; */

/*   field_ID_iterator find(std::string const name, */
/*                          Entity_kind const on_what) const { */
/*     auto it = known_fields_.begin(); */
/*     while (it != known_fields_.end()) { */
/*       if ((it->first == name) && (it->second == on_what)) */
/*         break; */
/*       else */
/*         ++it; */
/*     } */
/*     return it; */
/*   } */
/* }; */

}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_STATE_H_
