/*---------------------------------------------------------------------------~*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
 *---------------------------------------------------------------------------~*/
#ifndef LRE_PILE_INCLUDED
#define LRE_PILE_INCLUDED

#include <vector>
#include <valarray>
#include <cassert>

namespace Portage {
namespace swarm {
namespace Pairs {
  using std::vector;
  using std::valarray;

  //\///////////////////////////////////////////////////////////////////////////
  // pile typedef - definition 
  //\///////////////////////////////////////////////////////////////////////////

  /// basic pile of numbers, typedef for valarray<double>
  /** The pile collection of classes implements c-style tensors, e.g. 
      \code a[i][j][k] \endcode
      using STL valarrays, with a handy API which makes it easy to specifiy 
      sizes and initial values. 
  */
  typedef valarray<double> pile;

  //\///////////////////////////////////////////////////////////////////////////
  // vector pile class - definition 
  //\///////////////////////////////////////////////////////////////////////////
 
  /// vector of piles
  /** Class vpile provides a convenient wrapper for the STL class 
    \code valarray<valarray<double>> \endcode 
    which allows simple constructors and f90-style array operations. The code 
    \code 
    vpile a(3,10);
    a[0] = 2.;
    a[1] = 4.;
    a[2] = a[0]*a[1];
    for(size_t i=0;i<10;i++) cout<<a[2][i]; 
    \endcode
    produces the output:

    8888888888

    Each element of a vpile is a pile.
  */

  class vpile: public valarray< pile > {
  private:
    size_t nval;  ///< size of piles
    size_t nvec;  ///< number of elements 
 
  public:
    // constructors
    vpile();
    vpile(size_t, size_t);  
    vpile(size_t, const pile&);
    vpile(const vpile&);
    vpile(size_t, size_t, double**);

    // destructor
    ~vpile();

    // operators
    vpile& operator=(const vpile&);

    // size information
    vector<size_t> size() const;

    // change sizes
    void resize(const size_t, const size_t);
  };

  //\///////////////////////////////////////////////////////////////////////////
  // vector pile class - implementation 
  //\///////////////////////////////////////////////////////////////////////////

  /** default constructor */
  inline vpile::vpile(): valarray<pile>(){nval = nvec = 0;}
      
  /** \code vpile a(m,n) \endcode contains \a m piles of size \a n each 
   */
  inline vpile::vpile(size_t m, size_t n) : 
    valarray< pile >(pile(n), m) {nval=n; nvec=m;} 
	
  /** \code vpile a(m,v) \endcode contains \a m piles that are duplicates 
     of \a v 
  */
  inline vpile::vpile(size_t m, const pile& v) : 
    valarray< pile >(v, m) {nval=v.size(); nvec=m;}

  /** copy constructor 
   */
  inline vpile::vpile(const vpile& vp) : 
    valarray< pile >(vp) {nval=vp.nval; nvec=vp.nvec;}
	
  /** \code vpile a(m,n,data) \endcode converts a two-dimensional
      C-style array to a vpile with \code a[i][j]=data[i][j] \endcode
  */
  inline vpile::vpile(size_t m, size_t n, double** data) : 
		nvec(m), nval(n), valarray< pile >(pile(n), m)
    {
      for (size_t i=0; i<m; i++) for (size_t j=0; j<n; j++) {
	(*this)[i][j] = data[i][j];
      }
    }
	
  /** default destructor */
  inline vpile::~vpile(){}

  /** unsurprising operator= */
  inline vpile& vpile::operator=(const vpile& vp) {
    nval=vp.nval; nvec=vp.nvec; 
    vpile& vpnc=const_cast<vpile&>(vp);
    *(dynamic_cast< valarray<pile>* >(this))=
      dynamic_cast< valarray<pile>& >(vpnc); 
    return *this;
  }

  /** returns \code vector<size_t> s(2)\endcode where \code s(0) \endcode is 
      \c nvec and \code s(1) \endcode is \c nval */
  inline vector<size_t> vpile::size() const{
    vector<size_t> tmp(2); 
    tmp[0] = static_cast<size_t>(nvec);
    tmp[1] = static_cast<size_t>(nval);
    return tmp;
  }

  /** resize the arrays */
  inline void vpile::resize(const size_t s0, const size_t s1){
      nvec = s0;
      nval = s1;
      this->valarray<pile>::resize(s0);
      for (size_t i=0; i<s0; i++) this->operator[](i).resize(s1);
  }

  //\///////////////////////////////////////////////////////////////////////////
  // matrix pile class - definition 
  //\///////////////////////////////////////////////////////////////////////////

  /// matrix of piles
  /** Class vvpile provides a convenient wrapper for the STL class 
      \code valarray<valarray<valarray<double>>> \endcode 
      which allows simple constructors and f90-style array operations. The code 
      \code 
      vvpile a(2,2,10);
      a[0][0] = 2.;
      a[1][0] = 1.;
      a[0][1] = 3.;
      a[1][1] = 4.;
      pile p(10);
      p = a[0][0]*a[1][1]-a[0][1]*a[1][0];
      for(size_t i=0;i<10;i++) cout<<p[i]; 
      \endcode
      produces the output:

      5555555555

      Each element of an vvpile is a vpile. A element of an element of an vvpile 
      is a pile.
  */
  class vvpile: public valarray< vpile > {
  private:
    size_t nval;              ///< size of piles
    size_t nvec0;             ///< range of first index
    size_t nvec1;             ///< range of second index

  public:
    // constructors
    vvpile();
    vvpile(size_t, size_t, size_t);
    vvpile(size_t, size_t, const pile&);
    vvpile(const vvpile&);
    vvpile(size_t, size_t, size_t, double***);

    // destructor
    ~vvpile();

    // operators
    vvpile& operator=(const vvpile&);

    // size information
    vector<size_t> size() const;

    // change sizes
    void resize(const size_t, const size_t, const size_t);
  };

  //\///////////////////////////////////////////////////////////////////////////
  // matrix pile class - implementation 
  //\///////////////////////////////////////////////////////////////////////////

  /** default constructor */
  inline vvpile::vvpile(): valarray<vpile>() {nval = nvec0 = nvec1 = 0;}
      
  /** \code vvpile a(k,m,n) \endcode contains \a k by \a m piles of size \a n each 
   */
  inline vvpile::vvpile(size_t k, size_t m, size_t n) : 
    valarray<vpile>(vpile(m,n), k) 
    {nval=n; nvec1=m; nvec0=k;} 
		
  /** \code vvpile a(k,m,v) \endcode contains \a k by \a m piles that are duplicates 
     of \a v 
  */
  inline vvpile::vvpile(size_t k, size_t m, const pile& v) : 
    valarray<vpile>(vpile(m,v), k) 
    {nval=v.size(); nvec1=m; nvec0=k;}
	
  /** copy constructor */
  inline vvpile::vvpile(const vvpile& mp) : 
    valarray<vpile>(mp) 
    {nval=mp.nval; nvec1=mp.nvec1; nvec0=mp.nvec0;}

  /** \code vvpile a(k,m,n,data) \endcode converts a three-dimensional
      C-style array to an vvpile with \code a[i][j][k]=data[i][j][k] \endcode
  */
  inline vvpile::vvpile(size_t k, size_t m, size_t n, double*** data) : 
    valarray<vpile>(vpile(m,n), k)
    {
      nval=n; nvec1=m; nvec0=k;
      for (size_t h=0; h<k; h++) for (size_t i=0; i<m; i++) 
	for (size_t j=0; j<n; j++) {
	  (*this)[h][i][j] = data[h][i][j];
      }
    }
	
  /** default destructor */
  inline vvpile::~vvpile(){}

  /** unsurprising operator= */
  inline vvpile& vvpile::operator=(const vvpile& vvp) {
    nval=vvp.nval; nvec1=vvp.nvec1; nvec0=vvp.nvec0;
    vvpile& vvpnc=const_cast<vvpile&>(vvp);
    *(dynamic_cast< valarray<vpile>* >(this))=
      dynamic_cast< valarray<vpile>& >(vvpnc);
    return *this;
  }

  /** returns length 2 \code vector<size_t> v \endcode with size information:
      \code 
      vvpile a(2,3,4);
      vector<size_t> s(2);
      s=a.size();
      cout << s[0] << ',' << s[1];
      \endcode
      produces the result 
      2,3
  */
  inline vector<size_t> vvpile::size() const{
    vector<size_t> tmp(3); 
    tmp[0] = static_cast<size_t>(nvec0);
    tmp[1] = static_cast<size_t>(nvec1);
    tmp[2] = static_cast<size_t>(nval);
    return tmp;
  } 

  /** resize the arrays */
  inline void vvpile::resize(const size_t s0, const size_t s1, const size_t s2){
      nval = s2;
      nvec0 = s0;
      nvec1 = s1;
      this->valarray<vpile>::resize(s0);
      for (size_t i=0; i<s0; i++) this->operator[](i).resize(s1,s2);
  }

  //\///////////////////////////////////////////////////////////////////////////
  // 3D tensor pile class - definition 
  //\///////////////////////////////////////////////////////////////////////////

  /// tensor of piles
  /** Class vvvpile provides a convenient wrapper for the STL class 
      \code valarray<valarray<valarray<double>>> \endcode 
      which allows simple constructors and f90-style array operations. The code 
      \code 
      vvvpile a(2,2,2,10);
      a[0][0][0] = 2.;
      a[1][0] = 1.;
      a[0][0][0][1] = 3.;
      a[0][1][1] = 4.;
      pile p(10);
      p = a[0][0][0]*a[0][1][1]-a[0][0][1]*a[0][1][0];
      for(size_t i=0;i<10;i++) cout<<p[i]; 
      \endcode
      produces the output:

      5555555555

      Each element of an vvpile is a vpile. A element of a element of an vvpile is a pile.
  */
  class vvvpile: public valarray< vvpile > {
  private:
    size_t nval;              ///< size of piles
    size_t nvec0;             ///< range of first index
    size_t nvec1;             ///< range of second index
    size_t nvec2;             ///< range of third index

  public:
    // constructors
    vvvpile();
    vvvpile(size_t, size_t, size_t, size_t);
    vvvpile(size_t, size_t, size_t, const pile&);
    vvvpile(const vvvpile&);
    vvvpile(size_t, size_t, size_t, size_t, double****);

    // destructor
    ~vvvpile();

    // operators
    vvvpile& operator=(const vvvpile&);

    // other methods
    vector<size_t> size() const;
  };

  //\///////////////////////////////////////////////////////////////////////////
  // 3D tensor pile class - implementation 
  //\///////////////////////////////////////////////////////////////////////////

  /** default constructor */
  inline vvvpile::vvvpile(): valarray<vvpile>() {nval=nvec0=nvec1=nvec2=0;}
      
  /** \code vvvpile a(j,k,m,n) \endcode contains \a j by \a k by \a m piles of 
      size \a n each 
   */
  inline vvvpile::vvvpile(size_t j, size_t k, size_t m, size_t n) : 
    valarray<vvpile>(vvpile(k,m,n), j)
    {nval=n; nvec2=m; nvec1=k; nvec0=j;} 
		
  /** \code vvvpile a(j,k,m,v) \endcode contains \a j by \a k by \a m piles that are 
      duplicates 
     of \a v 
  */
  inline vvvpile::vvvpile(size_t j, size_t k, size_t m, const pile& v) : 
    valarray<vvpile>(vvpile(k,m,v), j)
    {nval=v.size(); nvec2=m; nvec1=k; nvec0=j;}
	
  /** copy constructor */
  inline vvvpile::vvvpile(const vvvpile& mp) : 
    valarray<vvpile>(mp) 
    {nval=mp.nval; nvec2=mp.nvec2; nvec1=mp.nvec1; nvec0=mp.nvec0;}

  /** \code vvvpile a(j,k,m,n,data) \endcode converts a four-dimensional
      C-style array to an vvvpile with \code a[i][j][k][m]=data[i][j][k][m] \endcode
  */
  inline vvvpile::vvvpile(size_t j, size_t k, size_t m, size_t n, double**** data) : 
    valarray<vvpile>(vvpile(k,m,n), j)
    {
      nval=n; nvec2=m; nvec1=k; nvec0=j;
      for (size_t g=0; g<j; g++) for (size_t h=0; h<k; h++) 
	for (size_t i=0; i<m; i++) for (size_t j=0; j<n; j++) {
	  (*this)[g][h][i][j] = data[g][h][i][j];
      }
    }
	
  /** default destructor */
  inline vvvpile::~vvvpile(){}

  /** unsurprising operator= */
  inline vvvpile& vvvpile::operator=(const vvvpile& vvvp) {
    nval=vvvp.nval; nvec1=vvvp.nvec1; nvec0=vvvp.nvec0;
    vvvpile& vvvpnc=const_cast<vvvpile&>(vvvp);
    *(dynamic_cast< valarray<vvpile>* >(this))=
      dynamic_cast< valarray<vvpile>& >(vvvpnc);
    return *this;
  }

  /** returns length 3 \code vector<size_t> v \endcode with size information:
      \code 
      vvvpile a(2,3,4);
      vector<size_t> s(2);
      s=a.size();
      cout << s[0] << " " << s[1];
      \endcode
      produces the result 
      2,3
  */
  inline vector<size_t> vvvpile::size() const{
    vector<size_t> tmp(3); 
    tmp[0] = static_cast<size_t>(nvec0);
    tmp[1] = static_cast<size_t>(nvec1);
    tmp[2] = static_cast<size_t>(nvec2);
    tmp[3] = static_cast<size_t>(nval);
    return tmp;
  }

  //\///////////////////////////////////////////////////////////////////////////
  // pile functions
  //\///////////////////////////////////////////////////////////////////////////

  /// dot product of a vpile and a vpile, result is a pile stored in \a r 
  inline void dot(pile& r, vpile &v1, vpile &v2) {
    assert(v1.size()[0]==v2.size()[0]);
    assert(v1.size()[1]==r.size() && v1.size()[1]==v2.size()[1]);
    r = 0.;
    for(size_t i=0; i<v1.size()[0]; i++) r+=v1[i]*v2[i];
  }

  /// dot product of a vpile and a vpile, function version, result is pile 
  inline pile dot(vpile &v1, vpile &v2) {
    pile r(v1.size()[0], v1.size()[1]);
    dot(r, v1, v2);
    return r;
  }

  /// dot product of an vvpile and a vpile, result is a vpile stored in \a r 
  inline void dot(vpile &r, vvpile &m, vpile &v) {
    assert(m.size()[0]==r.size()[0] && m.size()[1]==v.size()[0]);
    assert(m.size()[2]==r.size()[1] && m.size()[2]==v.size()[1]);
    for(size_t i=0; i<m.size()[0]; i++) {
      r[i] = 0.;
      for (size_t j=0; j<m.size()[1]; j++) r[i]+=m[i][j]*v[j];
    }
  }

  /// dot product of an vvpile and a vpile, function version, result is vpile
  inline vpile dot(vvpile &m, vpile &v) {
    vpile r(m.size()[0],m.size()[2]);
    dot(r, m, v);
    return r;
  }

  /// dot product of an vvpile and a vpile, result is a vpile stored in \a r 
  inline void dot(vpile &r, vpile &v, vvpile &m) {
    assert(m.size()[0]==v.size()[0] && m.size()[1]==r.size()[0]);
    assert(m.size()[2]==r.size()[1] && m.size()[2]==v.size()[1]);
    for(size_t i=0; i<m.size()[0]; i++) {
      r[i] = 0.;
      for (size_t j=0; j<m.size()[0]; j++) r[i]+=v[j]*m[j][i];
    }
  }

  /// dot product of an vvpile and a vpile, function version, result is vpile
  inline vpile dot(vpile &v, vvpile &m) {
    vpile r(m.size()[1],m.size()[2]);
    dot(r, v, m);
    return r;
  }

  /// dot product of a vvpile and a vvpile, result is a vvpile stored in \a r 
  inline void dot(vvpile& r, vvpile &m1, vvpile &m2) {
    assert(r.size()[2]==m1.size()[2] && r.size()[2]==m2.size()[2]);
    assert(r.size()[0]==m1.size()[0] && r.size()[1]==m2.size()[1]);
    assert(m1.size()[1]==m2.size()[0]);
    for(size_t i=0; i<m1.size()[0]; i++) for(size_t j=0; j<m2.size()[1]; j++) {
      r[i][j]=0.;
      for(size_t k=0; k<m1.size()[1]; k++) {
	r[i][j]+=m1[i][k]*m2[k][j];
      }
    }
  }

  /// dot product of a vvpile and a vvpile, result is a vvpile stored in \a r 
  inline vvpile dot(vvpile &m1, vvpile &m2) {
    vvpile r(m1.size()[0], m2.size()[1], m1.size()[2]);
    dot(r, m1, m2);
    return r;
  }
}
}
}

#endif


