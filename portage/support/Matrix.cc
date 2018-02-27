/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "Matrix.h"

#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <type_traits>
#include <string>
#include <sstream>

#ifdef HAVE_LAPACKE
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include "lapacke.h"
#endif  // LAPACKE

#include "Vector.h"   // Portage::Vector

/*!
  @file Matrix.cc
  @brief Matrix class for Portage 
*/ 

namespace Portage {


/*!
  @brief  solve a linear system A X = B with this matrix (A)
  @param[in] B  right-hand sides (multiple)
  @param[in] method what method to use for solution
  @param[in,out] error message, if any 
  @return the solution X
  
  method=="inverse" ==> use the  inverse operator
  method=="lapack-posv" ==> use lapack dposvx for symmetric positive definite A.
  method=="lapack-sysv" ==> use lapack dsysvx for symmetric A.
  method=="lapack-gesv" ==> use lapack dgesvx for general A.
  method=="lapack-sytr" ==> use lapack dsytrf & dsytrs for symmetric A.
  
  If \code error\endcode is not present or has value "ignore", no
  message will be returned.
  The value of  \code error\endcode will be "ignore" on return.
  If \code error\endcode is present and has a value other than
  "ignore", the value "none" will be returned if no error was
  generated, or contain the appropriate error message.
  */

Matrix Matrix::solve(Matrix const& B, std::string method, std::string &error) {
  assert(Rows_ == Columns_);
  assert(B.rows() == Columns_);
  method_ = method;
  std::stringstream infoword;
  
  Matrix X(B.rows(), B.columns(), 0.);

#ifdef HAVE_LAPACKE
  // stuff for lapack which is useful for solves where matrix is already
  // factored
  std::vector<double> AFactored_;
  std::vector<double> AEquilibrated_;
  std::vector<double> Scale_;
  std::vector<double> Scale2_;
  std::vector<double> Error_;
  std::vector<double> BackError_;
  std::vector<lapack_int> Pivot_;
  std::vector<double> Rpivot_;
#endif


  // built-in inverse
  if (method == "inverse") {
    auto inverse = this->inverse();
    X = inverse*B;
    
  }
#ifdef HAVE_LAPACKE    
  else if (method == "lapack-posv") {  // LAPACK positive-definite matrix
    
    // check symmetric
    bool symm = true;
    for (size_t i=0; i<Rows_; i++) for (size_t j=i; j<Columns_; j++) {
        symm = symm and ((*this)[i][j] - (*this)[j][i]) < 1.e-13;
      }
    if (not symm) std::cerr << "solve(dposv): matrix is not symmetric" << std::endl;
    assert(symm);
    
    // setup
    AEquilibrated_ = std::vector<double>(A_);
    AFactored_     = std::vector<double>(Rows_*Columns_,0.);
    Scale_         = std::vector<double>(Rows_,1.);
    Error_         = std::vector<double>(B.columns(),0.);
    BackError_     = std::vector<double>(B.columns(),0.);
    Matrix BT = B.transpose();
    Matrix XT = X.transpose();
    
    // The data for this matrix class is in row-major form.
    // LAPACKE creates temporaries and transposes the data, which does not
    // work on some systems. We transpose the input data ourselves and go
    // into LAPACKE in column-major form for direct access to the fortran
    // routines.
    
    int        matrix_layout =  LAPACK_COL_MAJOR;
    char       fact = 'E';
    char       uplo = 'U';
    lapack_int n = Rows_;
    lapack_int nrhs = B.columns();
    double    *a = &(AEquilibrated_[0]);
    lapack_int lda = Rows_;
    double    *af = &(AFactored_[0]);
    lapack_int ldaf = Rows_;
    char       equed[2] = {' ','\n'};
    double    *s = &(Scale_[0]);
    double    *b = const_cast<double *>(BT[0]);
    lapack_int ldb = Rows_;
    double    *x = XT[0];
    lapack_int ldx = Rows_;
    double     rcond;
    double    *ferr = &(Error_[0]);
    double    *berr = &(BackError_[0]);
    
    // solve it
    lapack_int info =  LAPACKE_dposvx ( matrix_layout,
        fact,     uplo,       n,          nrhs,        a,        lda,
        af,       ldaf,       equed,      s,           b,        ldb,
        x,        ldx,        &rcond,     ferr,        berr
      );

    X = XT.transpose();
    
    // checks
    if (error != "ignore") {
      error = "none";
      if (info <0) {
        infoword << -info;
        error = 
	    std::string("solve(posv): illegal value in ")+infoword.str()+"-th position";
      }
      if (info>0 and info<=n) {
        infoword << info;
        error = 
	    std::string("solve(posv): leading minor ")+infoword.str()+" is not positive definite";
      }
      if (info == n+1) {
        error = 
	    std::string("solve(posv): reciprocal condition number is less than machine precision");
      }
    }

  } else if (method == "lapack-sysv") {  // LAPACK symmetric matrix
    // check symmetric
    bool symm = true;
    for (size_t i=0; i<Rows_; i++) for (size_t j=i; j<Columns_; j++) {
        symm = symm and ((*this)[i][j] - (*this)[j][i]) < 1.e-13;
      }
    if (not symm) std::cerr << "solve(dposv): matrix is not symmetric" << std::endl;
    assert(symm);
    
    // setup
    AEquilibrated_ = std::vector<double>(A_);
    AFactored_     = std::vector<double>(Rows_*Columns_,0.);
    Error_         = std::vector<double>(B.columns(),0.);
    BackError_     = std::vector<double>(B.columns(),0.);
    Pivot_         = std::vector<lapack_int>(Rows_,0);
    Matrix BT = B.transpose();
    Matrix XT = X.transpose();
    
    // The data for this matrix class is in row-major form.
    // LAPACKE creates temporaries and transposes the data, which does not
    // work on some systems. We transpose the input data ourselves and go
    // into LAPACKE in column-major form for direct access to the fortran
    // routines.
    
    int        matrix_layout =  LAPACK_COL_MAJOR;
    char       fact = 'N';
    char       uplo = 'U';
    lapack_int n = Rows_;
    lapack_int nrhs = B.columns();
    double    *a = &(AEquilibrated_[0]);
    lapack_int lda = Rows_;
    double    *af = &(AFactored_[0]);
    lapack_int ldaf = Rows_;
    lapack_int *ipiv = &(Pivot_[0]);
    double    *b = const_cast<double *>(BT[0]);
    lapack_int ldb = Rows_;
    double    *x = XT[0];
    lapack_int ldx = Rows_;
    double     rcond;
    double    *ferr = &(Error_[0]);
    double    *berr = &(BackError_[0]);
    
    // solve it
    lapack_int info =  LAPACKE_dsysvx ( matrix_layout,
        fact,     uplo,       n,          nrhs,        a,        lda,
        af,       ldaf,       ipiv,       b,           ldb,
        x,        ldx,        &rcond,     ferr,        berr
      );

    X = XT.transpose();
    
    // checks
    if (error != "ignore") {
      error = "none";
      if (info <0) {
        infoword<<-info;
        error =  std::string("solve(sysv): illegal value in ")+infoword.str()+"-th position";
      }
      if (info>0 and info<=n) {
        infoword<<info;
        error = 
	    std::string("solve(sysv): diagonal entry ")+infoword.str()+" is zero";
      }
      if (info == n+1) {
        error = 
	    std::string("solve(sysv): reciprocal condition number is less than machine precision");
      }
    }
    
  } else if (method == "lapack-gesv") {  // LAPACK general matrix

    // setup
    AEquilibrated_ = std::vector<double>(A_);
    AFactored_     = std::vector<double>(Rows_*Columns_,0.);
    Scale_         = std::vector<double>(Rows_,1.);
    Scale2_        = std::vector<double>(Rows_,1.);
    Error_         = std::vector<double>(B.columns(),0.);
    BackError_     = std::vector<double>(B.columns(),0.);
    Pivot_         = std::vector<lapack_int>(Rows_,0.);
    Rpivot_        = std::vector<double>(Rows_,0.);
    Matrix BT = B.transpose();
    Matrix XT = X.transpose();
    
    // The data for this matrix class is in row-major form.
    // LAPACKE creates temporaries and transposes the data, which does not
    // work on some systems. We transpose the input data ourselves and go
    // into LAPACKE in column-major form for direct access to the fortran
    // routines.
    
    int        matrix_layout =  LAPACK_COL_MAJOR;
    char       fact = 'E';
    char       trans = 'N';
    lapack_int n = Rows_;
    lapack_int nrhs = B.columns();
    double    *a = &(AEquilibrated_[0]);
    lapack_int lda = Rows_;
    double    *af = &(AFactored_[0]);
    lapack_int ldaf = Rows_;
    lapack_int *ipiv = &(Pivot_[0]);
    char       equed[2] = {' ','\n'};
    double    *r = &(Scale_[0]);
    double    *c = &(Scale2_[0]);
    double    *b = const_cast<double *>(BT[0]);
    lapack_int ldb = Rows_;
    double    *x = XT[0];
    lapack_int ldx = Rows_;
    double     rcond;
    double    *ferr = &(Error_[0]);
    double    *berr = &(BackError_[0]);
    double    *rpivot = &(Rpivot_[0]);
    
    // solve it
    lapack_int info =  LAPACKE_dgesvx ( matrix_layout,
        fact,     trans,      n,          nrhs,        a,        lda,
        af,       ldaf,       ipiv,       equed,       r,        c,
        b,           ldb,     x,          ldx,        &rcond,
        ferr,        berr,    rpivot
      );

    X = XT.transpose();
    
    // checks
    if (error != "ignore") {
      error = "none";
      if (info <0) {
        infoword<<-info;
        error =  std::string("solve(gesv): illegal value in ")+infoword.str()+"-th position";
      }
      if (info>0 and info<=n) {
        infoword<<info;
        error = 
	    std::string("solve(gesv): upper triangle entry ")+infoword.str()+" is zero";
      }
      if (info == n+1) {
        error = 
	    std::string("solve(gesv): reciprocal condition number is less than machine precision");
      } 
    }
    
  } else if (method == "lapack-sytr") {  // LAPACK symmetric matrix
    
    // check symmetric
    bool symm = true;
    for (size_t i=0; i<Rows_; i++) for (size_t j=i; j<Columns_; j++) {
        symm = symm and ((*this)[i][j] - (*this)[j][i]) < 1.e-13;
      }
    if (not symm) std::cerr << "solve(sytr): matrix is not symmetric" << std::endl;
    assert(symm);
    
    // setup
    AEquilibrated_ = std::vector<double>(A_);
    Pivot_         = std::vector<lapack_int>(Rows_,0);
    Matrix XT(B.transpose());
    
    // The data for this matrix class is in row-major form.
    // LAPACKE creates temporaries and transposes the data, which does not
    // work on some systems. We transpose the input data ourselves and go
    // into LAPACKE in column-major form for direct access to the fortran
    // routines.
    
    int        matrix_layout =  LAPACK_COL_MAJOR;
    char       uplo = 'U';
    lapack_int n = Rows_;
    lapack_int nrhs = B.columns();
    double    *a = &(AEquilibrated_[0]);
    lapack_int lda = Rows_;
    lapack_int *ipiv = &(Pivot_[0]);
    double    *b = XT[0];
    lapack_int ldb = Rows_;
    lapack_int info;
    
    // factorize a
    info = LAPACKE_dsytrf(matrix_layout,uplo,n,a,lda,ipiv);
    
    // checks
    bool skip=false;
    if (error != "ignore") {
      error = "none";
      if (info < 0) {
        infoword<<-info;
        error = std::string("solve(sytr): illegal value in ")+infoword.str()+"-th position";
        skip = true;
      } else if (info > 0) {
        infoword<<info;
        error = std::string("solve(sytr): diagonal entry ")+infoword.str()+" is zero";
        skip = true;
      }
    }
    
    // solve it
    if (not skip) {
      info = LAPACKE_dsytrs(matrix_layout,uplo,n,nrhs,a,lda,ipiv,b,ldb);
      
      // checks
      if (error != "ignore") {
        error = "none";
        if (info < 0) {
          infoword<<-info;
          error = std::string("solve(sytr): illegal value in ")+infoword.str()+"-th position";
        }
      }
    }
    
    X = XT.transpose();
    
  } else {
    throw std::runtime_error(std::string("LAPACK solve requested but solver ")+method+" unrecognized");
  }
#else
  else {
    throw std::runtime_error("Invalid solver type or LAPACK solve requested but code not built with LAPACK support");
  }
#endif  // HAVE_LAPACKE

  return X;
}

}  // namespace Portage


