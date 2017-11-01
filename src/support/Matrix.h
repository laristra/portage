/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_MATRIX_H
#define PORTAGE_MATRIX_H

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
  @file Matrix.h
  @brief Matrix class for Portage 
*/ 

namespace Portage {

  static std::string ignore_msg="ignore";
  static std::string &ignore_msg_ref(ignore_msg);

class Matrix {
 public:
  Matrix() : Rows_(0), Columns_(0) {}
  
  Matrix(int const Rows, int const Columns) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);    // uninitialized
  }
  
  // Initialize to some constant value
  explicit Matrix(int const Rows, int const Columns,
                  double initval) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = initval;
  }
  
  // Initialize from vector of vectors - assume that each row vector
  // is of the same size - if its not, then some values will be
  // uninitialized
  
  explicit Matrix(std::vector<std::vector<double>> const& invals) {
    Rows_ = invals.size();
    Columns_ = invals[0].size();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = invals[i][j];
  }

  Matrix(Matrix const& M) : Rows_(M.rows()), Columns_(M.columns()) {
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = M[i][j];
  }

  Matrix& operator=(Matrix const& M) {
    Rows_ = M.rows();
    Columns_ = M.columns();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = M[i][j];
    return *this;
  }

  // Destructor
  ~Matrix() {}

  //! Number of rows
  int rows() const {return Rows_;}

  //! Number of columns
  int columns() const {return Columns_;}

  //! Return a row of values
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of the matrix
  
  double * operator[](int const RowIndex) {
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Return a row of values that cannot be modified
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of a const matrix
  
  double const * operator[](int const RowIndex) const {
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Get a transpose
  Matrix transpose() const {
    Matrix AT(Columns_, Rows_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        AT[j][i] = A_[i*Columns_+j];
    return AT;
  }

  //! Get Inverse - only if its a square matrix

  Matrix inverse() const {
    if (Rows_ != Columns_) {
      std::cerr << "Matrix is rectangular" << std::endl;
      throw std::exception();
    }
    
    Matrix Ainv(Rows_, Rows_, 0.0);

    // Create a temporary matrix with twice as many columns
    Matrix D(Rows_, 2*Rows_);
    
    // Initialize the reduction matrix
    int Rows2 = 2*Rows_;
    for (int i = 0; i < Rows_; i++) {
      for (int j = 0; j < Rows_; j++) {
        D[i][j] = A_[i*Columns_+j];
        D[i][Rows_+j] = 0.0;
      }
      D[i][Rows_+i] = 1.0;
    }
    
    // Do the reduction
    for (int i = 0; i < Rows_; i++) {
      double alpha = D[i][i];
      if (alpha == 0.0) {
        std::cerr << "invert: Singular Matrix" << std::endl;
        return Ainv;
      }
      
      for (int j = 0; j < Rows2; j++)
        D[i][j] = D[i][j]/alpha;
      
      for (int k = 0; k < Rows_; k++) {
        if ((k - i) == 0)
          continue;
        
        double beta = D[k][i];
        for (int j = 0; j < Rows2; j++)
          D[k][j] = D[k][j] - beta*D[i][j];
      }
    }
    
    // Copy result into output matrix
    for (int i = 0; i < Rows_; i++)
      for (int j = 0; j < Rows_; j++)
        Ainv[i][j] = D[i][j + Rows_];

    return Ainv;
  }
  
  /*!
    @brief  Matrix-Vector multiply with std::vector
    @param[in] X The vector to post-multiply with
    
  */
  
  std::vector<double> operator*(std::vector<double> const& X) {
    assert(Columns_ == X.size());

    std::vector<double> AX(Rows_);
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Vector multiply with Portage::Vector
    @param[in] X The vector to post-multiply with
    
  */
  
  template<long D>
  Vector<D> operator*(Vector<D> const& X) {
    assert(Rows_ == D && Columns_ == D);

    Vector<D> AX;
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Matrix multiply
    @param[in] B   matrix to post-multiply with
  */
  
  Matrix operator*(Matrix const& B) {
    assert(Columns_ == B.rows());
    
    Matrix AB(Rows_, B.columns(), 0.0);
    
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < B.columns(); ++j)
        for (int k = 0; k < Columns_; ++k)
          AB[i][j] += A_[i*Columns_+k]*B[k][j];
    
    return AB;
  }

  /*!
    @brief  Matrix-scalar multiply
    @param[in] b   scalar to post-multiply with
  */

  Matrix operator*(double b) {

    Matrix Ab(Rows_, Columns_);

    for (int i = 0; i < Rows_; ++i)
        for (int k = 0; k < Columns_; ++k)
          Ab[i][k] = A_[i*Columns_+k]*b;

    return Ab;
  }

  /*!
    @brief  solve a linear system A X = B with this matrix (A)
    @param[in] B  right-hand sides (multiple)
    @param[in] method what method to use for solution
    @param[in,out] error message, if any 
    @return the solution X

    method=="inverse" ==> use the  inverse operator
    method=="dposv" ==> use lapack dposvx for symmetric positive definite A.
    method=="dsysv" ==> use lapack dsysvx for symmetric  A.
    method=="dgesv" ==> use lapack dgesvx for general A.

    If \code error\endcode is not present or has value "ignore", no message will be returned.
      The value of  \code error\endcode will be "ignore" on return.
    If \code error\endcode is present and has a value other than "ignore", the value "none" 
      will be returned if no error was generated, or contain the appropriate error message.
  */
  Matrix solve(Matrix const& B,
               std::string method="inverse",
	       std::string &error=ignore_msg_ref)
  {
    assert(Rows_ == Columns_);
    assert(B.rows() == Columns_);
    method_ = method;
    std::stringstream infoword;

    Matrix X(B.rows(), B.columns(), 0.);

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

 private:
  int Rows_, Columns_;
  std::vector<double> A_;
  std::string method_;

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
};

template<long D>
Vector<D> matsolve(Matrix const& matrix, Vector<D> const& rhs) {
  auto inverse = matrix.inverse();
  Vector<D> result = inverse*rhs;
  return result;
}

}  // namespace Portage

#endif // PORTAGE_MATRIX_H
