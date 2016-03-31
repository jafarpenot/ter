#ifndef LINALG_FILE_SOLVE_MUMPS_CXX

#include "SolveMumps.hxx"

namespace linalg
{

  //! Mumps is called in double precision
  template<>
  void MatrixMumps<double>::CallMumps()
  {
    dmumps_c(&struct_mumps);
  }


  //! Mumps is called in complex double precision
  template<>
  void MatrixMumps<complex<double> >::CallMumps()
  {
    zmumps_c(&struct_mumps);
  }


  //! Default constructor
  template<class T>
  MatrixMumps<T>::MatrixMumps()
  {
    struct_mumps.comm_fortran = -987654;

    // parameters for mumps
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // 0 -> unsymmetric matrix

    struct_mumps.info[8] = 0;
    struct_mumps.info[9] = 0;

    // other parameters
    struct_mumps.n = 0;
    type_ordering = 7; // default : we let Mumps choose the ordering
    print_level = -1;
    info_facto = 0;
    out_of_core = false;
    coef_overestimate = 1.3;
    coef_increase_memory = 1.5;
    coef_max_overestimate = 50.0;
  }


  //! Function used to force factorisation when estimated space was too small
  template<class T>
  void MatrixMumps<T>::IterateFacto()
  {
    // if error -9 occurs, retrying with larger size
    int init_percentage = struct_mumps.icntl[13];
    int new_percentage = init_percentage;
    while (((struct_mumps.info[0] == -9) || (struct_mumps.info[0] == -8)
            || (struct_mumps.info[0] == -17) || (struct_mumps.info[0] == -20))
           && (new_percentage < coef_max_overestimate*100.0))
      {
        new_percentage = int(new_percentage*coef_increase_memory);
        struct_mumps.icntl[13] = new_percentage;
        struct_mumps.job = 2;
        CallMumps();
      }
    
    struct_mumps.icntl[13] = init_percentage;
    info_facto = struct_mumps.info[0];
  }
  
  
  //! Calls initialization routine provided by Mumps
  template<class T>
  void MatrixMumps<T>
  ::InitMatrix(bool sym, bool distributed)
  {
    // we clear previous factorization
    Clear();

    // symmetry is specified during the initialization stage
    struct_mumps.job = -1;
    if (sym)
      struct_mumps.sym = 2; // general symmetric matrix
    else
      struct_mumps.sym = 0; // unsymmetric matrix

    // mumps is called
    CallMumps();
    
    struct_mumps.icntl[13] = int(100.0*(coef_overestimate-1.0));
    struct_mumps.icntl[6] = type_ordering;
    
    // setting out of core parameters
    if (out_of_core)
      struct_mumps.icntl[21] = 1;
    else
      struct_mumps.icntl[21] = 0;

    struct_mumps.icntl[17] = 0;

    // the print level is set in mumps
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }


  //! Selects another ordering scheme.
  template<class T>
  void MatrixMumps<T>::SelectOrdering(int num_ordering)
  {
    type_ordering = num_ordering;
  }


  //! Destructor
  template<class T>
  MatrixMumps<T>::~MatrixMumps()
  {
    Clear();
  }


  //! Clears factorization
  template<class T>
  void MatrixMumps<T>::Clear()
  {
    if (struct_mumps.n > 0)
      {
	struct_mumps.job = -2;
	// Mumps variables are deleted.
        CallMumps();

	struct_mumps.n = 0;
      }
  }


  //! Informs Mumps that no message should be displayed
  template<class T>
  void MatrixMumps<T>::HideMessages()
  {
    print_level = -1;

    struct_mumps.icntl[0] = -1;
    struct_mumps.icntl[1] = -1;
    struct_mumps.icntl[2] = -1;
    struct_mumps.icntl[3] = 0;

  }


  //! Informs Mumps to display standard output
  template<class T>
  void MatrixMumps<T>::ShowMessages()
  {
    print_level = 0;

    struct_mumps.icntl[0] = 6;
    struct_mumps.icntl[1] = 0;
    struct_mumps.icntl[2] = 6;
    struct_mumps.icntl[3] = 2;

  }


  //! Enables writing on the disk (out of core).
  template<class T>
  void MatrixMumps<T>::EnableOutOfCore()
  {
    out_of_core = true;
  }


  //! Disables writing on the disk (incore).
  template<class T>
  void MatrixMumps<T>::DisableOutOfCore()
  {
    out_of_core = false;
  }

  
  //! Sets the coefficient used to overestimate the needed memory
  template<class T>
  void MatrixMumps<T>::SetCoefficientEstimationNeededMemory(double coef)
  {
    coef_overestimate = coef;
  }
  
  
  //! Sets the maximal allowed coefficient for the memory space multiplication
  template<class T>
  void MatrixMumps<T>::SetMaximumCoefficientEstimationNeededMemory(double coef)
  {
    coef_max_overestimate = coef;
  }
  
  
  //! Sets multiplication factor for each try to factorize the matrix
  template<class T>
  void MatrixMumps<T>::SetIncreaseCoefficientEstimationNeededMemory(double coef)
  {
    coef_increase_memory = coef;
  }
  

  //! Factorizes a given matrix
  /*!
    \param[in,out] mat matrix to factorize
    \param[in] sym symmetric matrix ?
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T>
  void MatrixMumps<T>::Factorize(SparseMatrix<T>& mat, bool sym, bool keep_matrix)
  {
    int n = mat.GetM();
    // conversion in coordinate format with fortran convention (1-index)
    Vector<int> num_row, num_col; Vector<T> values;
    if (sym)
      {
	int nnz = 0;
	for (int i = 0; i < mat.GetM(); i++)
	  for (int j = 0; j < mat.GetRowSize(i); j++)
	    if (i <= mat.Index(i, j))
	      nnz ++;
	
	num_row.Reallocate(nnz);
	num_col.Reallocate(nnz);
	values.Reallocate(nnz);
	nnz = 0;
	for (int i = 0; i < mat.GetM(); i++)
	  for (int j = 0; j < mat.GetRowSize(i); j++)
	    if (i <= mat.Index(i, j))
	      {
		num_row(nnz) = i+1;
		num_col(nnz) = mat.Index(i, j) + 1;
		values(nnz) = mat.Value(i, j);
		nnz++;
	      }
      }
    else
      {
	int nnz = 0;
	for (int i = 0; i < mat.GetM(); i++)
	  nnz += mat.GetRowSize(i);
	
	num_row.Reallocate(nnz);
	num_col.Reallocate(nnz);
	values.Reallocate(nnz);

	nnz = 0;
	for (int i = 0; i < mat.GetM(); i++)
	  for (int j = 0; j < mat.GetRowSize(i); j++)
	    {
	      num_row(nnz) = i+1;
	      num_col(nnz) = mat.Index(i, j) + 1;
	      values(nnz) = mat.Value(i, j);
	      nnz++;
	    }
	
      }

    if (!keep_matrix)
      mat.Clear();
    
    InitMatrix(sym);
    
    int nnz = values.GetM();
    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
    
    IterateFacto();
  }

  
  //! Returns memory used by the factorisation in bytes
  template<class T>
  size_t MatrixMumps<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
   
    if (struct_mumps.n <= 0)
      return taille;
    
    size_t nnz = struct_mumps.info[8];
    if (struct_mumps.info[8] < 0)
      nnz = abs(struct_mumps.info[8])*size_t(1024*1024);

    taille += sizeof(T)*nnz ;
    nnz = struct_mumps.info[9];
    if (struct_mumps.info[9] < 0)
      nnz = abs(struct_mumps.info[9])*int64_t(1024*1024);
    
    taille += sizeof(int)*nnz;
    return taille;
  }
  

  //! Returns information about factorization performed
  template<class T>
  int MatrixMumps<T>::GetInfoFactorization() const
  {
    return info_facto;
  }


  //! Solves a linear system using the computed factorization
  /*!
    \param[in,out] x right-hand-side on input, solution on output
    It is assumed that a call to Factorize has been done before
  */
  template<class T>
  void MatrixMumps<T>::Solve(Vector<T>& x)
  {
#ifdef LINALG_DEBUG
    if (x.GetM() != struct_mumps.n)
      throw WrongIndex("Mumps::Solve(TransA, c)",
		       string("The length of x is equal to ")
		       + to_string(x.GetM())
		       + " while the size of the matrix is equal to "
		       + to_string(struct_mumps.n) + ".");
#endif
    
    // no transpose
    struct_mumps.icntl[8] = 1;
    
    struct_mumps.nrhs = 1;
    struct_mumps.lrhs = struct_mumps.n;
    struct_mumps.rhs = reinterpret_cast<pointer>(x.GetData());
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }

  template class MatrixMumps<double>;
  template class MatrixMumps<complex<double> >;
  
}

#define LINALG_FILE_SOLVE_MUMPS_CXX
#endif
