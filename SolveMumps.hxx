#ifndef LINALG_FILE_SOLVE_MUMPS_HXX

// including Mumps headers
extern "C"
{
#include "dmumps_c.h"
#include "zmumps_c.h"
}

namespace linalg
{
  template<class T>
  class TypeMumps
  {
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<double>
  {
  public :
    typedef DMUMPS_STRUC_C data;
    typedef double* pointer;
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<complex<double> >
  {
  public :
    typedef ZMUMPS_STRUC_C data;
    typedef mumps_double_complex* pointer;
  };


  //! object used to solve linear system by calling mumps subroutines
  template<class T>
  class MatrixMumps
  {
  protected :
    //! object containing Mumps data structure
    typename TypeMumps<T>::data struct_mumps;    
    //! double* or complex<double>*
    typedef typename TypeMumps<T>::pointer pointer;
    int print_level;
    int type_ordering;
    int info_facto;
    bool out_of_core;
    double coef_overestimate;
    double coef_increase_memory;
    double coef_max_overestimate;

    // internal methods
    void CallMumps();
    void IterateFacto();
    
    void InitMatrix(bool sym, bool dist = false);

  public :
    MatrixMumps();
    ~MatrixMumps();

    void Clear();

    void HideMessages();
    void ShowMessages();

    void SelectOrdering(int num_ordering);
    void EnableOutOfCore();
    void DisableOutOfCore();
    
    void SetCoefficientEstimationNeededMemory(double);
    void SetMaximumCoefficientEstimationNeededMemory(double);
    void SetIncreaseCoefficientEstimationNeededMemory(double);

    size_t GetMemorySize() const;
    int GetInfoFactorization() const;

    void Factorize(SparseMatrix<T> & mat, bool sym,
		   bool keep_matrix = false);

    void Solve(Vector<T>& x);
  };

}

#define LINALG_FILE_SOLVE_MUMPS_HXX
#endif
