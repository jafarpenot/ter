#ifndef LINALG_FILE_COCG_HXX

namespace linalg
{

  //! Classe de base pour tout preconditionneur
  template<class T>
  class VirtualPreconditioner
  {
  public:
    // methode a surcharger dans les classes derivees
    virtual void Solve(const Vector<T>& r, Vector<T>& z) = 0;
    
  };


  template <class T>
  void ConjugateGradient(const VirtualMatrix<T>& A, Vector<T>& x, const Vector<T>& b,
                         VirtualPreconditioner<T>& prec, double epsilon = 1e-6, int nb_iter_max = 1000);
  
  //! Preconditioneur identite (M = I)
  template<class T>
  class IdentityPreconditioner : public VirtualPreconditioner<T>
  {
  public:
    void Solve(const Vector<T>& r, Vector<T>& z);
    
  };

}

#define LINALG_FILE_COCG_HXX
#endif
