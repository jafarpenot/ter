#ifndef LINALG_FILE_SPARSE_MATRIX_HXX

#include "SparseVector.hxx"

namespace linalg
{
  
  //! Classe de base pour toute matrice 
  template<class T>
  class VirtualMatrix
  {
  protected:
    //! nombre de lignes et colonnes
    int m_, n_;
    
  public:
    VirtualMatrix();
    
    int GetM() const;
    int GetN() const;
    
    virtual void Mlt(const Vector<T>& x, Vector<T>& y) const;
    virtual void MltAdd(const T& alpha, const Vector<T>& x, Vector<T>& y) const;
    
  };


  //! Matrice creuse (chaque ligne est un vecteur creux
  template<class T>
  class SparseMatrix : public VirtualMatrix<T>
  {
  protected:
    //! lignes de la matrice
    Vector<SparseVector<T> > val_;
    
  public:
    SparseMatrix();
    SparseMatrix(int m, int n);
    
    void Reallocate(int m, int n);
    void Clear();

    int GetRowSize(int i) const;
    void ReallocateRow(int i, int n);
    void ClearRow(int i);
    
    const SparseVector<T>& GetLine(int i) const;
    SparseVector<T>& GetLine(int i);
      
    int& Index(int i, int j);
    T& Value(int i, int j);
    
    int Index(int i, int j) const;
    const T& Value(int i, int j) const;
    
    const T operator()(int i, int j) const;
    
    void AddInteraction(int i, int j, const T& x);

    void Mlt(const Vector<T>& x, Vector<T>& y) const;
    void MltAdd(const T& alpha, const Vector<T>& x, Vector<T>& y) const;
    
    void AddM(const SparseMatrix<T>& B, SparseMatrix<T>& C) const;
    void MltConst(const T& val, SparseMatrix<T>& B);
      
    void Transpose(SparseMatrix<T>& B) const;
    void MltM(const SparseMatrix<T>& B, SparseMatrix<T>& AB) const;
    
    void ApplySsor(const Vector<T>& b, const Vector<T>& diag,
		   Vector<T>& x, double omega) const;

    void WriteText(const string&) const;
    
  };

  template<class T>
  ostream& operator<<(ostream& out, const SparseMatrix<T>& A);
  
}

#define LINALG_FILE_SPARSE_MATRIX_HXX
#endif
