#ifndef LINALG_FILE_SPARSE_MATRIX_CXX

#include "SparseMatrix.hxx"

namespace linalg
{
  
  //! constructeur par defaut
  template<class T>
  VirtualMatrix<T>::VirtualMatrix()
  {
    this->m_ = 0;
    this->n_ = 0;
  }

  //! Renvoie le nombre de lignes de la matrice
  template<class T>
  inline int VirtualMatrix<T>::GetM() const
  {
    return m_;
  }
    
  //! Renvoie le nombre de colonnes de la matrice
  template<class T>
  inline int VirtualMatrix<T>::GetN() const
  {
    return n_;
  }
    
  //! Effectue le produit matrice-vecteur y = A x
  template<class T>
  void VirtualMatrix<T>::Mlt(const Vector<T>& x, Vector<T>& y) const
  {
    cout << "Not implemented for any kind of matrix " << endl;
    abort();
  }
  
  //! Effectue le produit matrice-vecteur y = y + alpha A x
  template<class T>
  void VirtualMatrix<T>::MltAdd(const T& alpha, const Vector<T>& x, Vector<T>& y) const
  {
    cout << "Not implemented for any kind of matrix " << endl;
    abort();
  }
  
  //! Constructeur par defaut
  template<class T>
  SparseMatrix<T>::SparseMatrix()
  {
  }
  
  //! Constructeur avec le nombre de lignes et colonnes
  template<class T>
  SparseMatrix<T>::SparseMatrix(int m, int n)
  {
    this->m_ = m;
    this->n_ = n;
    val_.Reallocate(m);
  }

  //! Change le nombre de lignes et colonnes de la matrice
  template<class T>
  void SparseMatrix<T>::Reallocate(int m, int n)
  {
    this->m_ = m;
    this->n_ = n;
    val_.Reallocate(m);
  }
  
  //! Efface la matrice
  template<class T>
  void SparseMatrix<T>::Clear()
  {
    this->m_ = 0;
    this->n_ = 0;
    val_.Clear();
  }
  
  //! Retourne le nombre d'elements non-nuls de la ligne i
  template<class T>
  inline int SparseMatrix<T>::GetRowSize(int i) const
  {
    return val_(i).GetM();
  }

  //! Change le nombre d'elements non-nuls de la ligne i
  template<class T>
  void SparseMatrix<T>::ReallocateRow(int i, int n)
  {
    val_(i).Reallocate(n);
  }
  
  //! Efface la ligne i
  template<class T>
  void SparseMatrix<T>::ClearRow(int i)
  {
    val_(i).Clear();
  }
    
  //! Renvoie le numero de colonne de l'element non-nul j de la ligne i
  template<class T>
  inline int& SparseMatrix<T>::Index(int i, int j)
  {
    return val_(i).Index(j);
  }
  
  //! Renvoie la valeur de l'element non-nul j de la ligne i
  template<class T>
  inline T& SparseMatrix<T>::Value(int i, int j)
  {
    return val_(i).Value(j);
  }
    
  //! Renvoie le numero de colonne de l'element non-nul j de la ligne i
  template<class T>
  inline int SparseMatrix<T>::Index(int i, int j) const
  {
    return val_(i).Index(j);
  }

  //! Renvoie la valeur de l'element non-nul j de la ligne i  
  template<class T>
  inline const T& SparseMatrix<T>::Value(int i, int j) const
  {
    return val_(i).Value(j);
  }

  //! Renvoie la ligne i : le ieme SparseVector de la matrice
  template<class T>
  inline const SparseVector<T>& SparseMatrix<T>::GetLine(int i) const
  {
    return val_(i);
  }
    
  //! Renvoie la ligne i : le ieme SparseVector de la matrice
  template<class T>
  inline SparseVector<T>& SparseMatrix<T>::GetLine(int i)
  {
      return val_(i);
  }
    
  //! Renvoie A(i, j)
  template<class T>
  inline const T SparseMatrix<T>::operator()(int i, int j) const
  {
    return val_(i)(j);
  }
    
  //! Ajoute x a A(i, j)
  template<class T>
  void SparseMatrix<T>::AddInteraction(int i, int j, const T& x)
  {
#ifdef LINALG_DEBUG
    if ((i < 0) || (i >= this->m_) || (j < 0) || (j >= this->n_))
      throw WrongIndex("SparseMatrix::AddInteraction",
		       string("Trying to add a value to A(")
		       + to_string(i) + ", " + to_string(j) + ")"
		       + string(" but the size of the matrix is ")
		       + to_string(this->m_) + " x " + to_string(this->n_) + ".");
#endif
    
    val_(i).AddInteraction(j, x);
  }

  //! Effectue le produit matrice-vecteur y = A x
  template<class T>
  void SparseMatrix<T>::Mlt(const Vector<T>& x, Vector<T>& y) const
  {
    T val;
    for (int i = 0; i < this->GetM(); i++)
      {
	val = 0;
	for (int j = 0; j < this->GetRowSize(i); j++)
	  val += this->Value(i, j)*x(this->Index(i, j));
	
	y(i) = val;
      }
  }

    

    //! Effectue la somme des matrices A + B  = C
    template<class T>
    void SparseMatrix<T>::AddM(const SparseMatrix<T>& B, SparseMatrix<T>& C) const
    {
        if (this->GetM() != B.GetM() || this->GetN() != B.GetN())
        {
            cout << "les matrices doivent etre de meme taille" << endl;
            abort();
        }
        
        C.Clear();
        C.Reallocate(B.GetM(), B.GetN());
        
        for (int i=0; i < B.GetN(); i++)
        {
            for (int j=0; j<this->GetRowSize(i); j++)
                C.AddInteraction(i, this->Index(i,j), this->Value(i,j));
            
            for (int j=0; j<B.GetRowSize(i); j++)
                C.AddInteraction(i, B.Index(i,j), B.Value(i,j));
            
        }
    }

    //! Multiplie la matrice par val et stocke le resultat dans B
    template<class T>
    void SparseMatrix<T>::MltConst(const T& val, SparseMatrix<T>& B)
    {
        B.Clear();
        B.Reallocate(this->GetM(), this->GetN());
        for (int i=0; i < this->GetN(); i++)
        {
            for (int j=0; j<this->GetRowSize(i); j++)
                B.AddInteraction(i,j, this->Value(i,j)*val);
        }
    }
    
    

    
    //! Effectue le produit matrice-matrice M = A B
    template<class T>
    void SparseMatrix<T>::MltM(const SparseMatrix<T>& B, SparseMatrix<T>& AB) const
    {
        AB.Clear();
        AB.Reallocate(this->GetM(), B.GetN());
        T val;
        for (int k=0; k < B.GetN(); k++)
        {
            for (int i = 0; i < this->GetM(); i++)
            {
                val = 0;
                for (int j = 0; j < this->GetRowSize(i); j++)
                    val += this->Value(i, j)*B(this->Index(i,j), k);
                
                AB.AddInteraction(i, k, val);
            }
        }
    }
    
    //! Transpose la matrice
    template<class T>
    void SparseMatrix<T>::Transpose(SparseMatrix<T>& B) const
    {
        B.Clear();
        B.Reallocate(this->GetN(), this->GetM());
        for (int i = 0; i < this->GetM(); i++)
        {
            for (int j = 0; j < this->GetRowSize(i); j++)
                B.AddInteraction(this->Index(i,j), i, this->Value(i,j));
        }
    }
    
  
  //! Effectue le produit matrice-vecteur y = y + alpha A x
  template<class T>
  void SparseMatrix<T>::MltAdd(const T& alpha, const Vector<T>& x, Vector<T>& y) const
  {
    T val;
    for (int i = 0; i < this->GetM(); i++)
      {
	val = 0;
	for (int j = 0; j < this->GetRowSize(i); j++)
	  val += this->Value(i, j)*x(this->Index(i, j));
	
	y(i) += alpha*val;
      }
  }
  
  
  //! Effectue une iteration de SSOR (Symmetric Successive Over Relaxation)
  template<class T>
  void SparseMatrix<T>::ApplySsor(const Vector<T>& b, const Vector<T>& invDiag,
				  Vector<T>& x, double omega) const
  {
    T val;
    // etape de descente
    for (int i = 0; i < this->GetM(); i++)
      {
	val = b(i);
	for (int j = 0; j < this->GetRowSize(i); j++)
	  val -= this->Value(i, j)*x(this->Index(i, j));
	
	// on met a jour x
	x(i) += omega*val*invDiag(i);
      }
    
    // remontee
    for (int i = this->GetM()-1; i >= 0; i--)
      {
	val = b(i);
	for (int j = 0; j < this->GetRowSize(i); j++)
	  val -= this->Value(i, j)*x(this->Index(i, j));
	
	// on met a jour x
	x(i) += omega*val*invDiag(i);
      }
  }
  

  //! Writes the content of the matrix in a text file
  template<class T>
  void SparseMatrix<T>::WriteText(const string& FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("SparseMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
    
    FileStream << *this;
    
    // end of line to finish the file
    FileStream << '\n';
    
    FileStream.close();    
  }
  
  
  //! Ecrit la matrice A
  template<class T>
  ostream& operator<<(ostream& out, const SparseMatrix<T>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	out << i << " " << A.Index(i, j) << " " << A.Value(i, j) << '\n';
    
    return out;
  }
  
}

#define LINALG_FILE_SPARSE_MATRIX_CXX
#endif
