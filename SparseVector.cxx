#ifndef LINALG_FILE_SPARSE_VECTOR_CXX

#include "SparseVector.hxx"

namespace linalg
{
  
  //! Renvoie le nombre d'elements non-nuls
  template<class T>
  inline int SparseVector<T>::GetM() const
  {
    return values.GetM();
  }
  
  //! Change le nombre d'elements non nuls
  template<class T>
  inline void SparseVector<T>::Reallocate(int n)
  {
    values.Reallocate(n);
    index.Reallocate(n);
  }
  
  //! Change le nombre d'elements non nuls (les anciens elements sont conserves)
  template<class T>
  inline void SparseVector<T>::Resize(int n)
  {
    values.Resize(n);
    index.Resize(n);
  }
  
  //! Efface tous les elements non-nuls stockes
  template<class T>
  inline void SparseVector<T>::Clear()
  {
    values.Clear();
    index.Clear();
  }
  
  //! Retourne une reference du numero de colonne associe a l'element non-nul j
  template<class T>
  inline int& SparseVector<T>::Index(int j)
  {
    return index(j);
  }
  
  //! Retourne une reference de la valeur associee a l'element non-nul j
  template<class T>
  inline T& SparseVector<T>::Value(int j)
  {
    return values(j);
  }
  
  //! Retourne le numero de colonne associe a l'element non-nul j
  template<class T>
  inline int SparseVector<T>::Index(int j) const
  {
    return index(j);
  }
  
  //! Retourne la valeur associee a l'element non-nul j
  template<class T>
  inline const T& SparseVector<T>::Value(int j) const
  {
    return values(j);
  }
  
  //! Rajoute val a la colonne j du vecteur
  template<class T>
  void SparseVector<T>::AddInteraction(int j, const T& val)
  {
    int k = 0;
    while ((k < index.GetM()) && (index(k) < j))
      k++;
    
    if (k < index.GetM())
      {
	// on ajoute si on a trouve le bon numero
	if (index(k) == j)
	  values(k) += val;
	else
	  {
	    // sinon on insere la valeur au bon endroit
	    int n = values.GetM();
	    values.Resize(n+1);
	    index.Resize(n+1);
	    for (int i = n-1; i >= k; i--)
	      {
		values(i+1) = values(i);
		index(i+1) = index(i);
	      }
	    
	    index(k) = j;
	    values(k) = val;
	  }
      }
    else
      {
	// cas ou le numero de colonne j est a la fin
	index.PushBack(j);
	values.PushBack(val);
      }
  }

  //! Renvoie x(j)
  template<class T>
  const T SparseVector<T>::operator()(int j) const
  {
    int k = 0;
    while ((k < index.GetM()) && (index(k) < j))
      k++;
    
    if ((k < index.GetM()) && (index(k) == j))
      return values(k);
    
    return T(0);
  }
  
  //! Imprime le vecteur creux
  template<class T>
  ostream& operator<<(ostream& out, const SparseVector<T>& v)
  {
    for (int i = 0; i < v.GetM(); i++)
      out << v.Index(i) << " " << v.Value(i) << '\n';
    
    return out;
  }
  
}

#define LINALG_FILE_SPARSE_VECTOR_CXX
#endif
