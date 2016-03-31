#ifndef LINALG_FILE_VECTOR_EXPRESSION_CXX

namespace linalg
{
  
  //! returns the size of the associated vector
  template<class T, class E>
  inline size_t VectorExpression<T, E>::GetSize() const
  {
    return static_cast<const E&>(*this).GetSize();
  }

  
  //! returns the element i of expression
  template<class T, class E>
  inline const T VectorExpression<T, E>::operator()(size_t i) const
  {
    return static_cast<const E&>(*this)(i);
  }
    
  
  //! Constructor with two expressions u and v
  template<class T, class E1, class E2>
  inline VectorDifference<T, E1, E2>::
  VectorDifference(const VectorExpression<T, E1>& u, const VectorExpression<T, E2>& v)
    : u_(u), v_(v)
  {
#ifdef LINALG_DEBUG
    if (u_.GetSize() != v_.GetSize())
      throw WrongIndex("VectorDifference",
		       string("Cannot subtract u and v because the sizes ")
		       +to_string(u.GetSize()) + " and " + to_string(v.GetSize())
		       + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline size_t VectorDifference<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of the difference
  template<class T, class E1, class E2>
  inline const T VectorDifference<T, E1, E2>::operator()(size_t i) const
  {
    return u_(i) - v_(i);
  }
  
  
  //! Constructor with two expressions u and v
  template<class T, class E1, class E2>
  inline VectorSum<T, E1, E2>::VectorSum(const VectorExpression<T, E1>& u, const VectorExpression<T, E2>& v)
    : u_(u), v_(v)
  {
#ifdef LINALG_DEBUG
    if (u_.GetSize() != v_.GetSize())
      throw WrongIndex("VectorDifference",
		       string("Cannot Add u and v because the sizes ")
		       +to_string(u.GetSize()) + " and " + to_string(v.GetSize())
		       + " are different");
#endif
  }
  
  
  //! returns the size of the associated vectors
  template<class T, class E1, class E2>
  inline size_t VectorSum<T, E1, E2>::GetSize() const
  {
    return u_.GetSize();
  }
  

  //! returns the i-th element of the sum  
  template<class T, class E1, class E2>
  inline const T VectorSum<T, E1, E2>::operator()(size_t i) const
  {
    return u_(i) + v_(i);
  }
  

  //! Constructor with a scalar and an expression
  template<class T, class E>
  inline VectorScaled<T, E>::VectorScaled(const T& alpha,
					  const VectorExpression<T, E>& u)
    : alpha_(alpha), u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class E>
  inline size_t VectorScaled<T, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of alpha*u
  template<class T, class E>
  inline const T VectorScaled<T, E>::operator()(size_t i) const
  {
    return alpha_*u_(i);
  }
  

  //! Constructor -u with an expression u
  template<class T, class E>
  inline VectorOpposite<T, E>::VectorOpposite(const VectorExpression<T, E>& u)
    : u_(u)
  {
  }
  
    
  //! returns the size of the associated vector
  template<class T, class E>
  inline size_t VectorOpposite<T, E>::GetSize() const
  {
    return u_.GetSize();
  }
  
  
  //! returns the i-th element of -u
  template<class T, class E>
  inline const T VectorOpposite<T, E>::operator()(size_t i) const
  {
    return -u_(i);
  }

  
  /*************
   * Operators *
   *************/
  
  
  //! returns u+v
  template<class T, class E1, class E2>
  inline const VectorSum<T, E1, E2> operator+(const VectorExpression<T, E1>& u,
					      const VectorExpression<T, E2>& v)
  {
    return VectorSum<T, E1, E2>(u, v);
  }
  

  //! returns u-v
  template<class T, class E1, class E2>
  inline const VectorDifference<T, E1, E2> operator-(const VectorExpression<T, E1>& u,
					       const VectorExpression<T, E2>& v)
  {
    return VectorDifference<T, E1, E2>(u, v);
  }

  
  //! returns alpha*u
  template<class T, class E>
  inline const VectorScaled<T, E> operator*(const T& alpha,
				      const VectorExpression<T, E>& u)
  {
    return VectorScaled<T, E>(alpha, u);
  }

  
  //! returns u*alpha
  template<class T, class E>
  inline const VectorScaled<T, E> operator*(const VectorExpression<T, E>& u, const T& alpha)
  {
    return VectorScaled<T, E>(alpha, u);
  }


  //! returns -u
  template<class T, class E>
  inline const VectorOpposite<T, E> operator-(const VectorExpression<T, E>& u)
  {
    return VectorOpposite<T, E>(u);
  }
  
}

#define LINALG_FILE_VECTOR_EXPRESSION_HXX
#endif
