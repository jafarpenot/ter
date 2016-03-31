#ifndef LINALG_FILE_VECTOR_EXPRESSION_HXX

namespace linalg
{
  
  //! Expression between vectors
  template<class T, class E>
  class VectorExpression
  {
  public:
    size_t GetSize() const;
    const T operator()(size_t) const;
    
    inline operator E&() { return static_cast<E&>(*this); }
    inline operator E const&() const { return static_cast<const E&>(*this); }
    
  };


  //! Difference between two expressions
  template<class T, class E1, class E2>
  class VectorDifference : public VectorExpression<T, VectorDifference<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    VectorDifference(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);
    
    size_t GetSize() const;
    const T operator()(size_t) const;
    
  };


  //! Sum between two expressions
  template<class T, class E1, class E2>
  class VectorSum : public VectorExpression<T, VectorSum<T, E1, E2> >
  {
    const E1& u_;
    const E2& v_;
    
  public:
    VectorSum(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);
    
    size_t GetSize() const;
    const T operator()(size_t) const;
    
  };


  //! Multiplication between a scalar and a vector
  template<class T, class E>
  class VectorScaled : public VectorExpression<T, VectorScaled<T, E> >
  {
    const T alpha_;
    const E& u_;
    
  public:
    VectorScaled(const T&, const VectorExpression<T, E>&);
    
    size_t GetSize() const;
    const T operator()(size_t) const;
    
  };
  

  //! Opposite of a vector
  template<class T, class E>
  class VectorOpposite : public VectorExpression<T, VectorOpposite<T, E> >
  {
    const E& u_;
    
  public:
    VectorOpposite(const VectorExpression<T, E>&);
    
    size_t GetSize() const;
    const T operator()(size_t) const;
    
  };

  
  /*************
   * Operators *
   *************/
  
  
  template<class T, class E1, class E2>
  const VectorSum<T, E1, E2> operator+(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);


  template<class T, class E1, class E2>
  const VectorDifference<T, E1, E2> operator-(const VectorExpression<T, E1>&, const VectorExpression<T, E2>&);


  template<class T, class E>
  const VectorScaled<T, E> operator*(const T&, const VectorExpression<T, E>&);

  
  template<class T, class E>
  const VectorScaled<T, E> operator*(const VectorExpression<T, E>&, const T&);


  template<class T, class E>
  const VectorOpposite<T, E> operator-(const VectorExpression<T, E>&);
  
}

#define LINALG_FILE_VECTOR_EXPRESSION_HXX
#endif
