#ifndef LINALG_FILE_TINY_VECTOR_INLINE_CXX

#include "TinyVector.hxx"
#include "TinyVectorExpressionInline.cxx"

namespace linalg
{
  
  /**************
   * TinyVector *
   **************/
  
  
  //! all components are set to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector()
  {
    this->Zero();
  }
  
  
  //! first component is initialized to a, others to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector(const T& a)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::FillGreater(*this, T(0));
  }
  
  
  //! first two components are initialized to a and b, others to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector(const T& a, const T& b)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(*this, b);
    TinyVectorLoop<2*(m>2)>::FillGreater(*this, T(0));
  }
  
  
  //! first three components are initialized to a, b and c, others to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector(const T& a, const T& b ,const T& c)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(*this, b);
    TinyVectorLoop<2*(m>2)>::InitValue(*this, c);
    TinyVectorLoop<3*(m>3)>::FillGreater(*this, T(0));
  }

  
  //! first four components are initialized to a, b, c, and d, others to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector(const T& a, const T& b ,const T& c, const T& d)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(*this, b);
    TinyVectorLoop<2*(m>2)>::InitValue(*this, c);
    TinyVectorLoop<3*(m>3)>::InitValue(*this, d);
    TinyVectorLoop<4*(m>4)>::FillGreater(*this, T(0));
  }


  //! first four components are initialized to a, b, c, and d, others to 0
  template <class T, int m>
  inline TinyVector<T, m>::TinyVector(const T& a, const T& b ,const T& c, const T& d, const T& e)
  {
    this->data_[0] = a;
    TinyVectorLoop<(m>1)>::InitValue(*this, b);
    TinyVectorLoop<2*(m>2)>::InitValue(*this, c);
    TinyVectorLoop<3*(m>3)>::InitValue(*this, d);
    TinyVectorLoop<4*(m>4)>::InitValue(*this, e);
    TinyVectorLoop<5*(m>5)>::FillGreater(*this, T(0));
  }
  

  //! initializing a vector with an expression
  template <class T, int m> template<class E>
  inline TinyVector<T, m>::TinyVector(const TinyVectorExpression<T, m, E> & u)
  {
    TinyVectorLoop<m>::Copy(u, *this);
  }


  //! sets first component to a
  template <class T, int m>
  inline void TinyVector<T, m>::Init(const T& a) 
  {
#ifdef LINALG_DEBUG
    if (m < 1)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 1,")+
		       "but is equal to " + to_string(m) + ".");
#endif

    this->data_[0] = a;
  }

  
  //! sets two first components of the vector
  template <class T, int m>
  inline void TinyVector<T, m>::Init(const T& a, const T& b) 
  {
#ifdef LINALG_DEBUG
    if (m < 2)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 2,")+
		       "but is equal to " + to_string(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
  }
  
  
  //! sets three first components of the vector
  template <class T, int m>
  inline void TinyVector<T, m>::Init(const T& a, const T& b, const T& c) 
  {
#ifdef LINALG_DEBUG
    if (m < 3)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 3,")+
		       "but is equal to " + to_string(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
  }
  
  
  //! sets four first components of the vector
  template <class T, int m>
  inline void TinyVector<T, m>::Init(const T& a, const T& b, const T& c, const T& d) 
  {
#ifdef LINALG_DEBUG
    if (m < 4)
      throw WrongIndex("TinyVector::TinyVector",
                       string("size should be greater than 4,")
		       +"but is equal to " + to_string(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
    this->data_[3] = d;
  }
  

  //! sets five first components of the vector
  template <class T, int m>
  inline void TinyVector<T, m>::Init(const T& a, const T& b,
				     const T& c, const T& d, const T& e) 
  {
#ifdef LINALG_DEBUG
    if (m < 5)
      throw WrongIndex("TinyVector::TinyVector",
		       string("size should be greater than 5,")+
		       "but is equal to " + to_string(m) + ".");
#endif

    this->data_[0] = a;
    this->data_[1] = b;
    this->data_[2] = c;
    this->data_[3] = d;
    this->data_[4] = e;
  }
  

  //! Returns the number of stored elements.
  template<class T, int m_>
  inline int TinyVector<T, m_>::GetM() const
  {
    return m_;
  }
  
  
  //! Returns the number of stored elements.
  template<class T, int m_>
  inline int TinyVector<T, m_>::GetSize() const
  {
    return m_;
  }


  //! returns x(i)
  template <class T, int m_>
  inline const T& TinyVector<T, m_>::operator()(int i) const
  {
#ifdef LINALG_DEBUG
    if ((i < 0) || (i >= m_))
      throw WrongIndex("TinyVector::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_string(int(m_)-1) + "], but is equal to "
		       + to_string(i) + ".");
#endif

    return this->data_[i];
  }
  
  
  //! returns x(i)
  template <class T, int m_>
  inline T& TinyVector<T, m_>::operator()(int i)
  {
#ifdef LINALG_DEBUG
    if ((i < 0) || (i >= m_))
      throw WrongIndex("TinyVector::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_string(int(m_)-1) + "], but is equal to "
		       + to_string(i) + ".");
#endif
    
    return this->data_[i];
  }
  
  
  //! Sets all components to 0
  template <class T, int m>
  inline void TinyVector<T, m>::Zero()
  {
    TinyVectorLoop<m>::Zero(*this);
  } 
  

  //! Sets all components to a given value
  template<class T, int m>
  inline void TinyVector<T, m>::Fill(const T& x)
  {
    TinyVectorLoop<m>::Fill(*this, x);
  }
  
  
  //! Writes the current vector in output stream out without brackets ()
  template<class T, int m>
  inline void TinyVector<T, m>::WriteText(ostream& out)
  {
    TinyVectorLoop<m>::Print(*this, out);
  }


  //! *this = u
  template <class T, int m> template<class T1, class E>
  inline TinyVector<T, m>& TinyVector<T, m>::
  operator =(const TinyVectorExpression<T1, m, E> & u)
  {
    TinyVectorLoop<m>::Copy(u, *this);
    
    return *this;
  }


  //! *this = *this + u
  template <class T, int m> template<class T1, class E>
  inline TinyVector<T, m>& TinyVector<T, m>::
  operator +=(const TinyVectorExpression<T1, m, E> & u)
  {
    TinyVectorLoop<m>::AddCopy(u, *this);
    
    return *this;
  }


  //! *this = *this - u
  template <class T, int m> template<class T1, class E>
  inline TinyVector<T, m>& TinyVector<T, m>::
  operator -=(const TinyVectorExpression<T1, m, E> & u)
  {
    TinyVectorLoop<m>::DiffCopy(u, *this);
    
    return *this;
  }
  
  
  //! sets all components of the vector to a given value
  template<class T, int m>
  inline TinyVector<T, m>& TinyVector<T, m>::operator =(const T& x)
  {
    this->Fill(x);
    return *this;
  }
  
  
  //! operation this = this*a where a is scalar
  template <class T, int m> template<class T1>
  inline TinyVector<T, m> & TinyVector<T, m>::operator *=( const T1& a )
  {
    TinyVectorLoop<m>::Mlt(a, *this);    
    return *this;
  }
  

  //! operation this = this*a where a is scalar
  template <class T, int m> template<class T1>
  inline TinyVector<T, m> & TinyVector<T, m>::operator /=( const T1& a )
  {
    T1 one(1);
    TinyVectorLoop<m>::Mlt(one/a, *this);    
    return *this;
  }

  
  /*************
   * Operators *
   *************/


  //! returns true if u == v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator==(const TinyVectorExpression<T, m, E1>& u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::Equal(u, v);
  }
  
  
  //! returns false if u == v, true otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator!=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::Equal(u, v);
  }
  
    
  //! returns true if u <= v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator<=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::LessOrEqual(u, v);
  }

  //! returns *this < u
  template<class T, int m, class E1, class E2>
  inline bool operator<(const TinyVectorExpression<T, m, E1> & u,
			const TinyVectorExpression<T, m, E2>& v)
  {
    return TinyVectorLoop<m>::Less(u, v);
  }
  
    
  //! returns true if u >= v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator>=(const TinyVectorExpression<T, m, E1> & u,
			 const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::Less(u, v);
  }

    
  //! returns true if u > v, false otherwise
  template<class T, int m, class E1, class E2>
  inline bool operator>(const TinyVectorExpression<T, m, E1> & u,
			const TinyVectorExpression<T, m, E2>& v)
  {
    return !TinyVectorLoop<m>::LessOrEqual(u, v);
  }


  //! displays V
  template <class T, int m, class E>
  inline ostream& operator <<(ostream& out, const TinyVectorExpression<T, m, E>& V)
  {
    out<<"(";
    return TinyVectorLoop<m>::WriteText(V, out);
  }
  
  
  //! reads V
  template <class T, int m>
  inline istream& operator >>(istream& in, TinyVector<T,m>& V)
  {
    return TinyVectorLoop<m>::ReadText(V, in);
  }
    
    
  /******************
   * TinyVectorLoop *
   ******************/


  //! Sets all components to 0
  template<int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::Zero(TinyVector<T0, m>& x)
  {
    x(n-1) = 0;
    TinyVectorLoop<n-1>::Zero(x);
  }
  

  //! Fills with a in ascendant order (instead of descendant order for Fill)
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::FillGreater(TinyVector<T0, m>& x, const T0& a)
  {
    x.data_[n] = a;
    TinyVectorLoop<(n+1)*(m>n+1)>::FillGreater(x, a);
  }

  
  //! initializes a single value
  template <int n> template<int m, class T0>
  inline void TinyVectorLoop<n>::InitValue(TinyVector<T0, m>& x, const T0& a)
  {
    x.data_[n] = a;
  }


  //! Sets all components to a given value
  template <int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::Fill(TinyVector<T0, m>& x, const T1& a)
  {
    x.data_[n-1] = a;
    TinyVectorLoop<n-1>::Fill(x, a);
  }
  

  //! y = x
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) = x(n-1);
    TinyVectorLoop<n-1>::Copy(x, y);
  }
  
  
  //! y += x
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) += x(n-1);
    TinyVectorLoop<n-1>::AddCopy(x, y); 
  }
  

  //! y -= x  
  template<int n> template<int m, class T1, class E, class T0>
  inline void TinyVectorLoop<n>::DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y)
  {
    y(n-1) -= x(n-1);
    TinyVectorLoop<n-1>::DiffCopy(x, y); 
  }
  

  //! returns x == y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::
  Equal(const TinyVectorExpression<T0, m, E0>& x,
	const TinyVectorExpression<T1, m, E1>& y)
  {
    if (abs(x(m-n)-y(m-n)) > TinyVector<T0, m>::threshold)
      return false;
    
    return TinyVectorLoop<n-1>::Equal(x, y);
  }

  
  //! returns x <= y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
                                             const TinyVectorExpression<T1, m, E1>& y)
  {
    if (x(m-n) < (y(m-n)-TinyVector<T0, m>::threshold) )
      return true;
    else if ( x(m-n) > (y(m-n) + TinyVector<T0, m>::threshold))
      return false;
    
    return TinyVectorLoop<n-1>::LessOrEqual(x, y);
  }
  

  //! returns x < y
  template <int n> template<int m, class T0, class T1, class E0, class E1>
  inline bool TinyVectorLoop<n>::Less(const TinyVectorExpression<T0, m, E0>& x,
				      const TinyVectorExpression<T1, m, E1>& y)
  {
    if (x(m-n) < (y(m-n)-TinyVector<T0, m>::threshold) )
      return true;
    else if (x(m-n) > (y(m-n) + TinyVector<T0, m>::threshold))
      return false;

    return TinyVectorLoop<n-1>::Less(x, y);
  }
  

  //! Multiplies x by alpha
  template <int n> template<int m, class T0, class T1>
  inline void TinyVectorLoop<n>::Mlt(const T1& alpha, TinyVector<T0, m>& x)
  {
    x(n-1) *= alpha;
    TinyVectorLoop<n-1>::Mlt(alpha, x);
  }
  

  //! scal = u.v
  template <int n> template<int m, class T0, class T1, class T2, class E0, class E1>
  inline void TinyVectorLoop<n>::
  DotProd(const TinyVectorExpression<T0, m, E0>& u,
	  const TinyVectorExpression<T1, m, E1>& v, T2& scal)
  {
    scal += u(n)*v(n);  
    TinyVectorLoop<n-1>::DotProd(u, v, scal);
  }
  
  
  //! computes scal = \sum |u_i|^2
  template <int n> template<int m, class T0, class T1, class E0>
  inline void TinyVectorLoop<n>::AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal)
  {
    scal += absSquare(u(n));
    TinyVectorLoop<n-1>::AbsSquare(u, scal);
  }


  //! displays V without brackets
  template<int n> template<int m, class T, class E>
  inline void TinyVectorLoop<n>::Print(const TinyVectorExpression<T, m, E>& V, ostream& out)
  {
    out << V(m-n) << " ";      
    return TinyVectorLoop<n-1>::Print(V, out);
  }
  
  
  //! displays V
  template<int n> template<int m, class T, class E>
  inline ostream& TinyVectorLoop<n>::WriteText(const TinyVectorExpression<T, m, E>& V, ostream& out)
  {
    if (n==1)
      out << V(m-n) << ')';  
    else
      out << V(m-n) << ", ";  
    
    return TinyVectorLoop<n-1>::WriteText(V, out);
  }
  

  //! reads V
  template<int n> template<int m, class T>
  inline istream& TinyVectorLoop<n>::ReadText(TinyVector<T, m>& V, istream& in)
  {
    in >> V(m-n);  
    return TinyVectorLoop<n-1>::ReadText(V, in);
  }
  

  /**********************
   * Exterior functions *
   **********************/
  

  template<class T>
  T absSquare(const T& x)
  {
    return x*x;
  }

  
  template<class T>
  T absSquare(const complex<T>& x)
  {
    return real(x)*real(x) + imag(x)*imag(x);
  }
  
  
  //! returns scalar product (u,v)
  template<class T, int m, class E1, class E2>
  inline T DotProd(const TinyVectorExpression<T, m, E1> & u,
		   const TinyVectorExpression<T, m, E2> & v)
  {
    T scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }

  
  //! returns scalar product (u,v)
  template<class T, int m, class E1, class E2>
  inline complex<T> DotProd(const TinyVectorExpression<T, m, E1> & u,
			    const TinyVectorExpression<complex<T>, m, E2> & v)
  {
    complex<T> scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }
  
  
  //! returns (u,v) where u is complex and v real
  template<class T, int m, class E1, class E2>
  inline complex<T> DotProd(const TinyVectorExpression<complex<T>, m, E1> & u,
			    const TinyVectorExpression<T, m, E2> & v)
  {
    complex<T> scal = u(0)*v(0);
    TinyVectorLoop<m-1>::DotProd(u, v, scal);
    
    return scal;
  }
  
  
  //! returns || u ||^2
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal AbsSquare(const TinyVectorExpression<T, m, E> & u)
  {
    typename ClassComplexType<T>::Treal scal = absSquare(u(0));
    TinyVectorLoop<m-1>::AbsSquare(u, scal);
    
    return scal;
  }


  //! returns || v ||_2
  template<class T, int m, class E>
  inline typename ClassComplexType<T>::Treal Norm2(const TinyVectorExpression<T, m, E> & p)
  {
    return sqrt(AbsSquare(p));
  }
  
} // end namespace

#define LINALG_FILE_TINY_VECTOR_INLINE_CXX
#endif
