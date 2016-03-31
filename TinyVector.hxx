#ifndef LINALG_FILE_TINY_VECTOR_HXX

#include "TinyVectorExpression.hxx"

namespace linalg
{  
 
  template<int n>
  class TinyVectorLoop;

  //! Class storing a tiny vector whose size is known at compilation time
  template <class T, int m_>
  class TinyVector : public TinyVectorExpression<T, m_, TinyVector<T, m_> >
  {
  public :
    typedef T value_type;
    
    template<int n>
    friend class TinyVectorLoop;
    
    //! threshold used to compare points
    static constexpr double threshold = 1e-12;
    
  protected: 
    // pointer data_;
    T data_[m_];
    
  public :
    TinyVector();
    
    explicit TinyVector(const T& a);
    explicit TinyVector(const T& a, const T& b);
    TinyVector(const T& a, const T& b, const T& c);
    TinyVector(const T& a, const T& b, const T& c, const T& d);
    TinyVector(const T& a, const T& b, const T& c, const T& d, const T& e);
    
    
    template<class E>
    TinyVector(const TinyVectorExpression<T, m_, E> & u);

    void Init(const T & a);
    void Init(const T & a, const T& b);
    void Init(const T & a, const T& b, const T& c);
    void Init(const T & a, const T& b, const T& c, const T& d);
    void Init(const T & a, const T& b, const T& c, const T& d, const T& e);
    
    // Basic Functions
    int GetM() const;
    int GetSize() const;
    
    // Element access and affection
    const T& operator()(int i) const;
    T& operator()(int i);

    // Convenient Functions
    void Zero() ;  // we put zero on coordinates
    void Fill(const T& a); // all values are set to a
    void WriteText(ostream& out); // the vector is written in the stream without brackets and commas
    
    // operation this = u
    template<class T1, class E>
    TinyVector<T, m_>& operator =(const TinyVectorExpression<T1, m_, E> & u);

    // this += u
    template<class T1, class E>
    TinyVector<T, m_>& operator +=(const TinyVectorExpression<T1, m_, E> & u);

    // this -= u
    template<class T1, class E>
    TinyVector<T, m_>& operator -=(const TinyVectorExpression<T1, m_, E> & u);

    // fills vector with x
    TinyVector<T, m_>& operator =(const T& x);

    // operation this = a*this where a is scalar
    template<class T1>
    TinyVector<T, m_>& operator *=( const T1& a );
    
    // operation this = 1/a*this where a is scalar
    template<class T1>
    TinyVector<T, m_>& operator /=( const T1& a );
    
  };


  template<class T, int m, class E1, class E2>
  bool operator==(const TinyVectorExpression<T, m, E1>& u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator!=(const TinyVectorExpression<T, m, E1> & u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator<=(const TinyVectorExpression<T, m, E1> & u,
		  const TinyVectorExpression<T, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator<(const TinyVectorExpression<T, m, E1> & u,
		 const TinyVectorExpression<T, m, E2>& v);
    
  template<class T, int m, class E1, class E2>
  bool operator>=(const TinyVectorExpression<int, m, E1> & u,
		  const TinyVectorExpression<int, m, E2>& v);

  template<class T, int m, class E1, class E2>
  bool operator>(const TinyVectorExpression<T, m, E1> & u,
		 const TinyVectorExpression<T, m, E2>& v);
  
  template <class T, int m, class E>
  ostream& operator <<(ostream& out, const TinyVectorExpression<T, m, E>& V);

  template <class T, int m>
  istream& operator >>(istream& out, TinyVector<T, m>& V);
  
  //! class used for unrolling loops for operators/functions of TinyVector
  template<int n>
  class TinyVectorLoop
  {
  public :
    template<int m, class T0>
    static void Zero(TinyVector<T0, m>& x);
    
    template<int m, class T0>
    static void FillGreater(TinyVector<T0, m>& x, const T0& a);

    template<int m, class T0>
    static void InitValue(TinyVector<T0, m>& x, const T0& a);

    template<int m, class T0, class T1>
    static void Fill(TinyVector<T0, m>& x, const T1& a);

    template<int m, class T1, class E, class T0>
    static void Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T1, class E, class T0>
    static void AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T1, class E, class T0>
    static void DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y);

    template<int m, class T0, class T1, class E0, class E1>
    static bool Equal(const TinyVectorExpression<T0, m, E0>& x,
		      const TinyVectorExpression<T1, m, E1>& y);
    
    template<int m, class T0, class T1, class E0, class E1>
    static bool LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
			    const TinyVectorExpression<T1, m, E1>& y);
    
    template<int m, class T0, class T1, class E0, class E1>
    static bool Less(const TinyVectorExpression<T0, m, E0>& x,
		     const TinyVectorExpression<T1, m, E1>& y);

    template<int m, class T0, class T1>
    static void Mlt(const T1& alpha, TinyVector<T0, m>& x);
    
    template<int m, class T0, class T1, class T2, class E0, class E1>
    static void DotProd(const TinyVectorExpression<T0, m, E0>& u,
			const TinyVectorExpression<T1, m, E1>& v, T2& scal);
    
    template<int m, class T0, class T1, class E0>
    static void AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal);
    
    template<int m, class T, class E>
    static void Print(const TinyVectorExpression<T, m, E>& V, ostream& out);

    template<int m, class T, class E>
    static ostream& WriteText(const TinyVectorExpression<T, m, E>& V, ostream& out);

    template<int m, class T>
    static istream& ReadText(TinyVector<T, m>& V, istream& in);
    
  };
  
  //! Class used to terminate the loop
  template<>
  class TinyVectorLoop<0>
  {
  public :
    template<int m, class T0>
    static inline void Zero(TinyVector<T0, m>& x) {}

    template<int m, class T0>
    static inline bool IsZero(const TinyVector<T0, m>& x) { return true; }
    
    template<int m, class T0>
    static inline void FillGreater(TinyVector<T0, m>& x, const T0& a){}

    template<int m, class T0>
    static inline void InitValue(TinyVector<T0, m>& x, const T0& a){}

    template<int m, class T0, class T1>
    static inline void Fill(TinyVector<T0, m>& x, const T1& a){}

    template<int m, class T1, class E, class T0>
    static inline void Copy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 

    template<int m, class T1, class E, class T0>
    static inline void AddCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 

    template<int m, class T1, class E, class T0>
    static inline void DiffCopy(const TinyVectorExpression<T1, m, E>& x, TinyVector<T0, m>& y) {} 
    
    template<int m, class T0, class T1, class E0, class E1>
    static inline bool Equal(const TinyVectorExpression<T0, m, E0>& x,
			     const TinyVectorExpression<T1, m, E1>& y)
    {
      return true;
    }

    template<int m, class T0, class T1, class E0, class E1>
    static inline bool LessOrEqual(const TinyVectorExpression<T0, m, E0>& x,
				   const TinyVectorExpression<T1, m, E1>& y)
    {
      return true;
    }

    template<int m, class T0, class T1, class E0, class E1>
    static inline bool Less(const TinyVectorExpression<T0, m, E0>& x,
			    const TinyVectorExpression<T1, m, E1>& y)
    {
      return false;
    }

    template<int m, class T0, class T1>
    static inline void Mlt(const T1& alpha, TinyVector<T0, m>& x) {}
    
    template<int m, class T0, class T1, class T2, class E0, class E1>
    static inline void DotProd(const TinyVectorExpression<T0, m, E0>& u,
			       const TinyVectorExpression<T1, m, E1>& v, T2& scal) {}

    template<int m, class T0, class T1, class E0>
    static inline void AbsSquare(const TinyVectorExpression<T0, m, E0>& u, T1& scal) {}
    
    template<int m, class T, class E>
    static inline void Print(const TinyVectorExpression<T, m, E>& V, ostream& out) { }

    template<int m, class T, class E>
    static inline ostream& WriteText(const TinyVectorExpression<T, m, E>& V, ostream& out) { return out; }

    template<int m, class T>
    static inline istream& ReadText(const TinyVector<T, m>& V, istream& in) { return in; }
    
  };

  template<class T>
  T absSquare(const T& x);
  
  template<class T>
  T absSquare(const complex<T>& x);

  // returns (u,v)
  template<class T, int m, class E1, class E2>
  T DotProd(const TinyVectorExpression<T, m, E1> & u,
	    const TinyVectorExpression<T, m, E2> & v);
    
  // returns (u,v) where v is complex and u real
  template<class T, int m, class E1, class E2>
  complex<T> DotProd(const TinyVectorExpression<T, m, E1> & u,
		     const TinyVectorExpression<complex<T>, m, E2> & v);
  
  template<class T, int m, class E1, class E2>
  complex<T> DotProd(const TinyVectorExpression<complex<T>, m, E1> & u,
		     const TinyVectorExpression<T, m, E2> & v);

  template<class T>
  class ClassComplexType
  {
  public:
    typedef T Treal;
  };


  template<class T>
  class ClassComplexType<complex<T> >
  {
  public:
    typedef T Treal;
  };

  // returns || u ||^2 
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal AbsSquare(const TinyVectorExpression<T, m, E> & u);

  // returns || u || 
  template<class T, int m, class E>
  typename ClassComplexType<T>::Treal Norm2(const TinyVectorExpression<T, m, E> & v);
  
} // end namespace

#define LINALG_FILE_TINY_VECTOR_HXX
#endif
