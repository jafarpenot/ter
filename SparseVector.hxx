#ifndef LINALG_FILE_SPARSE_VECTOR_HXX

namespace linalg
{
  
  //! classe stockant un vecteur creux
  template<class T>
  class SparseVector
  {
  private:
    Vector<T> values;
    Vector<int> index;
    
  public:
    int GetM() const;
    
    void Reallocate(int n);
    void Resize(int n);
    void Clear();
    
    int& Index(int j);
    T& Value(int j);
    
    int Index(int j) const;
    const T& Value(int j) const;

    const T operator()(int i) const;
    
    void AddInteraction(int j, const T& val);
    
  };

}

#define LINALG_FILE_SPARSE_VECTOR_HXX
#endif
