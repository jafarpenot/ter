#ifndef LINALG_FILE_COMMON_OUTPUT_HXX

namespace linalg
{

  void WriteMedit(const string& nom, const Vector<complex<double> >& u, int type);
  
  void WriteBinaryInteger(Vector<int>&, ostream& file_out,
			  bool with_size = true, bool swap = false);

  template<class T, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<T, Allocator>& output_vector,
				ostream& file_out, bool double_prec = false,
				bool with_size = true, bool swap = false);
  
  double swapEndian(double d);
  float swapEndian(float d);
  int swapEndian(int d);

}

#define LINALG_FILE_COMMON_OUTPUT_HXX
#endif
