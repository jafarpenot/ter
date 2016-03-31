#ifndef LINALG_FILE_COMMON_OUTPUT_CXX

#include "CommonOutput.hxx"

namespace linalg
{
  
  //! writes u in a file readable by Medit (.bb file)
  void WriteMedit(const string& nom, const Vector<complex<double> >& u, int type)
  {
    ofstream file_out(nom.data());
    file_out.precision(8);
    file_out << "2 1 " << u.GetM() << " 2 \n";
    if (type == 0)
      for (int i = 0; i < u.GetM(); i++)
	file_out << real(u(i)) << '\n';
    else if (type == 1)
      for (int i = 0; i < u.GetM(); i++)
	file_out << imag(u(i)) << '\n';
    else
      for (int i = 0; i < u.GetM(); i++)
	file_out << abs(u(i)) << '\n';
  }
  

  //! writes a vector of integers in binary with eventual swapping
  void WriteBinaryInteger(Vector<int>& output_vector,
			  ostream& file_out, bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (swap)
      for (int i = 0; i < N; i++)
	output_vector(i) = swapEndian(output_vector(i));
    
    output_vector.Write(file_out, with_size);
  }
  
  
  //! writing a real vector in single or double precision
  template<class T, class Allocator>
  void WriteBinaryDoubleOrFloat(const Vector<T, Allocator>& output_vector,
				ostream& file_out, bool double_prec, bool with_size, bool swap)
  {
    int N = output_vector.GetM();
    if (double_prec)
      {
	Vector<double> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(double(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = double(output_vector(i));
	
	output.Write(file_out, with_size);
      }
    else
      {
	Vector<float> output(N);
	if (swap)
	  for (int i = 0; i < N; i++)
	    output(i) = swapEndian(float(output_vector(i)));
	else
	  for (int i = 0; i < N; i++)
	    output(i) = float(output_vector(i));
	
	output.Write(file_out, with_size);
      }
  }


  //! swapping octets of a real number (conversion big endian <-> little endian)
  double swapEndian(double d)
  {
    double a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];
    
    return a;
  }


  //! swapping octets of a real number (conversion big endian <-> little endian)
  float swapEndian(float d)
  {
    float a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    
    return a;
  }


  //! swapping octets of an integer (conversion big endian <-> little endian)
  int swapEndian(int d)
  {
    int a;
    char *dst = reinterpret_cast<char*>(&a);
    char *src = reinterpret_cast<char*>(&d);
    
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    
    return a;
  }
  
}

#define LINALG_FILE_COMMON_OUTPUT_CXX
#endif

