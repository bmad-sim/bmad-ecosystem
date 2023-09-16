#ifndef CONVERTER_TEMPLATES

#include <string>
#include <valarray>
#include <complex>
#include "bmad_std_typedef.h"

//---------------------------------------------------------------------------

template <class T> void operator<< (valarray<T>& arr, const T* ptr) {
  unsigned int n = arr.size();
  for (unsigned int i = 0; i < n; i++) arr[i] = ptr[i];
}

template <class T> void operator<< (valarray< valarray<T> >& mat, const T* ptr) {
  unsigned int n1 = mat.size();
  if (n1 == 0) return;
  unsigned int n2 = mat[0].size();
  for (unsigned int i = 0; i < n1; i++) {
    for (unsigned int j = 0; j < n2; j++) {
      mat[i][j] = ptr[i*n2+j];
    }
  }
}

template <class T> void operator<< (valarray< valarray< valarray<T> > >& tensor, const T* ptr) {
  unsigned int n1 = tensor.size();
  if (n1 == 0) return;
  unsigned int n2 = tensor[0].size();
  unsigned int n3 = tensor[0][0].size();
  for (unsigned int i = 0; i < n1; i++) {
    for (unsigned int j = 0; j < n2; j++) {
      for (unsigned int k = 0; k < n3; k++) {
        tensor[i][j][k] = ptr[i*n2*n3 + j*n3 + k];
      }
    }
  }
}

template <class T> void operator<< (valarray<T>& arr1, const valarray<T>& arr2) {
  unsigned int n1 = arr1.size(), n2 = arr2.size();
  if (n1 != n2) arr1.resize(n2);
  arr1 = arr2;
}

template <class T> void operator<< (valarray< valarray<T> >& mat1, 
                              const valarray< valarray<T> >& mat2) {
  unsigned int n1_1 = mat1.size(), n2_1 = mat2.size();
  unsigned int n1_2 = 0, n2_2 = 0;
  if (n1_1 > 0) n1_2 = mat1[0].size();
  if (n2_1 > 0) n2_2 = mat2[0].size();
  if (n1_1 != n2_1) mat1.resize(n2_1);
  if (n1_2 != n2_2) {for (unsigned int i = 0; i < n1_1; i++) mat1[i].resize(n2_2);}
  mat1 = mat2;
}

template <class T> void matrix_to_vec (const valarray< valarray<T> >& mat, T* vec) {
  unsigned int n1 = mat.size();
  if (n1 == 0) return;
  unsigned int n2 = mat[0].size();
  for (unsigned int i = 0; i < n1; i++) {
    for (unsigned int j = 0; j < n2; j++) {
      vec[i*n2+j] = mat[i][j];
    }
  }
}

template <class T> void tensor_to_vec (const valarray< valarray< valarray<T> > >& tensor, T* vec) {
  unsigned int n1 = tensor.size();
  if (n1 == 0) return;
  unsigned int n2 = tensor[0].size();
  unsigned int n3 = tensor[0][0].size();
  for (unsigned int i = 0; i < n1; i++) {
    for (unsigned int j = 0; j < n2; j++) {
      for (unsigned int k = 0; k < n3; k++) {
        vec[i*n2*n3 + j*n3 + k] = tensor[i][j][k];
      }
    }
  }
}

//---------------------------------------------------------------------------
// Instantiate instances for conversion from array to C++ structure.

template void operator<< (Bool_ARRAY&,  c_Bool*);
template void operator<< (Bool_MATRIX&, c_Bool*);

template void operator<< (Real_ARRAY&,  c_Real*);
template void operator<< (Real_MATRIX&, c_Real*);
template void operator<< (Real_TENSOR&, c_Real*);

template void operator<< (Complex_ARRAY&,  c_Complex*);
template void operator<< (Complex_MATRIX&, c_Complex*);
template void operator<< (Complex_TENSOR&, c_Complex*);

template void operator<< (Int_ARRAY&,  c_Int*);
template void operator<< (Int_MATRIX&, c_Int*);
template void operator<< (Int_TENSOR&, c_Int*);

//---------------------------------------------------------------------------
// Instantiate instances for transfer

template void operator<< (Real_ARRAY&,  const Real_ARRAY&);
template void operator<< (Real_MATRIX&, const Real_MATRIX&);
template void operator<< (Real_TENSOR&, const Real_TENSOR&);

template void operator<< (Complex_ARRAY&,  const Complex_ARRAY&);
template void operator<< (Complex_MATRIX&, const Complex_MATRIX&);
template void operator<< (Complex_TENSOR&, const Complex_TENSOR&);

template void operator<< (Int_ARRAY&,  const Int_ARRAY&);
template void operator<< (Int_MATRIX&, const Int_MATRIX&);
template void operator<< (Int_TENSOR&, const Int_TENSOR&);

//---------------------------------------------------------------------------

template void matrix_to_vec (const Bool_MATRIX&,     Bool*);
template void matrix_to_vec (const Complex_MATRIX&,  Complex*);
template void matrix_to_vec (const Real_MATRIX&,     Real*);
template void matrix_to_vec (const Int_MATRIX&,      Int*);

template void tensor_to_vec (const Complex_TENSOR&,  Complex*);
template void tensor_to_vec (const Real_TENSOR&,     Real*);
template void tensor_to_vec (const Int_TENSOR&,      Int*);

#define CONVERTER_TEMPLATES
#endif
