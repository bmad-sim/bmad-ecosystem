#ifndef BMAD_STD_TYPEDEF

using namespace std;

typedef bool               Bool;
typedef complex<double>    Complex;
typedef double             Real;
typedef int                Int;
typedef long int           Int8;
typedef char*              Char;

typedef const bool               c_Bool;
typedef const Complex            c_Complex;
typedef const double             c_Real;
typedef const int                c_Int;
typedef const long int           c_Int8;
typedef const string             c_String;
typedef const char*              c_Char;

typedef const bool*              c_BoolArr;
typedef const Complex*           c_ComplexArr;
typedef const double*            c_RealArr;
typedef const int*               c_IntArr;
typedef const long int*          c_Int8Arr;

typedef valarray<bool>           Bool_ARRAY;
typedef valarray<Complex>        Complex_ARRAY;
typedef valarray<double>         Real_ARRAY;
typedef valarray<int>            Int_ARRAY;
typedef valarray<string>         String_ARRAY;

typedef valarray<Bool_ARRAY>     Bool_MATRIX;
typedef valarray<Complex_ARRAY>  Complex_MATRIX;
typedef valarray<Real_ARRAY>     Real_MATRIX;
typedef valarray<Int_ARRAY>      Int_MATRIX;

typedef valarray<Bool_MATRIX>      Bool_TENSOR;
typedef valarray<Complex_MATRIX>   Complex_TENSOR;
typedef valarray<Real_MATRIX>      Real_TENSOR;
typedef valarray<Int_MATRIX>       Int_TENSOR;

#define BMAD_STD_TYPEDEF
#endif
