#include <stdio.h>

typedef const double    Re;
typedef const int       Int;
typedef const char*     Char;
typedef const bool      Bool;
typedef const double*   ReArr;
typedef const int*      IntArr;

struct zzz_struct {};

struct C_zzz {
  int i, j;
};

//--------------------------------------------------------------------
//--------------------------------------------------------------------

extern "C" void zzz_to_c (zzz_struct*, C_zzz&);
extern "C" void zzz_to_f2 (zzz_struct*, IntArr);

extern "C" void zzz_to_f (C_zzz& c_zzz, zzz_struct* f_zzz) {
  int c_int[2];
  c_int[0] = c_zzz.i;
  c_int[1] = c_zzz.j;
  zzz_to_f2 (f_zzz, c_int);
} 

extern "C" void zzz_to_c2 (C_zzz& c_zzz, IntArr c_vec) {
  c_zzz.i = c_vec[0];
  c_zzz.j = c_vec[1];
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

extern "C" void test2_f_zzz (C_zzz&);

//--------------------------------------------------------------------

extern "C" void test_c_zzz (zzz_struct* f_zzz, int& c_ok) {

  C_zzz c_zzz;

  zzz_to_c (f_zzz, c_zzz);

  printf ("C_side_convert F->C: %d  %d\n", c_zzz.i, c_zzz.j);

  //---------------

  c_zzz.i = 44;
  c_zzz.j = 66;
  test2_f_zzz (c_zzz);

  printf ("F_side_convert F->C: %d  %d\n", c_zzz.i, c_zzz.j);

  //---------------

  c_zzz.i = 101;
  c_zzz.j = 666;
  zzz_to_f (c_zzz, f_zzz);

}

