#include <iostream>
#include <stdlib.h>
#include "cpp_and_bmad.h"

using namespace std;

void ele_comp (const C_ele& x, const C_ele& y);

//

extern "C" void test_f_coord_(C_coord& c1, C_coord& c2);
extern "C" void test_f_twiss_(C_twiss& c1, C_twiss& c2);
extern "C" void test_f_xy_disp_(C_xy_disp& c1, C_xy_disp& c2);
extern "C" void test_f_floor_position_(C_floor_position& c1, C_floor_position& c2);
extern "C" void test_f_wig_term_(C_wig_term& c1, C_wig_term& c2);
extern "C" void test_f_taylor_term_(C_taylor_term& c1, C_taylor_term& c2);
extern "C" void test_f_taylor_(C_taylor& c1, C_taylor& c2);
extern "C" void test_f_sr_table_wake_(C_sr_table_wake& c1, C_sr_table_wake& c2);
extern "C" void test_f_sr_mode_wake_(C_sr_mode_wake& c1, C_sr_mode_wake& c2);
extern "C" void test_f_lr_wake_(C_lr_wake& c1, C_lr_wake& c2);
extern "C" void test_f_wake_(C_wake& c1, C_wake& c2);
extern "C" void test_f_control_(C_control& c1, C_control& c2);
extern "C" void test_f_param_(C_param& c1, C_param& c2);
extern "C" void test_f_amode_(C_amode& c1, C_amode& c2);
extern "C" void test_f_linac_mode_(C_linac_mode& c1, C_linac_mode& c2);
extern "C" void test_f_modes_(C_modes& c1, C_modes& c2);
extern "C" void test_f_bmad_com_(C_bmad_com& c1, C_bmad_com& c2);
extern "C" void test_f_em_field_(C_em_field& c1, C_em_field& c2);
extern "C" void test_f_ele_(C_ele& c1, C_ele& c2);
extern "C" void test_f_mode_info_(C_mode_info& c1, C_mode_info& c2);
extern "C" void test_f_lat_(C_lat& c1, C_lat& c2);

C_coord           c_coord_in(101, 102, 103, 104, 105, 106);
C_coord           c_coord_out(201, 202, 203, 204, 205, 206);
C_twiss           c_twiss_in(1, 2, 3, 4, 5, 6, 7, 8);
C_twiss           c_twiss_out(8, 7, 6, 5, 4, 3, 2, 1);
C_xy_disp         c_xy_disp_in(1, 2);
C_xy_disp         c_xy_disp_out(2, 1);
C_floor_position  c_floor_position_in(1, 2, 3, 4, 5, 6);
C_floor_position  c_floor_position_out(6, 5, 4, 3, 2, 1);
C_wig_term        c_wig_term_in(1, 2, 3, 4, 5, 6);
C_wig_term        c_wig_term_out(6, 5, 4, 3, 2, 1);
C_taylor_term     c_taylor_term_in(1, 2, 3, 4, 5, 6, 7);
C_taylor_term     c_taylor_term_out(7, 6, 5, 4, 3, 2, 1);
C_taylor          c_taylor_in(2, -1);
C_taylor          c_taylor_out(2, 1);
C_sr_table_wake   c_sr_table_wake_in(1, 2, 3);
C_sr_table_wake   c_sr_table_wake_out(3, 2, 1);
C_sr_mode_wake    c_sr_mode_wake_in(21, 22, 23, 24, 25, 26, 27, 28);
C_sr_mode_wake    c_sr_mode_wake_out(31, 32, 33, 34, 35, 36, 37, 38);
C_lr_wake         c_lr_wake_in(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1);
C_lr_wake         c_lr_wake_out(10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
C_wake            c_wake_in("ABCD", "XYZZY", 2, 2, 0, 1);
C_wake            c_wake_out("abcd", "xyzzy", 0, 0, 2, 2);
C_control         c_control_in(1, 2, 3, 4);
C_control         c_control_out(4, 3, 2, 1);
C_param           c_param_in;
C_param           c_param_out;
C_amode           c_amode_in(1, 2, 3, 4, 5, 6, 7);
C_amode           c_amode_out(11, 12, 13, 14, 15, 16, 17);
C_linac_mode      c_linac_mode_in(1, 2, 3, 4, 5, 6, 7);
C_linac_mode      c_linac_mode_out(11, 12, 13, 14, 15, 16, 17);
C_modes           c_modes_in(1, 2, 3, 4, 5, 6, 7,
                  c_amode_in, c_amode_out, c_amode_in, c_linac_mode_in);
C_modes           c_modes_out(11, 12, 13, 14, 15, 16, 17, 
                  c_amode_out, c_amode_out, c_amode_in, c_linac_mode_out);
C_bmad_com        c_bmad_com_in;
C_bmad_com        c_bmad_com_out;
C_em_field        c_em_field_in;
C_em_field        c_em_field_out;
C_ele             c_ele_in;
C_ele             c_ele_out;
C_mode_info       c_mode_info_in(1, 2, 3);
C_mode_info       c_mode_info_out(-1, -2, -3);
C_lat             c_lat_in;
C_lat             c_lat_out;

//---------------------------------------------------------------------------

void init_all_c_structs () {
 
  Real_Matrix m6a(M6_mat), m6b(M6_mat);
  Real_Matrix m3a(M3_mat), m3b(M3_mat), m3c(M3_mat);
  Real_Matrix m2a(M2_mat);
  Real_Array v6a(V6_array), v6b(V6_array), v3a(V3_array), v3b(V3_array), v3c(V3_array);
  Int_Array i6a(0, 6);

  for (int i = 0; i < 6; i++) {
    v6a[i] = i + 101;
    v6b[i] = i + 201;
    i6a[i] = i + 101;
    for (int j = 0; j < 6; j++) {
      m6a[i][j] = 10*i + j + 1; 
      m6b[i][j] = 10*i + j + 101;
    }
  }

  for (int i = 0; i < 3; i++) {
    v3a[i] = i + 101;
    v3b[i] = i + 201;
    v3c[i] = i + 301;
    for (int j = 0; j < 3; j++) {
      m3a[i][j] = 10*i + j + 1; 
      m3b[i][j] = 10*i + j + 101;
      m3c[i][j] = 10*i + j + 201;
    }
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      m2a[i][j] = 10*i + j + 1; 
    }
  }

  bool T = true, F = false;

  c_taylor_in.term[0] = C_taylor_term(1, 2, 3, 4, 5, 6, 7);
  c_taylor_in.term[1] = C_taylor_term(8, 9, 10, 11, 12, 13, 14); 

  c_taylor_out.term[0] = C_taylor_term(-1, -2, -3, -4, -5, -6, -7);
  c_taylor_out.term[1] = C_taylor_term(-8, -9, -10, -11, -12, -13, -14);

  c_wake_in.sr_table[0] = c_sr_table_wake_in;
  c_wake_in.sr_table[1] = c_sr_table_wake_out;
  c_wake_in.sr_mode_long[0] = c_sr_mode_wake_in;
  c_wake_in.sr_mode_long[1] = c_sr_mode_wake_out;
  c_wake_in.z_sr_mode_max = 100;
  c_wake_in.lr[0] = c_lr_wake_in;

  c_wake_out.z_sr_mode_max = 101;
  c_wake_out.sr_mode_trans[0] = c_sr_mode_wake_in;
  c_wake_out.sr_mode_trans[1] = c_sr_mode_wake_out;
  c_wake_out.lr[0] = C_lr_wake(-1, -2, -3, -4, 5, -6, -7, -8, -9, 10, 1);
  c_wake_out.lr[1] = C_lr_wake(-11, -12, -13, -14, -15, -16, -17, -18, -19, 20, 0);

  c_param_in  = C_param(1, 2, 3, m6a, m6b, 11, 12, 13, 14, 15, 1, 0, 1);
  c_param_out = C_param(11, 12, 13, m6b, m6a, 111, 112, 113, 114, 115, 0, 0, 0);

  c_bmad_com_in  = C_bmad_com(2, v6a, 3, 4, 5, 6, 7, 8, 
                                  1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1);
  c_bmad_com_out = C_bmad_com(12, v6b, 13, 14, 15, 16, 17, 18, 
                                  1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0);

  c_em_field_in  = C_em_field(v3a, v3b, v3c, m3a, m3b, m3c, 77);
  c_em_field_out = C_em_field(v3c, v3b, v3a, m3c, m3b, m3a, -77);

  c_ele_in.name = "1234";
  c_ele_in.type = "abcd";
  c_ele_in.alias = "5678";
  c_ele_in.attribute_name = "efgh";
  c_ele_in.value[1] = -1;
  c_ele_in.value[Bmad::N_ATTRIB_MAXX] = -Bmad::N_ATTRIB_MAXX;
  c_ele_in.gen0 = v6a;
  c_ele_in.vec0 = v6b;
  c_ele_in.mat6 = m6a;
  c_ele_in.c_mat = m2a;
  c_ele_in.gamma_c = 51;
  c_ele_in.s = 52;
  c_ele_in.key = 54;
  c_ele_in.sub_key = 55;
  c_ele_in.control_type = 56;
  c_ele_in.ix_value = 57;
  c_ele_in.n_slave = 58;
  c_ele_in.ix1_slave = 59;
  c_ele_in.ix2_slave = 60;
  c_ele_in.n_lord = 61;
  c_ele_in.ic1_lord = 62;
  c_ele_in.ic2_lord = 63;
  c_ele_in.ix_pointer = 64;
  c_ele_in.ixx = 65;
  c_ele_in.ix_ele = 66;
  c_ele_in.mat6_calc_method = 67;
  c_ele_in.tracking_method = 68;
  c_ele_in.field_calc = 69;
  c_ele_in.num_steps = 70;
  c_ele_in.integrator_order = 71;
  c_ele_in.ptc_kind = 72;
  c_ele_in.taylor_order = 73;
  c_ele_in.aperture_at = 74;
  c_ele_in.coupler_at = 75;
  c_ele_in.symplectify = T;
  c_ele_in.mode_flip = F;
  c_ele_in.multipoles_on = T;
  c_ele_in.map_with_offsets = F;
  c_ele_in.field_master = T;
  c_ele_in.is_on = F;
  c_ele_in.internal_logic = T;
  c_ele_in.logic = F;
  c_ele_in.on_an_i_beam = T;
  c_ele_in.csr_calc_on = F;
  c_ele_in.floor = c_floor_position_in;
  c_ele_in.x = c_xy_disp_in;
  c_ele_in.y = c_xy_disp_out;
  c_ele_in.a = c_twiss_in;
  c_ele_in.b = c_twiss_out;
  c_ele_in.z = c_twiss_in;

  c_ele_out = c_ele_in;

  c_ele_in.descrip = "descrip";
  c_ele_in.wake = c_wake_in;
  c_ele_in.const_arr.resize(6);
  c_ele_in.const_arr = -v6b;
  c_ele_in.taylor[0].term.resize(2);
  c_ele_in.taylor[1].term.resize(2);
  c_ele_in.taylor[0] = c_taylor_in;
  c_ele_in.taylor[1] = c_taylor_out;
  c_ele_in.r.resize(2);
  for (int i = 0; i < 2; i++) {
    c_ele_in.r[i].resize(3);
    for (int j = 0; j < 3; j++) {
      c_ele_in.r[i][j] = m3c[i][j];
    }
  }
  c_ele_in.wig_term.resize(2);
  c_ele_in.wig_term[0] = c_wig_term_in;
  c_ele_in.wig_term[1] = c_wig_term_out;
  c_ele_in.a_pole.resize(Bmad::N_POLE_MAXX+1);
  c_ele_in.b_pole.resize(Bmad::N_POLE_MAXX+1);
  for (int i = 0; i < Bmad::N_POLE_MAXX + 1; i++) {
    c_ele_in.a_pole[i] = -i - 100;
    c_ele_in.b_pole[i] = -i - 200;
  }



  c_lat_in.name = "abc";
  c_lat_in.lattice = "def";
  c_lat_in.input_file_name = "123";
  c_lat_in.title = "title";
  c_lat_in.version = 21;
  c_lat_in.n_ele_use = 22;
  c_lat_in.n_ele_max = 1;
  c_lat_in.n_control_max = 2;
  c_lat_in.n_ic_max = 6;
  c_lat_in.input_taylor_order = 34;

  c_lat_in.ele.resize(2);
  c_lat_in.ele[0] = c_ele_in;
  c_lat_in.ele[0].ix_ele = 0;
  c_lat_in.ele[1] = c_ele_out;
  c_lat_in.ele[1].ix_ele = 1;
  c_lat_in.control.resize(2);
  c_lat_in.control[0] = c_control_in;
  c_lat_in.control[1] = c_control_out;
  c_lat_in.ic.resize(6);
  c_lat_in.ic = i6a;

  c_lat_in.ele_init = c_ele_in;
  c_lat_in.ele_init.name = "ele_init";

  c_lat_out = c_lat_in;

}

//---------------------------------------------------------------------------

extern "C" void test_c_coord_(coord_struct* f1, coord_struct* f2, int& c_ok) {

  C_coord c1, c2;

  init_all_c_structs();

  coord_struct* f3;
  init_coord_struct_(f3);

  f1 >> c1;

  if (c1 == c_coord_in) {
    cout << " C_side_convert: Coord F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: COORD F to C: FAILED!!" << endl;
    c_ok = 0;
  }

  test_f_coord_(c1, c2);

  if (c2 == c_coord_out) {
    cout << " F_side_convert: Coord F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: COORD F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_twiss_(twiss_struct* f1, twiss_struct* f2, int& c_ok) {

  C_twiss c1, c2;

  f1 >> c1; 

  if (c1 == c_twiss_in) {
    cout << " C_side_convert: Twiss F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: TWISS F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_twiss_(c1, c2);

  if (c2 == c_twiss_out) {
    cout << " F_side_convert: Twiss F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: TWISS F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_xy_disp_(xy_disp_struct* f1, xy_disp_struct* f2, int& c_ok) {

  C_xy_disp c1, c2;

  f1 >> c1; 

  if (c1 == c_xy_disp_in) {
    cout << " C_side_convert: Xy_Disp F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: XY_DISP F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_xy_disp_(c1, c2);

  if (c2 == c_xy_disp_out) {
    cout << " F_side_convert: Xy_Disp F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: XY_DISP F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_floor_position_(floor_position_struct* f1, floor_position_struct* f2, int& c_ok) {

  C_floor_position c1, c2;

  f1 >> c1; 

  if (c1 == c_floor_position_in) {
    cout << " C_side_convert: floor_position F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: FLOOR_POSITION F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_floor_position_(c1, c2);

  if (c2 == c_floor_position_out) {
    cout << " F_side_convert: floor_position F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: FLOOR_POSITION F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_wig_term_(wig_term_struct* f1, wig_term_struct* f2, int& c_ok) {

  C_wig_term c1, c2;

  f1 >> c1; 

  if (c1 == c_wig_term_in) {
    cout << " C_side_convert: wig_term F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: WIG_TERM F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_wig_term_(c1, c2);

  if (c2 == c_wig_term_out) {
    cout << " F_side_convert: wig_term F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: WIG_TERM F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_taylor_term_(taylor_term_struct* f1, taylor_term_struct* f2, int& c_ok) {

  C_taylor_term c1, c2;

  f1 >> c1; 

  if (c1 == c_taylor_term_in) {
    cout << " C_side_convert: taylor_term F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: TAYLOR_TERM F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_taylor_term_(c1, c2);

  if (c2 == c_taylor_term_out) {
    cout << " F_side_convert: taylor_term F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: TAYLOR_TERM F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_taylor_(taylor_struct* f1, taylor_struct* f2, int& c_ok) {

  C_taylor c1, c2;
  
  f1 >> c1; 

  if (c1 == c_taylor_in) {
    cout << " C_side_convert: taylor F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: TAYLOR F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_taylor_(c1, c2);

  if (c2 == c_taylor_out) {
    cout << " F_side_convert: taylor F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: TAYLOR F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_sr_table_wake_(sr_table_wake_struct* f1, sr_table_wake_struct* f2, int& c_ok) {

  C_sr_table_wake c1, c2;

  f1 >> c1; 

  if (c1 == c_sr_table_wake_in) {
    cout << " C_side_convert: sr_table_wake F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: sr_table_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_sr_table_wake_(c1, c2);

  if (c2 == c_sr_table_wake_out) {
    cout << " F_side_convert: sr_table_wake F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: sr_table_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_sr_mode_wake_(sr_mode_wake_struct* f1, sr_mode_wake_struct* f2, int& c_ok) {

  C_sr_mode_wake c1, c2;

  f1 >> c1; 

  if (c1 == c_sr_mode_wake_in) {
    cout << " C_side_convert: sr_mode_wake F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: sr_mode_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_sr_mode_wake_(c1, c2);

  if (c2 == c_sr_mode_wake_out) {
    cout << " F_side_convert: sr_mode_wake F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: sr_mode_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_lr_wake_(lr_wake_struct* f1, lr_wake_struct* f2, int& c_ok) {

  C_lr_wake c1, c2;

  f1 >> c1; 

  if (c1 == c_lr_wake_in) {
    cout << " C_side_convert: lr_wake F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: LR_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_lr_wake_(c1, c2);

  if (c2 == c_lr_wake_out) {
    cout << " F_side_convert: lr_wake F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: LR_WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_wake_(wake_struct* f1, wake_struct* f2, int& c_ok) {

  C_wake c1, c2;
  
  f1 >> c1; 

  if (c1 == c_wake_in) {
    cout << " C_side_convert: wake F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_wake_(c1, c2);

  if (c2 == c_wake_out) {
    cout << " F_side_convert: wake F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: WAKE F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_control_(control_struct* f1, control_struct* f2, int& c_ok) {

  C_control c1, c2;

  f1 >> c1; 

  if (c1 == c_control_in) {
    cout << " C_side_convert: control F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: control F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_control_(c1, c2);

  if (c2 == c_control_out) {
    cout << " F_side_convert: control F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: control F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_param_(lat_param_struct* f1, lat_param_struct* f2, int& c_ok) {

  C_param c1, c2;

  f1 >> c1; 

  if (c1 == c_param_in) {
    cout << " C_side_convert: param F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: param F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_param_(c1, c2);

  if (c2 == c_param_out) {
    cout << " F_side_convert: param F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: param F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_amode_(anormal_mode_struct* f1, anormal_mode_struct* f2, int& c_ok) {

  C_amode c1, c2;

  f1 >> c1; 

  if (c1 == c_amode_in) {
    cout << " C_side_convert: amode F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: amode F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_amode_(c1, c2);

  if (c2 == c_amode_out) {
    cout << " F_side_convert: amode F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: amode F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_linac_mode_(linac_normal_mode_struct* f1, linac_normal_mode_struct* f2, int& c_ok) {

  C_linac_mode c1, c2;

  f1 >> c1; 

  if (c1 == c_linac_mode_in) {
    cout << " C_side_convert: linac_mode F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: linac_mode F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_linac_mode_(c1, c2);

  if (c2 == c_linac_mode_out) {
    cout << " F_side_convert: linac_mode F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: linac_mode F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_modes_(normal_modes_struct* f1, normal_modes_struct* f2, int& c_ok) {

  C_modes c1, c2;

  f1 >> c1; 

  if (c1 == c_modes_in) {
    cout << " C_side_convert: modes F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: modes F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_modes_(c1, c2);

  if (c2 == c_modes_out) {
    cout << " F_side_convert: modes F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: modes F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_bmad_com_(int& c_ok) {

  C_bmad_com c1, c2;

  bmad_com_to_c_(c1); 

  if (c1 == c_bmad_com_in) {
    cout << " C_side_convert: bmad_com F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: bmad_com F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_bmad_com_(c1, c2);

  if (c2 == c_bmad_com_out) {
    cout << " F_side_convert: bmad_com F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: bmad_com F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  bmad_com_to_f_(c2);

}

//---------------------------------------------------------------------------

extern "C" void test_c_em_field_(em_field_struct* f1, em_field_struct* f2, int& c_ok) {

  C_em_field c1, c2;

  f1 >> c1; 

  if (c1 == c_em_field_in) {
    cout << " C_side_convert: em_field F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: em_field F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_em_field_(c1, c2);

  if (c2 == c_em_field_out) {
    cout << " F_side_convert: em_field F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: em_field F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_ele_(ele_struct* f1, ele_struct* f2, int& c_ok) {

  C_ele c1, c2;

  f1 >> c1; 

  if (c1 == c_ele_in) {
    cout << " C_side_convert: ele F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: ele F TO C: FAILED!!" << endl;
    ele_comp (c1, c_ele_in);
    c_ok = 0;
  }


  test_f_ele_(c1, c2);

  if (c2 == c_ele_out) {
    cout << " F_side_convert: ele F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: ele F TO C: FAILED!!" << endl;
    ele_comp (c2, c_ele_out);
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_mode_info_(mode_info_struct* f1, mode_info_struct* f2, int& c_ok) {

  C_mode_info c1, c2;

  f1 >> c1; 

  if (c1 == c_mode_info_in) {
    cout << " C_side_convert: mode_info F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: mode_info F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_mode_info_(c1, c2);

  if (c2 == c_mode_info_out) {
    cout << " F_side_convert: mode_info F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: mode_info F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}

//---------------------------------------------------------------------------

extern "C" void test_c_lat_(lat_struct* f1, lat_struct* f2, int& c_ok) {

  C_lat c1, c2;

  f1 >> c1; 

  if (c1 == c_lat_in) {
    cout << " C_side_convert: lat F to C: OK" << endl;
  } else {
    cout << " C_SIDE_CONVERT: lat F TO C: FAILED!!" << endl;
    c_ok = 0;
  }


  test_f_lat_(c1, c2);

  if (c2 == c_lat_out) {
    cout << " F_side_convert: lat F to C: OK" << endl;
  } else {
    cout << " F_SIDE_CONVERT: lat F TO C: FAILED!!" << endl;
    c_ok = 0;
  }

  c2 >> f2;

}


//---------------------------------------------------------------------------
