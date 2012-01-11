#include "cpp_and_bmad.h"

//---------------------------------------------------------------------------

template <class T> void operator<< (valarray<T>& arr, const T* ptr) {
  int n = arr.size();
  for (int i = 0; i < n; i++) arr[i] = ptr[i];
}

template <class T> void operator<< (valarray< valarray<T> >& mat, const T* ptr) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      mat[i][j] = ptr[i*n2+j];
    }
  }
}

template <class T> void operator<< (valarray<T>& arr1, const valarray<T>& arr2) {
  int n1 = arr1.size(), n2 = arr2.size();
  if (n1 != n2) arr1.resize(n2);
  arr1 = arr2;
}

template <class T> void operator<< (valarray< valarray<T> >& mat1, 
                              const valarray< valarray<T> >& mat2) {
  int n1_1 = mat1.size(), n2_1 = mat2.size();
  int n1_2 = 0, n2_2 = 0;
  if (n1_1 > 0) n1_2 = mat1[0].size();
  if (n2_1 > 0) n2_2 = mat2[0].size();
  if (n1_1 != n2_1) mat1.resize(n2_1);
  if (n1_2 != n2_2) {for (int i = 0; i < n1_1; i++) mat1[i].resize(n2_2);}
  mat1 = mat2;
}

//---------------------------------------------------------------------------
 template void operator<< (Real_Array&, const double*);
 template void operator<< (Real_Matrix&, const double*);
 template void operator<< (Int_Array&, const int*);

 template void operator<< (Real_Array&, const Real_Array&);
 template void operator<< (Real_Matrix&, const Real_Matrix&);
 template void operator<< (Int_Array&, const Int_Array&);

 template void operator<< (C_taylor_array&, const C_taylor_array&);
 template void operator<< (C_wig_term_array&, const C_wig_term_array&);

//---------------------------------------------------------------------------

void matrix_to_array (const Real_Matrix& mat, double* arr) {
  int n1 = mat.size();
  if (n1 == 0) return;
  int n2 = mat[0].size();
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      arr[i*n2+j] = mat[i][j];
    }
  }
}

//---------------------------------------------------------------------------
// Coord

extern "C" void coord_to_f2_(coord_struct*, ReArr, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void coord_to_f_(C_coord& c, coord_struct* f) {
  coord_to_f2_(f, &c.vec[0], c.s, c.t, c.spin1.real(), c.spin1.imag(), c.spin2.real(), c.spin2.imag(), 
                c.e_field_x, c.e_field_y, c.phase_x, c.phase_y);
}

extern "C" void coord_to_c2_(C_coord& c, ReArr vec, Re& s, Re& t, Re& sp1_re, Re& sp1_im, Re& sp2_re, Re& sp2_im, 
                             Re& field_x, Re& field_y, Re& p_x, Re& p_y) {
  c = C_coord(vec, s, t, Complx(sp1_re, sp1_im), Complx(sp2_re, sp2_im), field_x, field_y, p_x, p_y);
}

void operator>> (C_coord& c, coord_struct* f) {
  coord_to_f_(c, f);
}

void operator>> (coord_struct* f, C_coord& c) {
  coord_to_c_(f, c);
}

//---------------------------------------------------------------------------
// Twiss

extern "C" void twiss_to_f2_(twiss_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void twiss_to_f_(C_twiss& c, twiss_struct* f) {
  twiss_to_f2_(f, c.beta, c.alpha, c.gamma, c.phi, c.eta, c.etap, 
                  c.sigma, c.sigma_p, c.emit, c.norm_emit);
}

extern "C" void twiss_to_c2_(C_twiss& c, Re& beta, Re& alpha, Re& gamma, Re& phi, 
                          Re& eta, Re& etap, Re& sigma, Re& sigma_p, Re& emit, Re& norm_emit) {
  c = C_twiss(beta, alpha, gamma, phi, eta, etap, sigma, sigma_p, emit, norm_emit);
}

void operator>> (C_twiss& c, twiss_struct* f) {
  twiss_to_f_(c, f);
}

void operator>> (twiss_struct* f, C_twiss& c) {
  twiss_to_c_(f, c);
}


//---------------------------------------------------------------------------
// XY_Disp

extern "C" void xy_disp_to_f2_(xy_disp_struct*, Re&, Re&);

extern "C" void xy_disp_to_f_(C_xy_disp& c, xy_disp_struct* f) {
  xy_disp_to_f2_(f, c.eta, c.etap);
}

extern "C" void xy_disp_to_c2_(C_xy_disp& c, Re& eta, Re& etap) {
  c = C_xy_disp(eta, etap);
}

void operator>> (C_xy_disp& c, xy_disp_struct* f) {
  xy_disp_to_f_(c, f);
}

void operator>> (xy_disp_struct* f, C_xy_disp& c) {
  xy_disp_to_c_(f, c);
}


//-----------------------------------------------------------------------------
// Floor_position

extern "C" void floor_position_to_f2_(floor_position_struct*, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void floor_position_to_f_(C_floor_position& c, floor_position_struct* f) {
  floor_position_to_f2_(f, c.x, c.y, c.z, c.theta, c.phi, c.psi);
}

extern "C" void floor_position_to_c2_(C_floor_position& c, Re& x, Re& y, Re& z, 
                                                      Re& theta, Re& phi, Re& psi) {
 c = C_floor_position(x, y, z, theta, phi, psi);
}

void operator>> (C_floor_position& c, floor_position_struct* f) {
  floor_position_to_f_(c, f);
}

void operator>> (floor_position_struct* f, C_floor_position& c) {
  floor_position_to_c_(f, c);
}

//-----------------------------------------------------------------------------
// Wig_term

extern "C" void wig_term_to_f2_(wig_term_struct*, Re&, Re&, Re&, Re&, Re&, Int&);

extern "C" void wig_term_to_f_(C_wig_term& c, wig_term_struct* f) {
  wig_term_to_f2_(f, c.coef, c.kx, c.ky, c.kz, c.phi_z, c.type);
}

extern "C" void wig_term_to_c2_(C_wig_term& c, Re& coef, Re& kx, Re& ky, 
                                                      Re& kz, Re& phi_z, Int& type) {
 c = C_wig_term(coef, kx, ky, kz, phi_z, type);
}

void operator>> (C_wig_term& c, wig_term_struct* f) {
  wig_term_to_f_(c, f);
}

void operator>> (wig_term_struct* f, C_wig_term& c) {
  wig_term_to_c_(f, c);
}

//---------------------------------------------------------------------------
// Taylor_term

extern "C" void taylor_term_to_f2_(taylor_term_struct*, Re&, IntArr);

extern "C" void taylor_term_to_f_(C_taylor_term& c, taylor_term_struct* f) {
  taylor_term_to_f2_(f, c.coef, &c.expn[0]);
}

extern "C" void taylor_term_to_c2_(C_taylor_term& c, Re& coef, int expn[]) {
  c = C_taylor_term(coef, expn);
}

void operator>> (C_taylor_term& c, taylor_term_struct* f) {
  taylor_term_to_f_(c, f);
}

void operator>> (taylor_term_struct* f, C_taylor_term& c) {
  taylor_term_to_c_(f, c);
}

//---------------------------------------------------------------------------
// Taylor

extern "C" void taylor_to_f2_(taylor_struct*, Int&, Re&);
extern "C" void taylor_term_in_taylor_to_f2_(taylor_struct*, Int&, Re&, IntArr);

extern "C" void taylor_to_f_(C_taylor& c, taylor_struct* f) {
  int n_term = c.term.size();
  taylor_to_f2_(f, n_term, c.ref);
  for (int i = 0; i < n_term; i++) {
    taylor_term_in_taylor_to_f2_(f, i+1, c.term[i].coef, &c.term[i].expn[0]);
  }
}

extern "C" void taylor_to_c2_(C_taylor& c, Int& n_term, Re& ref) {
  if (c.term.size() != n_term) c.term.resize(n_term);
  c.ref = ref;
}

extern "C" void taylor_term_in_taylor_to_c2_(C_taylor& c, Int& it, Re& coef, int expn[]) {
  c.term[it-1] = C_taylor_term(coef, expn);
}


void operator>> (C_taylor& c, taylor_struct* f) {
  taylor_to_f_(c, f);
}

void operator>> (taylor_struct* f, C_taylor& c) {
  taylor_to_c_(f, c);
}

C_taylor& C_taylor::operator= (const C_taylor& c) {
  if (term.size() != c.term.size()) term.resize(c.term.size());
  ref = c.ref;
  term = c.term;
  return *this;
}

//---------------------------------------------------------------------------
// rf_wake_sr_table

extern "C" void rf_wake_sr_table_to_f2_(rf_wake_sr_table_struct*, Re&, Re&, Re&);

extern "C" void rf_wake_sr_table_to_f_(C_rf_wake_sr_table& c, rf_wake_sr_table_struct* f) {
  rf_wake_sr_table_to_f2_(f, c.z, c.longitudinal, c.transverse);
}

extern "C" void rf_wake_sr_table_to_c2_(C_rf_wake_sr_table& c, Re& z, Re& lw, Re& tw) {
  c = C_rf_wake_sr_table(z, lw, tw);
}

void operator>> (C_rf_wake_sr_table& c, rf_wake_sr_table_struct* f) {
  rf_wake_sr_table_to_f_(c, f);
}

void operator>> (rf_wake_sr_table_struct* f, C_rf_wake_sr_table& c) {
  rf_wake_sr_table_to_c_(f, c);
}

//---------------------------------------------------------------------------
// rf_wake_sr_mode

extern "C" void rf_wake_sr_mode_to_f2_(rf_wake_sr_mode_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void rf_wake_sr_mode_to_f_(C_rf_wake_sr_mode& c, rf_wake_sr_mode_struct* f) {
  rf_wake_sr_mode_to_f2_(f, c.amp, c.damp, c.k, c.phi,
                          c.b_sin, c.b_cos, c.a_sin, c.a_cos);
}

extern "C" void rf_wake_sr_mode_to_c2_(C_rf_wake_sr_mode& c, Re& amp, Re& damp, Re& k, Re& phi,
                          Re& b_sin, Re& b_cos, Re& a_sin, Re& a_cos) {
  c = C_rf_wake_sr_mode(amp, damp, k, phi, b_sin, b_cos, a_sin, a_cos);
}

void operator>> (C_rf_wake_sr_mode& c, rf_wake_sr_mode_struct* f) {
  rf_wake_sr_mode_to_f_(c, f);
}

void operator>> (rf_wake_sr_mode_struct* f, C_rf_wake_sr_mode& c) {
  rf_wake_sr_mode_to_c_(f, c);
}

//---------------------------------------------------------------------------
// rf_wake_lr

extern "C" void rf_wake_lr_to_f2_(rf_wake_lr_struct*, Re&, Re&, Re&, Re&, Re&,
                                      Re&, Re&, Re&, Re&, Re&, Int&, Int&);

extern "C" void rf_wake_lr_to_f_(C_rf_wake_lr& c, rf_wake_lr_struct* f) {
  rf_wake_lr_to_f2_(f, c.freq, c.freq_in, c.R_over_Q, c.Q, c.angle, 
          c.b_sin, c.b_cos, c.a_sin, c.a_cos, c.t_ref, c.m, c.polarized);
}

extern "C" void rf_wake_lr_to_c2_(C_rf_wake_lr& c, Re& freq, Re& freq_in, 
                  Re& R_over_Q, Re& Q, Re& ang, Re& n_sin, Re& n_cos, 
                  Re& s_sin, Re& s_cos, Re& t_ref, Int& m, Int& pol) {
  c = C_rf_wake_lr(freq, freq_in, R_over_Q, Q, ang, 
                    n_sin, n_cos, s_sin, s_cos, t_ref, m, pol);
}

void operator>> (C_rf_wake_lr& c, rf_wake_lr_struct* f) {
  rf_wake_lr_to_f_(c, f);
}

void operator>> (rf_wake_lr_struct* f, C_rf_wake_lr& c) {
  rf_wake_lr_to_c_(f, c);
}

//---------------------------------------------------------------------------
// rf_wake

extern "C" void rf_wake_to_f2_(rf_wake_struct*, Char, Int&, Char, Int&, Re&, Int&, Int&, Int&, Int&);
extern "C" void rf_wake_sr_table_in_rf_wake_to_f2_(rf_wake_struct*, Int&, Re&, Re&, Re&);
extern "C" void rf_wake_sr_mode_long_in_rf_wake_to_f2_(rf_wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);
extern "C" void rf_wake_sr_mode_trans_in_rf_wake_to_f2_(rf_wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);
extern "C" void rf_wake_lr_in_rf_wake_to_f2_(rf_wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&,
                                                         Re&, Re&, Re&, Re&, Re&, Int&, Int&);

extern "C" void rf_wake_to_f_(C_rf_wake& c, rf_wake_struct* f) {
  int n_lr = c.lr.size();
  int n_sr_table = c.sr_table.size(); 
  int n_sr_mode_long = c.sr_mode_long.size(); 
  int n_sr_mode_trans = c.sr_mode_trans.size(); 
  const char* srf = c.sr_file.data();     int n_srf = c.sr_file.length();
  const char* lrf = c.lr_file.data();     int n_lrf = c.lr_file.length();
  rf_wake_to_f2_(f, srf, n_srf, lrf, n_lrf, c.z_sr_mode_max, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr);
  for (int i = 0; i < n_sr_table; i++) {
    rf_wake_sr_table_in_rf_wake_to_f2_(f, i, c.sr_table[i].z, c.sr_table[i].longitudinal, c.sr_table[i].transverse);
  }
  for (int i = 0; i < n_sr_mode_long; i++) {
    rf_wake_sr_mode_long_in_rf_wake_to_f2_(f, i+1, c.sr_mode_long[i].amp, c.sr_mode_long[i].damp, 
        c.sr_mode_long[i].k, c.sr_mode_long[i].phi, c.sr_mode_long[i].b_sin, 
        c.sr_mode_long[i].b_cos, c.sr_mode_long[i].a_sin, c.sr_mode_long[i].a_cos);
  }
  for (int i = 0; i < n_sr_mode_trans; i++) {
    rf_wake_sr_mode_trans_in_rf_wake_to_f2_(f, i+1, c.sr_mode_trans[i].amp, c.sr_mode_trans[i].damp, 
        c.sr_mode_trans[i].k, c.sr_mode_trans[i].phi, c.sr_mode_trans[i].b_sin, 
        c.sr_mode_trans[i].b_cos, c.sr_mode_trans[i].a_sin, c.sr_mode_trans[i].a_cos);
  }
  for (int i = 0; i < n_lr; i++) {
    rf_wake_lr_in_rf_wake_to_f2_(f, i+1, c.lr[i].freq, c.lr[i].freq_in, 
        c.lr[i].R_over_Q, c.lr[i].Q, c.lr[i].angle, c.lr[i].b_sin, 
        c.lr[i].b_cos, c.lr[i].a_sin, c.lr[i].a_cos, c.lr[i].t_ref, c.lr[i].m, c.lr[i].polarized);
  }
}

extern "C" void rf_wake_to_c2_(C_rf_wake& c, char* srf, char* lrf, Re& z_cut, Int& n_sr_table, Int& n_sr_mode_long,
                            Int& n_sr_mode_trans, Int& n_lr) {
  if (c.sr_table.size() != n_sr_table) c.sr_table.resize(n_sr_table);
  if (c.sr_mode_long.size() != n_sr_mode_long) c.sr_mode_long.resize(n_sr_mode_long);
  if (c.sr_mode_trans.size() != n_sr_mode_trans) c.sr_mode_trans.resize(n_sr_mode_trans);
  c.sr_file = srf;
  if (c.lr.size() != n_lr) c.lr.resize(n_lr);
  c.lr_file = lrf;
  c.z_sr_mode_max = z_cut;
}

extern "C" void rf_wake_sr_table_in_rf_wake_to_c2_(C_rf_wake& c, Int& it, Re& z, Re& l, Re& t) {
  c.sr_table[it] = C_rf_wake_sr_table(z, l, t);
}

extern "C" void rf_wake_sr_mode_long_in_rf_wake_to_c2_(C_rf_wake& c, Int& it, Re& a, Re& d, 
                      Re& f, Re& p, Re& ns, Re& nc, Re& ss, Re& sc) {
  c.sr_mode_long[it-1] = C_rf_wake_sr_mode(a, d, f, p, ns, nc, ss, sc);
}

extern "C" void rf_wake_sr_mode_trans_in_rf_wake_to_c2_(C_rf_wake& c, Int& it, Re& a, Re& d, 
                      Re& f, Re& p, Re& ns, Re& nc, Re& ss, Re& sc) {
  c.sr_mode_trans[it-1] = C_rf_wake_sr_mode(a, d, f, p, ns, nc, ss, sc);
}

extern "C" void rf_wake_lr_in_rf_wake_to_c2_(C_rf_wake& c, Int& it, Re& f, Re& k, 
           Re& i, Re& q, Re& ang, Re& ns, Re& nc, Re& ss, Re& sc, Re& t_ref, Int& m, Int& pol) {
  c.lr[it-1] = C_rf_wake_lr(f, k, i, q, ang, ns, nc, ss, sc, t_ref, m, pol);
}


void operator>> (C_rf_wake& c, rf_wake_struct* f) {
  rf_wake_to_f_(c, f);
}

void operator>> (rf_wake_struct* f, C_rf_wake& c) {
  rf_wake_to_c_(f, c);
}

C_rf_wake& C_rf_wake::operator= (const C_rf_wake& c) {
  if (sr_table.size() != c.sr_table.size()) sr_table.resize(c.sr_table.size());
  sr_table = c.sr_table;
  if (sr_mode_long.size() != c.sr_mode_long.size()) sr_mode_long.resize(c.sr_mode_long.size());
  sr_mode_long = c.sr_mode_long;
  if (sr_mode_trans.size() != c.sr_mode_trans.size()) sr_mode_trans.resize(c.sr_mode_trans.size());
  sr_mode_trans = c.sr_mode_trans;
  sr_file = c.sr_file;
  if (lr.size() != c.lr.size()) lr.resize(c.lr.size());
  lr = c.lr;
  lr_file = c.lr_file;
  return *this;
}

//---------------------------------------------------------------------------
// control

extern "C" void control_to_f2_(control_struct*, Re&, Int&, Int&, Int&, Int&);

extern "C" void control_to_f_(C_control& c, control_struct* f) {
  control_to_f2_(f, c.coef, c.ix_lord, c.ix_slave, c.ix_branch, c.ix_attrib);
}

extern "C" void control_to_c2_(C_control& c, Re& coef, Int& il, Int& is, Int& pl, Int& ia) {
  c = C_control(coef, il, is, pl, ia);
}

void operator>> (C_control& c, control_struct* f) {
  control_to_f_(c, f);
}

void operator>> (control_struct* f, C_control& c) {
  control_to_c_(f, c);
}

//---------------------------------------------------------------------------
// lat_param

extern "C" void lat_param_to_f2_(lat_param_struct*, Re&, Re&, Re&, ReArr, ReArr, 
                               Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&);

extern "C" void lat_param_to_f_(C_lat_param& c, lat_param_struct* f) {
  double arr1[36], arr2[36];
  matrix_to_array (c.t1_with_RF, arr1);
  matrix_to_array (c.t1_no_RF, arr2);
  lat_param_to_f2_(f, c.n_part, c.total_length, c.unstable_factor,
      arr1, arr2, c.particle, c.ix_lost, c.end_lost_at, c.plane_lost_at,
      c.lattice_type, c.ixx, c.stable, c.aperture_limit_on, c.lost);
}

extern "C" void lat_param_to_c2_(C_lat_param& c, Re& np, Re& total_l, 
      Re& growth_r, ReArr t1_with, ReArr t1_no, Int& part, Int& ixl, Int& end_lost,
      Int& plane_lost, Int& lattice_type, Int& ixx, Int& stable, Int& ap_lim, Int& lost) {
  static Real_Matrix m1(M6_mat), m2(M6_mat);
  m1 << t1_with;
  m2 << t1_no;
  c = C_lat_param(np, total_l, growth_r, m1, m2, part, ixl, 
              end_lost, plane_lost, lattice_type, ixx, stable, ap_lim, lost);
}

void operator>> (C_lat_param& c, lat_param_struct* f) {
  lat_param_to_f_(c, f);
}

void operator>> (lat_param_struct* f, C_lat_param& c) {
  lat_param_to_c_(f, c);
}

//---------------------------------------------------------------------------
// anormal_mode

extern "C" void anormal_mode_to_f2_(anormal_mode_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void anormal_mode_to_f_(C_anormal_mode& c, anormal_mode_struct* f) {
  anormal_mode_to_f2_(f, c.emittance, c.synch_int4, c.synch_int5, c.synch_int6, c.j_damp, 
                                               c.alpha_damp, c.chrom, c.tune);
}

extern "C" void anormal_mode_to_c2_(C_anormal_mode& c, Re& emit, Re& i4, Re& i5, Re& i6,
                                        Re& jd, Re& ad, Re& chrom, Re& tune) {
  c = C_anormal_mode(emit, i4, i5, i6, jd, ad, chrom, tune);
}

void operator>> (C_anormal_mode& c, anormal_mode_struct* f) {
  anormal_mode_to_f_(c, f);
}

void operator>> (anormal_mode_struct* f, C_anormal_mode& c) {
  anormal_mode_to_c_(f, c);
}

//---------------------------------------------------------------------------
// linac_normal_mode

extern "C" void linac_normal_mode_to_f2_(linac_normal_mode_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void linac_normal_mode_to_f_(C_linac_normal_mode& c, linac_normal_mode_struct* f) {
  linac_normal_mode_to_f2_(f, c.i2_E4, c.i3_E7, c.i5a_E6, c.i5b_E6, c.sig_E1, 
                                               c.a_emittance_end, c.b_emittance_end);
}

extern "C" void linac_normal_mode_to_c2_(C_linac_normal_mode& c, Re& i2, Re& i3, Re& i5a, 
                                             Re& i5b, Re& sig_e, Re& ea, Re& eb) {
  c = C_linac_normal_mode(i2, i3, i5a, i5b, sig_e, ea, eb);
}

void operator>> (C_linac_normal_mode& c, linac_normal_mode_struct* f) {
  linac_normal_mode_to_f_(c, f);
}

void operator>> (linac_normal_mode_struct* f, C_linac_normal_mode& c) {
  linac_normal_mode_to_c_(f, c);
}

//---------------------------------------------------------------------------
// normal_modes

extern "C" void normal_modes_to_f2_(normal_modes_struct*, ReArr, Re&, Re&, Re&, Re&, Re&, 
                              C_anormal_mode&, C_anormal_mode&, C_anormal_mode&, C_linac_normal_mode&);

extern "C" void normal_modes_to_f_(C_normal_modes& c, normal_modes_struct* f) {
  normal_modes_to_f2_(f, &c.synch_int[0],
            c.sigE_E, c.sig_z, c.e_loss, c.rf_voltage, c.pz_aperture, c.a, c.b, c.z, c.lin);
}

extern "C" void normal_modes_to_c2_(C_normal_modes& c, ReArr synch_int, Re& sige, 
       Re& sig_z, Re& e_loss, Re& rf_volt, Re& pz, 
       anormal_mode_struct* a, anormal_mode_struct* b, anormal_mode_struct* z, linac_normal_mode_struct* lin) {
  a >> c.a; b >> c.b; z >> c.z; lin >> c.lin;
  c = C_normal_modes(synch_int, sige, sig_z, e_loss, rf_volt, pz, c.a, c.b, c.z, c.lin);
}

void operator>> (C_normal_modes& c, normal_modes_struct* f) {
  normal_modes_to_f_(c, f);
}

void operator>> (normal_modes_struct* f, C_normal_modes& c) {
  normal_modes_to_c_(f, c);
}

//---------------------------------------------------------------------------
// bmad_com

extern "C" void bmad_com_to_f2_(Re&, ReArr, Re&, Re&, Re&, Re&, Re&, Re&, Re&, 
     Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&);

extern "C" void bmad_com_to_f_(C_bmad_com& c) {
  bmad_com_to_f2_(c.max_aperture_limit, &c.d_orb[0],
    c.default_ds_step, c.significant_longitudinal_length, c.rel_tolerance, 
    c.abs_tolerance, c.rel_tol_adaptive_tracking, c.abs_tol_adaptive_tracking, 
    c.taylor_order, c.default_integ_order, c.canonical_coords, 
    c.sr_wakes_on, c.lr_wakes_on, c.mat6_track_symmetric,
    c.auto_bookkeeper, c.trans_space_charge_on, c.coherent_synch_rad_on, 
    c.spin_tracking_on, c.radiation_damping_on, c.radiation_fluctuations_on, 
    c.conserve_taylor_maps);
}

extern "C" void bmad_com_to_c2_(C_bmad_com& c, 
              Re& ap, ReArr orb, Re& ds, Re& significant, Re& rel, Re& abs,  
              Re& rel_adapt, Re& abs_adapt, 
              Int& to, Int& dflt_ord, Int& cc, 
              Int& sr, Int& lr, Int& sym, Int& a_book, Int& tsc_on, Int& csr_on, 
              Int& st_on, Int& rad_d, Int& rad_f, Int& con_t) {
  c.max_aperture_limit               = ap;
  c.d_orb                            << orb;
  c.default_ds_step                  = ds;
  c.significant_longitudinal_length  = significant;
  c.rel_tolerance                    = rel;
  c.abs_tolerance                    = abs;
  c.rel_tol_adaptive_tracking        = rel_adapt;
  c.abs_tol_adaptive_tracking        = abs_adapt;
  c.taylor_order                     = to;
  c.default_integ_order              = dflt_ord;
  c.canonical_coords                 = cc;
  c.sr_wakes_on                      = sr;
  c.lr_wakes_on                      = lr;
  c.mat6_track_symmetric             = sym;
  c.auto_bookkeeper                  = a_book;
  c.trans_space_charge_on            = tsc_on;
  c.coherent_synch_rad_on            = csr_on;
  c.spin_tracking_on                 = st_on;
  c.radiation_damping_on             = rad_d;
  c.radiation_fluctuations_on        = rad_f;
  c.conserve_taylor_maps             = con_t;
}

//---------------------------------------------------------------------------
// em_field

extern "C" void em_field_to_f2_(em_field_struct*, ReArr, ReArr, ReArr, ReArr);

extern "C" void em_field_to_f_(C_em_field& c, em_field_struct* f) {
  double de[9], db[9];
  matrix_to_array (c.dE, de);
  matrix_to_array (c.dB, db);
  em_field_to_f2_(f, &c.E[0], &c.B[0], de, db);
}

extern "C" void em_field_to_c2_(C_em_field& c, ReArr e, ReArr b, ReArr de, ReArr db) {
  c.E << e; c.B << b; 
  c.dE << de; c.dB << db;
}

void operator>> (C_em_field& c, em_field_struct* f) {
  em_field_to_f_(c, f);
}

void operator>> (em_field_struct* f, C_em_field& c) {
  em_field_to_c_(f, c);
}

//---------------------------------------------------------------------------
// ele

extern "C" void ele_to_f2_(ele_struct*, Char, Int&, Char, Int&, Char, Int&, Char, Int&, Char, Int&, 
  C_twiss&, C_twiss&, C_twiss&, C_xy_disp&, C_xy_disp&, 
  C_floor_position&, C_coord&, C_coord&, void*, 
  C_taylor&, C_taylor&, C_taylor&, C_taylor&, C_taylor&, C_taylor&, 
  Int&, ReArr, ReArr, ReArr, ReArr, ReArr,
  Re&, Re&, Re&, ReArr, Int&, Int&, 
  ReArr, ReArr, Int&, ReArr, Int&, 
  Int&, Int&, Int&, Int&, Int&,                        // key
  Int&, Int&, Int&, Int&,                              // slave_status
  Int&, Int&, Int&, Int&,                              // lord_status
  Int&, Int&, Int&, Int&, Int&, Int&, Int&,            // ix_pointer
  Int&, Int&, Int&, Int&, Int&, Int&,                  // symp
  Int&, Int&, Int&, Int&, Int&, Int&,                  // map_with_off
  Int&, Int&, Int&, Int&);

extern "C" void wig_term_in_ele_to_f2_(ele_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Int&);

extern "C" void ele_to_f_(C_ele& c, ele_struct* f) {
  const char* nam = c.name.data();       int n_nam = c.name.length();
  const char* typ = c.type.data();       int n_typ = c.type.length();
  const char* ali = c.alias.data();      int n_ali = c.alias.length();
  const char* descrip = c.descrip.data();    int n_descrip = c.descrip.length();
  const char* component = c.component_name.data(); 
                                    int n_component = c.component_name.length();
  int n_ab = c.a_pole.size(), n_const = c.const_arr.size();
  int n_wig = c.wig_term.size();
  int nr1 = c.r.size(), nr2 = 0;
  if (nr1) {
    nr2 = c.r[0].size();
  }
  double mat6[36], c_mat[4];
  double *r_arr = new double[nr1*nr2];  
  matrix_to_array (c.mat6, mat6);
  matrix_to_array (c.c_mat, c_mat);
  matrix_to_array (c.r, r_arr);
  ele_to_f2_(f, nam, n_nam, typ, n_typ, ali, n_ali, component, n_component, descrip, n_descrip, 
    c.a, c.b, c.z, c.x, c.y, 
    c.floor, c.map_ref_orb_in, c.map_ref_orb_out, c.gen_field,  
    c.taylor[0], c.taylor[1], c.taylor[2], c.taylor[3], c.taylor[4], c.taylor[5], 
    n_wig, &c.value[1], &c.gen0[0], &c.vec0[0], mat6, c_mat,
    c.gamma_c, c.s, c.ref_time, r_arr, nr1, nr2, 
    &c.a_pole[0], &c.b_pole[0], n_ab, &c.const_arr[0], n_const, 
    c.key, c.sub_key, c.ix_ele, c.ix_branch, c.ix_value, 
    c.slave_status, c.n_slave, c.ix1_slave, c.ix2_slave, 
    c.lord_status, c.n_lord, c.ic1_lord, c.ic2_lord, 
    c.ix_pointer, c.ixx, c.mat6_calc_method, c.tracking_method, c.spin_tracking_method, c.field_calc, c.ref_orbit,
    c.aperture_at, c.aperture_type, c.symplectify, c.mode_flip, c.multipoles_on, c.scale_multipoles,
    c.map_with_offsets, c.field_master, c.reversed, c.is_on, c.old_is_on, c.logic, c.bmad_logic,
    c.on_a_girder, c.csr_calc_on, c.offset_moves_aperture);
  for (int i = 0; i < n_wig; i++) {
    wig_term_in_ele_to_f2_(f, i+1, c.wig_term[i].coef, 
            c.wig_term[i].kx, c.wig_term[i].ky, c.wig_term[i].kz, 
            c.wig_term[i].phi_z, c.wig_term[i].type);
  }
  delete[] r_arr;
}

extern "C" void ele_to_c2_(C_ele& c, char* name, char* type, char* alias, 
    char* component, char* descrip, 
    twiss_struct* a, twiss_struct* b, twiss_struct* z, xy_disp_struct* x, xy_disp_struct* y, 
    floor_position_struct* floor, coord_struct* ref_orb_in, coord_struct* ref_orb_out, void* gen,
    taylor_struct* tlr0, taylor_struct* tlr1, taylor_struct* tlr2,
    taylor_struct* tlr3, taylor_struct* tlr4, taylor_struct* tlr5, 
    Int& n_wig, ReArr val, ReArr gen0, ReArr vec0,
    ReArr mat6, ReArr c_mat, Re& gamma_c, Re& s, Re& ref_t, 
    ReArr r_arr, Int& nr1, Int& nr2, ReArr a_pole, ReArr b_pole, Int& n_ab, ReArr const_arr, Int& n_const, 
    Int& key, Int& sub_key, Int& ix_ele, Int& ix_branch, Int& ix_value,  
    Int& slave_status, Int& n_slave, Int& ix1_s, Int& ix2_s, 
    Int& lord_status, Int& n_lord, Int& ic1_l, Int& ic2_l, 
    Int& ix_p, Int& ixx, Int& mat6_calc, Int& tracking, Int& spin_meth, Int& field_calc, Int& ref_orbit,
    Int& aperture_at, Int& aperture_type, 
    Int& symp, Int& mode_flip, Int& multi_on, Int& scale_multi, Int& map_with_off, 
    Int& field_master, Int& reversed, Int& is_on, Int& old_is_on, Int& logic, Int& bmad_logic, Int& on_a_gird, 
    Int& csr_calc, Int& offset_moves_ap) {

  c.name                  = name;
  c.type                  = type;
  c.alias                 = alias; 
  c.component_name        = component;
  c.descrip               = descrip;
  a >> c.a;  b >> c.b;  z >> c.z;
  x >> c.x;  y >> c.y;
  floor >> c.floor;
  ref_orb_in >> c.map_ref_orb_in;
  ref_orb_out >> c.map_ref_orb_out;
  c.gen_field = gen;
  tlr0 >> c.taylor[0];   tlr1 >> c.taylor[1];   tlr2 >> c.taylor[2]; 
  tlr3 >> c.taylor[3];   tlr4 >> c.taylor[4];   tlr5 >> c.taylor[5]; 
  if (c.wig_term.size() != n_wig) c.wig_term.resize(n_wig);
  c.value                 << val;
  c.gen0                  << gen0; 
  c.vec0                  << vec0; 
  c.mat6 << mat6;
  c.c_mat << c_mat;
  c.gamma_c               = gamma_c; 
  c.s                     = s; 
  c.ref_time              = ref_t;

  if (nr1*nr2 == 0) {
    if (!c.r.size()) c.r.resize(0);
  } else {
    if (!(c.r.size() == nr1)) c.r.resize(nr1);
    if (!(c.r[0].size() == nr2)) {
      for (int i = 0; i < nr1; i++) c.r[i].resize(nr2);
    }
    c.r << r_arr;
  }

  if (n_ab > 0) {
    c.a_pole.resize(Bmad::N_POLE_MAXX+1);
    c.b_pole.resize(Bmad::N_POLE_MAXX+1);
    c.a_pole = Real_Array(a_pole, Bmad::N_POLE_MAXX+1);
    c.b_pole = Real_Array(b_pole, Bmad::N_POLE_MAXX+1);
  }

  if (c.const_arr.size() != n_const) c.const_arr.resize(n_const);
  c.const_arr  << const_arr;

  c.key                   = key;
  c.sub_key               = sub_key;
  c.ix_ele                = ix_ele;
  c.ix_branch             = ix_branch;
  c.ix_value              = ix_value;
  c.slave_status          = slave_status;
  c.n_slave               = n_slave;
  c.ix1_slave             = ix1_s;
  c.ix2_slave             = ix2_s;
  c.lord_status           = lord_status;
  c.n_lord                = n_lord;
  c.ic1_lord              = ic1_l;
  c.ic2_lord              = ic2_l;
  c.ix_pointer            = ix_p;
  c.ixx                   = ixx;

  c.mat6_calc_method      = mat6_calc;
  c.tracking_method       = tracking;
  c.spin_tracking_method  = spin_meth;
  c.field_calc            = field_calc;
  c.ref_orbit             = ref_orbit;
  c.aperture_at           = aperture_at;
  c.aperture_type         = aperture_type;
  c.symplectify           = symp;
  c.mode_flip             = mode_flip;
  c.multipoles_on         = multi_on;
  c.scale_multipoles      = scale_multi;
  c.map_with_offsets      = map_with_off;
  c.field_master          = field_master;
  c.reversed              = reversed;
  c.is_on                 = is_on;
  c.old_is_on             = old_is_on;
  c.logic                 = logic;
  c.bmad_logic            = bmad_logic;
  c.on_a_girder           = on_a_gird;
  c.csr_calc_on           = csr_calc;
  c.offset_moves_aperture = offset_moves_ap;
}

extern "C" void wig_term_in_ele_to_c2_(C_ele& c, Int& it, 
                  Re& coef, Re& kx, Re& ky, Re& kz, Re& phi_z, Int& type) {
  c.wig_term[it-1] = C_wig_term(coef, kx, ky, kz, phi_z, type);
}

void operator>> (C_ele& c, ele_struct* f) {
  ele_to_f_(c, f);
}

void operator>> (ele_struct* f, C_ele& c) {
  ele_to_c_(f, c);
}

C_ele& C_ele::operator= (const C_ele& c) {
  name                  = c.name;
  type                  = c.type;
  alias                 = c.alias;
  component_name        = c.component_name;
  descrip               = c.descrip;
  a                     = c.a;
  b                     = c.b;
  z                     = c.z;
  x                     = c.x;
  y                     = c.y;
  floor                 = c.floor;
  map_ref_orb_in        = c.map_ref_orb_in;
  map_ref_orb_out       = c.map_ref_orb_out;
  gen_field             = c.gen_field;
  taylor                << c.taylor;
  wig_term              << c.wig_term;
  value                 << c.value;
  gen0                  << c.gen0;
  vec0                  << c.vec0;
  mat6                  << c.mat6;
  c_mat                 << c.c_mat;
  gamma_c               = c.gamma_c;
  s                     = c.s;
  ref_time              = c.ref_time;
  r                     << c.r;
  a_pole                << c.a_pole;
  b_pole                << c.b_pole;
  const_arr             << c.const_arr;
  key                   = c.key;
  sub_key               = c.sub_key;
  ix_ele                = c.ix_ele;
  ix_branch             = c.ix_branch;
  ix_value              = c.ix_value;
  slave_status          = c.slave_status;
  n_slave               = c.n_slave;
  ix1_slave             = c.ix1_slave;
  ix2_slave             = c.ix2_slave;
  lord_status           = c.lord_status;
  n_lord                = c.n_lord;
  ic1_lord              = c.ic1_lord;
  ic2_lord              = c.ic2_lord;
  ix_pointer            = c.ix_pointer;
  ixx                   = c.ixx;

  mat6_calc_method      = c.mat6_calc_method;
  tracking_method       = c.tracking_method;
  field_calc            = c.field_calc;
  ref_orbit             = c.ref_orbit;
  aperture_at           = c.aperture_at;
  aperture_type        = c.aperture_type;
  symplectify           = c.symplectify;
  mode_flip             = c.mode_flip;
  multipoles_on         = c.multipoles_on;
  scale_multipoles      = c.scale_multipoles;
  map_with_offsets      = c.map_with_offsets;
  field_master          = c.field_master;
  reversed              = c.reversed;
  is_on                 = c.is_on;
  old_is_on             = c.old_is_on;
  logic                 = c.logic;
  on_a_girder           = c.on_a_girder;
  csr_calc_on           = c.csr_calc_on;
  offset_moves_aperture = c.offset_moves_aperture;

  return *this;
}

//---------------------------------------------------------------------------
// mode_info

extern "C" void mode_info_to_f2_(mode_info_struct*, Re&, Re&, Re&, Re&, Re&);

extern "C" void mode_info_to_f_(C_mode_info& c, mode_info_struct* f) {
  mode_info_to_f2_(f, c.tune, c.emit, c.chrom, c.sigma, c.sigmap);
}

extern "C" void mode_info_to_c2_(C_mode_info& c, Re& tune, Re& emit, Re& chrom,
                               Re& sigma, Re& sigmap) {
  c = C_mode_info(tune, emit, chrom, sigma, sigmap);
}

void operator>> (C_mode_info& c, mode_info_struct* f) {
  mode_info_to_f_(c, f);
}

void operator>> (mode_info_struct* f, C_mode_info& c) {
  mode_info_to_c_(f, c);
}

//---------------------------------------------------------------------------
// lat

extern "C" void lat_to_f2_(lat_struct*, Char, Int&, Char, Int&, Char, Int&,
    Char, Int&, C_mode_info&, C_mode_info&, C_mode_info&, C_lat_param&, 
    Int&, Int&, Int&, Int&, Int&, Int&, Int&,
    C_ele&, Int&, IntArr, Int&);

extern "C" void ele_from_lat_to_f2_(lat_struct*, Int&, C_ele&);
extern "C" void control_from_lat_to_f2_(lat_struct*, Int&, C_control&);


extern "C" void lat_to_f_(C_lat& c, lat_struct* f) {
  const char* use_name  = c.use_name.data();    int n_name = c.use_name.size();
  const char* lat   = c.lattice.data();         int n_lat = c.lattice.size();
  const char* file  = c.input_file_name.data(); int n_file = c.input_file_name.size();
  const char* title = c.title.data();           int n_title = c.title.size();
  int n_con = c.control.size();
  int n_ic  = c.ic.size();
  int n_ele_max = c.ele.size() - 1;
  lat_to_f2_(f, use_name, n_name, lat, n_lat, file, n_file, title, n_title,
      c.a, c.b, c.z, c.param, 
      c.version, c.n_ele_track, c.n_ele_max, n_ele_max, c.n_control_max, 
      c.n_ic_max, c.input_taylor_order, c.ele_init, n_con, &c.ic[0], n_ic);
  for (int i = 0; i < n_ele_max+1; i++) {
    ele_from_lat_to_f2_(f, i, c.ele[i]);
  }
  for (int i = 0; i < n_con; i++) {
    control_from_lat_to_f2_(f, i+1, c.control[i]);
  }
}

extern "C" void lat_to_c2_(C_lat& c, char* use_name, char* lat, char* file,
    char* title, mode_info_struct* a, mode_info_struct* b, mode_info_struct* z,
    lat_param_struct* param, Int& ver, Int& n_track, Int& n_max, Int& n_maxx, 
    Int& n_con_max, Int& n_ic_max, Int& n_taylor, 
    ele_struct* ele_init, Int& n_con_array, IntArr ic, Int& n_ic_array) {
  c.use_name            = use_name;
  c.lattice             = lat;
  c.input_file_name     = file;
  c.title               = title;
  c.version             = ver;
  c.n_ele_track         = n_track;
  c.n_ele_max           = n_max;
  c.n_control_max       = n_con_max;
  c.n_ic_max            = n_ic_max;
  c.input_taylor_order  = n_taylor;

  if (c.control.size() != n_con_array) c.control.resize(n_con_array);
  if (!(c.ele.size() == n_maxx)) c.ele.resize(n_maxx);
  if (!(c.ic.size() == n_ic_array)) c.ic.resize(n_ic_array);

  c.ic << ic;
  a >> c.a;  b >> c.b;  z >> c.z;
  param >> c.param;
  ele_init >> c.ele_init;
}


extern "C" void ele_from_lat_to_c2_(C_lat& lat, Int& it, ele_struct* ele) {
  ele >> lat.ele[it];
}

extern "C" void control_from_lat_to_c2_(C_lat& c, Int& it, control_struct* control) {
  control >> c.control[it-1];
}

void operator>> (C_lat& c, lat_struct* f) {
  lat_to_f_(c, f);
}

void operator>> (lat_struct* f, C_lat& c) {
  lat_to_c_(f, c);
}

C_lat& C_lat::operator= (const C_lat& c) {
  use_name           = c.use_name;
  lattice            = c.lattice;
  input_file_name    = c.input_file_name;
  title              = c.title;
  a                  = c.a;
  b                  = c.b;
  z                  = c.z;
  param              = c.param;
  version            = c.version;
  n_ele_track        = c.n_ele_track;
  n_ele_max          = c.n_ele_max;    
  n_control_max      = c.n_control_max;        
  n_ic_max           = c.n_ic_max;   
  input_taylor_order = c.input_taylor_order;             
  ele_init           = c.ele_init;   

  if (ele.size() < n_ele_max+1) ele.resize(n_ele_max+1);
  for (int i = 0; i < n_ele_max+1; i++) ele[i] = c.ele[i];

  int n_control = c.control.size();
  if (control.size() < n_control) control.resize(n_control);
  for (int i = 0; i < n_control; i++) control[i] = c.control[i];
  
  int n_ic = c.ic.size();
  if (ic.size() < n_ic) ic.resize(n_ic);
  for (int i = 0; i < n_ic; i++) ic[i] = c.ic[i];

  return *this;  
}
