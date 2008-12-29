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
 void operator<< (Real_Array&, const double*);
 void operator<< (Real_Matrix&, const double*);
 void operator<< (Int_Array&, const int*);

 void operator<< (Real_Array&, const Real_Array&);
 void operator<< (Real_Matrix&, const Real_Matrix&);
 void operator<< (Int_Array&, const Int_Array&);

 void operator<< (C_taylor_array&, const C_taylor_array&);
 void operator<< (C_wig_term_array&, const C_wig_term_array&);

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

extern "C" void coord_to_f2_(coord_struct*, ReArr);

extern "C" void coord_to_f_(C_coord& c, coord_struct* f) {
  coord_to_f2_(f, &c.vec[0]);
}

extern "C" void coord_to_c2_(C_coord& c, double vec[]) {
  c.vec << vec;
}

void operator>> (C_coord& c, coord_struct* f) {
  coord_to_f_(c, f);
}

void operator>> (coord_struct* f, C_coord& c) {
  coord_to_c_(f, c);
}

//---------------------------------------------------------------------------
// Twiss

extern "C" void twiss_to_f2_(twiss_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void twiss_to_f_(const C_twiss& c, twiss_struct* f) {
  twiss_to_f2_(f, c.beta, c.alpha, c.gamma, c.phi, c.eta, c.etap, c.sigma, c.sigma_p, c.emit);
}

extern "C" void twiss_to_c2_(C_twiss& c, Re& beta, Re& alpha, Re& gamma, Re& phi, 
                          Re& eta, Re& etap, Re& sigma, Re& sigma_p, Re& emit) {
  c = C_twiss(beta, alpha, gamma, phi, eta, etap, sigma, sigma_p, emit);
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

extern "C" void xy_disp_to_f_(const C_xy_disp& c, xy_disp_struct* f) {
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
  taylor_term_to_f2_(f, c.coef, &c.exp[0]);
}

extern "C" void taylor_term_to_c2_(C_taylor_term& c, Re& coef, int exp[]) {
  c = C_taylor_term(coef, exp);
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
    taylor_term_in_taylor_to_f2_(f, i+1, c.term[i].coef, &c.term[i].exp[0]);
  }
}

extern "C" void taylor_to_c2_(C_taylor& c, Int& n_term, Re& ref) {
  if (c.term.size() != n_term) c.term.resize(n_term);
  c.ref = ref;
}

extern "C" void taylor_term_in_taylor_to_c2_(C_taylor& c, Int& it, 
                                                       Re& coef, int exp[]) {
  c.term[it-1] = C_taylor_term(coef, exp);
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
// sr_table_wake

extern "C" void sr_table_wake_to_f2_(sr_table_wake_struct*, Re&, Re&, Re&);

extern "C" void sr_table_wake_to_f_(C_sr_table_wake& c, sr_table_wake_struct* f) {
  sr_table_wake_to_f2_(f, c.z, c.longitudinal, c.transverse);
}

extern "C" void sr_table_wake_to_c2_(C_sr_table_wake& c, Re& z, Re& lw, Re& tw) {
  c = C_sr_table_wake(z, lw, tw);
}

void operator>> (C_sr_table_wake& c, sr_table_wake_struct* f) {
  sr_table_wake_to_f_(c, f);
}

void operator>> (sr_table_wake_struct* f, C_sr_table_wake& c) {
  sr_table_wake_to_c_(f, c);
}

//---------------------------------------------------------------------------
// sr_mode_wake

extern "C" void sr_mode_wake_to_f2_(sr_mode_wake_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void sr_mode_wake_to_f_(C_sr_mode_wake& c, sr_mode_wake_struct* f) {
  sr_mode_wake_to_f2_(f, c.amp, c.damp, c.k, c.phi,
                          c.norm_sin, c.norm_cos, c.skew_sin, c.skew_cos);
}

extern "C" void sr_mode_wake_to_c2_(C_sr_mode_wake& c, Re& amp, Re& damp, Re& k, Re& phi,
                          Re& norm_sin, Re& norm_cos, Re& skew_sin, Re& skew_cos) {
  c = C_sr_mode_wake(amp, damp, k, phi, norm_sin, norm_cos, skew_sin, skew_cos);
}

void operator>> (C_sr_mode_wake& c, sr_mode_wake_struct* f) {
  sr_mode_wake_to_f_(c, f);
}

void operator>> (sr_mode_wake_struct* f, C_sr_mode_wake& c) {
  sr_mode_wake_to_c_(f, c);
}

//---------------------------------------------------------------------------
// lr_wake

extern "C" void lr_wake_to_f2_(lr_wake_struct*, Re&, Re&, Re&, Re&, Re&,
                                      Re&, Re&, Re&, Re&, Int&, Int&);

extern "C" void lr_wake_to_f_(C_lr_wake& c, lr_wake_struct* f) {
  lr_wake_to_f2_(f, c.freq, c.freq_in, c.R_over_Q, c.Q, c.angle, 
          c.norm_sin, c.norm_cos, c.skew_sin, c.skew_cos, c.m, c.polarized);
}

extern "C" void lr_wake_to_c2_(C_lr_wake& c, Re& freq, Re& freq_in, 
                  Re& R_over_Q, Re& Q, Re& ang, Re& n_sin, Re& n_cos, 
                  Re& s_sin, Re& s_cos, Int& m, Int& pol) {
  c = C_lr_wake(freq, freq_in, R_over_Q, Q, ang, 
                    n_sin, n_cos, s_sin, s_cos, m, pol);
}

void operator>> (C_lr_wake& c, lr_wake_struct* f) {
  lr_wake_to_f_(c, f);
}

void operator>> (lr_wake_struct* f, C_lr_wake& c) {
  lr_wake_to_c_(f, c);
}

//---------------------------------------------------------------------------
// wake

extern "C" void wake_to_f2_(wake_struct*, Char, Int&, Char, Int&, Re&, Int&, Int&, Int&, Int&);
extern "C" void sr_table_wake_in_wake_to_f2_(wake_struct*, Int&, Re&, Re&, Re&);
extern "C" void sr_mode_long_wake_in_wake_to_f2_(wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);
extern "C" void sr_mode_trans_wake_in_wake_to_f2_(wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Re&, Re&, Re&);
extern "C" void lr_wake_in_wake_to_f2_(wake_struct*, Int&, Re&, Re&, Re&, Re&, Re&,
                                                         Re&, Re&, Re&, Re&, Int&, Int&);

extern "C" void wake_to_f_(C_wake& c, wake_struct* f) {
  int n_lr = c.lr.size();
  int n_sr_table = c.sr_table.size(); 
  int n_sr_mode_long = c.sr_mode_long.size(); 
  int n_sr_mode_trans = c.sr_mode_trans.size(); 
  const char* srf = c.sr_file.data();     int n_srf = c.sr_file.length();
  const char* lrf = c.lr_file.data();     int n_lrf = c.lr_file.length();
  wake_to_f2_(f, srf, n_srf, lrf, n_lrf, c.z_sr_mode_max, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr);
  for (int i = 0; i < n_sr_table; i++) {
    sr_table_wake_in_wake_to_f2_(f, i, c.sr_table[i].z, c.sr_table[i].longitudinal, c.sr_table[i].transverse);
  }
  for (int i = 0; i < n_sr_mode_long; i++) {
    sr_mode_long_wake_in_wake_to_f2_(f, i+1, c.sr_mode_long[i].amp, c.sr_mode_long[i].damp, 
        c.sr_mode_long[i].k, c.sr_mode_long[i].phi, c.sr_mode_long[i].norm_sin, 
        c.sr_mode_long[i].norm_cos, c.sr_mode_long[i].skew_sin, c.sr_mode_long[i].skew_cos);
  }
  for (int i = 0; i < n_sr_mode_trans; i++) {
    sr_mode_trans_wake_in_wake_to_f2_(f, i+1, c.sr_mode_trans[i].amp, c.sr_mode_trans[i].damp, 
        c.sr_mode_trans[i].k, c.sr_mode_trans[i].phi, c.sr_mode_trans[i].norm_sin, 
        c.sr_mode_trans[i].norm_cos, c.sr_mode_trans[i].skew_sin, c.sr_mode_trans[i].skew_cos);
  }
  for (int i = 0; i < n_lr; i++) {
    lr_wake_in_wake_to_f2_(f, i+1, c.lr[i].freq, c.lr[i].freq_in, 
        c.lr[i].R_over_Q, c.lr[i].Q, c.lr[i].angle, c.lr[i].norm_sin, 
        c.lr[i].norm_cos, c.lr[i].skew_sin, c.lr[i].skew_cos, c.lr[i].m, c.lr[i].polarized);
  }
}

extern "C" void wake_to_c2_(C_wake& c, char* srf, char* lrf, Re& z_cut, Int& n_sr_table, Int& n_sr_mode_long,
                            Int& n_sr_mode_trans, Int& n_lr) {
  if (c.sr_table.size() != n_sr_table) c.sr_table.resize(n_sr_table);
  if (c.sr_mode_long.size() != n_sr_mode_long) c.sr_mode_long.resize(n_sr_mode_long);
  if (c.sr_mode_trans.size() != n_sr_mode_trans) c.sr_mode_trans.resize(n_sr_mode_trans);
  c.sr_file = srf;
  if (c.lr.size() != n_lr) c.lr.resize(n_lr);
  c.lr_file = lrf;
  c.z_sr_mode_max = z_cut;
}

extern "C" void sr_table_wake_in_wake_to_c2_(C_wake& c, Int& it, Re& z, Re& l, Re& t) {
  c.sr_table[it] = C_sr_table_wake(z, l, t);
}

extern "C" void sr_mode_long_wake_in_wake_to_c2_(C_wake& c, Int& it, Re& a, Re& d, 
                      Re& f, Re& p, Re& ns, Re& nc, Re& ss, Re& sc) {
  c.sr_mode_long[it-1] = C_sr_mode_wake(a, d, f, p, ns, nc, ss, sc);
}

extern "C" void sr_mode_trans_wake_in_wake_to_c2_(C_wake& c, Int& it, Re& a, Re& d, 
                      Re& f, Re& p, Re& ns, Re& nc, Re& ss, Re& sc) {
  c.sr_mode_trans[it-1] = C_sr_mode_wake(a, d, f, p, ns, nc, ss, sc);
}

extern "C" void lr_wake_in_wake_to_c2_(C_wake& c, Int& it, Re& f, Re& k, 
           Re& i, Re& q, Re& ang, Re& ns, Re& nc, Re& ss, Re& sc, Int& m, Int& pol) {
  c.lr[it-1] = C_lr_wake(f, k, i, q, ang, ns, nc, ss, sc, m, pol);
}


void operator>> (C_wake& c, wake_struct* f) {
  wake_to_f_(c, f);
}

void operator>> (wake_struct* f, C_wake& c) {
  wake_to_c_(f, c);
}

C_wake& C_wake::operator= (const C_wake& c) {
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
  control_to_f2_(f, c.coef, c.ix_lord, c.ix_slave, c.ix_photon_line, c.ix_attrib);
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
// param

extern "C" void param_to_f2_(lat_param_struct*, Re&, Re&, Re&, ReArr, ReArr, 
                               Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&);

extern "C" void param_to_f_(C_param& c, lat_param_struct* f) {
  double arr1[36], arr2[36];
  matrix_to_array (c.t1_with_RF, arr1);
  matrix_to_array (c.t1_no_RF, arr2);
  param_to_f2_(f, c.n_part, c.total_length, c.growth_rate,
      arr1, arr2, c.particle, c.ix_lost, c.end_lost_at, c.plane_lost_at,
      c.lattice_type, c.ixx, c.stable, c.aperture_limit_on, c.lost);
}

extern "C" void param_to_c2_(C_param& c, Re& np, Re& total_l, 
      Re& growth_r, ReArr t1_with, ReArr t1_no, Int& part, Int& ixl, Int& end_lost,
      Int& plane_lost, Int& lattice_type, Int& ixx, Int& stable, Int& ap_lim, Int& lost) {
  static Real_Matrix m1(M6_mat), m2(M6_mat);
  m1 << t1_with;
  m2 << t1_no;
  c = C_param(np, total_l, growth_r, m1, m2, part, ixl, 
              end_lost, plane_lost, lattice_type, ixx, stable, ap_lim, lost);
}

void operator>> (C_param& c, lat_param_struct* f) {
  param_to_f_(c, f);
}

void operator>> (lat_param_struct* f, C_param& c) {
  param_to_c_(f, c);
}

//---------------------------------------------------------------------------
// amode

extern "C" void amode_to_f2_(anormal_mode_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void amode_to_f_(C_amode& c, anormal_mode_struct* f) {
  amode_to_f2_(f, c.emittance, c.synch_int4, c.synch_int5, c.j_damp, 
                                               c.alpha_damp, c.chrom, c.tune);
}

extern "C" void amode_to_c2_(C_amode& c, Re& emit, Re& i4, Re& i5, 
                                        Re& jd, Re& ad, Re& chrom, Re& tune) {
  c = C_amode(emit, i4, i5, jd, ad, chrom, tune);
}

void operator>> (C_amode& c, anormal_mode_struct* f) {
  amode_to_f_(c, f);
}

void operator>> (anormal_mode_struct* f, C_amode& c) {
  amode_to_c_(f, c);
}

//---------------------------------------------------------------------------
// linac_mode

extern "C" void linac_mode_to_f2_(linac_normal_mode_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&);

extern "C" void linac_mode_to_f_(C_linac_mode& c, linac_normal_mode_struct* f) {
  linac_mode_to_f2_(f, c.i2_E4, c.i3_E7, c.i5a_E6, c.i5b_E6, c.sig_E1, 
                                               c.a_emittance_end, c.b_emittance_end);
}

extern "C" void linac_mode_to_c2_(C_linac_mode& c, Re& i2, Re& i3, Re& i5a, 
                                             Re& i5b, Re& sig_e, Re& ea, Re& eb) {
  c = C_linac_mode(i2, i3, i5a, i5b, sig_e, ea, eb);
}

void operator>> (C_linac_mode& c, linac_normal_mode_struct* f) {
  linac_mode_to_f_(c, f);
}

void operator>> (linac_normal_mode_struct* f, C_linac_mode& c) {
  linac_mode_to_c_(f, c);
}

//---------------------------------------------------------------------------
// modes

extern "C" void modes_to_f2_(normal_modes_struct*, Re&, Re&, Re&, Re&, Re&, Re&, Re&,
                              C_amode, C_amode, C_amode, C_linac_mode);

extern "C" void modes_to_f_(C_modes& c, normal_modes_struct* f) {
  modes_to_f2_(f, c.synch_int1, c.synch_int2, c.synch_int3,
            c.sigE_E, c.sig_z, c.e_loss, c.pz_aperture, c.a, c.b, c.z, c.lin);
}

extern "C" void modes_to_c2_(C_modes& c, Re& i1, Re& i2, Re& i3, Re& sige, 
       Re& sig_z, Re& e_loss, Re& pz, 
       C_amode a, C_amode b, C_amode z, C_linac_mode lin) {
  c = C_modes(i1, i2, i3, sige, sig_z, e_loss, pz, a, b, z, lin);
}

void operator>> (C_modes& c, normal_modes_struct* f) {
  modes_to_f_(c, f);
}

void operator>> (normal_modes_struct* f, C_modes& c) {
  modes_to_c_(f, c);
}

//---------------------------------------------------------------------------
// bmad_com

extern "C" void bmad_com_to_f2_(Re&, ReArr, Re&, Re&, Re&, Re&, Re&, Re&, 
     Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, 
     Int&, Int&, Int&, Int&);

extern "C" void bmad_com_to_f_(C_bmad_com& c) {
  bmad_com_to_f2_(c.max_aperture_limit, &c.d_orb[0], c.grad_loss_sr_wake, 
    c.default_ds_step, c.rel_tolerance, c.abs_tolerance, 
    c.rel_tol_adaptive_tracking, c.abs_tol_adaptive_tracking, 
    c.taylor_order, 
    c.default_integ_order, c.canonical_coords, 
    c.use_liar_lcavity, c.sr_wakes_on, c.lr_wakes_on, c.mat6_track_symmetric,
    c.auto_bookkeeper, c.trans_space_charge_on, c.coherent_synch_rad_on, 
    c.spin_tracking_on, c.radiation_damping_on, c.radiation_fluctuations_on, 
    c.compute_ref_energy, c.conserve_taylor_maps);
}

extern "C" void bmad_com_to_c2_(C_bmad_com& c, 
              Re& ap, ReArr orb, Re& kl, Int& ds, Re& rel, Re& abs, Re& rel_adapt, Re& abs_adapt, 
              Int& to, Int& dflt_ord, Int& cc, Int& liar, 
              Int& sr, Int& lr, Int& sym, Int& a_book, Int& tsc_on, Int& csr_on, 
              Int& st_on, Int& rad_d, Int& rad_f, Int& ref_e, Int& con_t) {
  c = C_bmad_com (ap, orb, kl, ds, rel, abs, rel_adapt, abs_adapt, to, dflt_ord, cc, liar, sr, 
             lr, sym, a_book, tsc_on, csr_on, st_on, rad_d, rad_f, ref_e, con_t);
}

//---------------------------------------------------------------------------
// em_field

extern "C" void em_field_to_f2_(em_field_struct*, ReArr, ReArr, ReArr, ReArr, 
                ReArr, ReArr, Int&);

extern "C" void em_field_to_f_(C_em_field& c, em_field_struct* f) {
  double de[9], db[9], dk[9];
  matrix_to_array (c.dE, de);
  matrix_to_array (c.dB, db);
  matrix_to_array (c.dkick, dk);
  em_field_to_f2_(f, &c.E[0], &c.B[0], &c.kick[0], de, db, dk, c.type);
}

extern "C" void em_field_to_c2_(C_em_field& c, ReArr e, ReArr b, ReArr k, 
                                    ReArr de, ReArr db, ReArr dk, Int& tp) {
  c.E << e; c.B << b; c.kick << k;
  c.dE << de; c.dB << db; c.dkick << dk; c.type = tp;
}

void operator>> (C_em_field& c, em_field_struct* f) {
  em_field_to_f_(c, f);
}

void operator>> (em_field_struct* f, C_em_field& c) {
  em_field_to_c_(f, c);
}

//---------------------------------------------------------------------------
// ele

extern "C" void ele_to_f2_(ele_struct*, Char, Int&, Char, Int&, Char, Int&, Char, 
  Int&, C_xy_disp, C_xy_disp, C_twiss&, C_twiss&, C_twiss&, C_floor_position&, 
  ReArr, ReArr, ReArr, 
  ReArr, ReArr, Re&, Re&, ReArr, Int&, Int&, ReArr, ReArr, Int&, ReArr, Int&, 
  Char, Int&, void*, C_taylor&, C_taylor&, C_taylor&, C_taylor&, C_taylor&, 
  C_taylor&, C_wake&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, 
  Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&,
  Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&, Int&);

extern "C" void wig_term_in_ele_to_f2_(ele_struct*, Int&, Re&, Re&, Re&, Re&, Re&, Int&);

extern "C" void ele_to_f_(C_ele& c, ele_struct* f) {
  const char* nam = c.name.data();       int n_nam = c.name.length();
  const char* typ = c.type.data();       int n_typ = c.type.length();
  const char* ali = c.alias.data();      int n_ali = c.alias.length();
  const char* des = c.descrip.data();    int n_des = c.descrip.length();
  const char* attrib = c.attribute_name.data(); 
                                    int n_attrib = c.attribute_name.length();
  int n_ab = c.a_pole.size(), n_const = c.const_arr.size();
  int n_sr_table = c.wake.sr_table.size(), n_sr_mode_long = c.wake.sr_mode_long.size();
  int n_sr_mode_trans = c.wake.sr_mode_trans.size(), n_lr = c.wake.lr.size();
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
  ele_to_f2_(f, nam, n_nam, typ, n_typ, ali, n_ali, attrib, n_attrib,
    c.x, c.y, c.a, c.b, c.z, c.floor, &c.value[1], &c.gen0[0], &c.vec0[0], mat6, c_mat, 
    c.gamma_c, c.s, r_arr, nr1, nr2, &c.a_pole[0], &c.b_pole[0], n_ab, &c.const_arr[0], 
    n_const, des, n_des, c.gen_field, c.taylor[0], c.taylor[1], c.taylor[2], 
    c.taylor[3], c.taylor[4], c.taylor[5], c.wake, n_sr_table, n_sr_mode_long, 
    n_sr_mode_trans, n_lr, n_wig, c.key, 
    c.sub_key, c.control_type, c.ix_value, c.n_slave, c.ix1_slave, 
    c.ix2_slave, c.n_lord, c.ic1_lord, c.ic2_lord, c.ix_pointer, 
    c.ixx, c.ix_ele, c.ix_photon_line, c.mat6_calc_method, c.tracking_method, 
    c.field_calc, c.num_steps, c.integrator_order, c.ref_orbit, c.taylor_order, 
    c.aperture_at, c.coupler_at, c.symplectify, c.mode_flip, c.multipoles_on, 
    c.map_with_offsets, c.field_master, c.is_on, c.old_is_on, c.logic, c.on_a_girder, 
    c.csr_calc_on, c.offset_moves_aperture);
  for (int i = 0; i < n_wig; i++) {
    wig_term_in_ele_to_f2_(f, i+1, c.wig_term[i].coef, 
            c.wig_term[i].kx, c.wig_term[i].ky, c.wig_term[i].kz, 
            c.wig_term[i].phi_z, c.wig_term[i].type);
  }
  delete[] r_arr;
}

extern "C" void ele_to_c2_(C_ele& c, char* name, char* type, char* alias,
    char* attrib, xy_disp_struct* x, xy_disp_struct* y, 
    twiss_struct* a, twiss_struct* b, twiss_struct* z,
    floor_position_struct* floor, ReArr val, ReArr gen0, ReArr vec0,
    ReArr mat6, ReArr c_mat, Re& gamma_c, Re& s, ReArr r_arr, Int& nr1, 
    Int& nr2, ReArr a_pole, ReArr b_pole, 
    Int& n_ab, ReArr const_arr, Int& n_const, char* descrip, 
    void* gen, taylor_struct* tlr0, taylor_struct* tlr1, taylor_struct* tlr2,
    taylor_struct* tlr3, taylor_struct* tlr4, taylor_struct* tlr5, 
    wake_struct* wake, Int& wake_here, Int& n_wig, Int& key, Int& sub_key, 
    Int& control, Int& ix_v, Int& n_s, Int& ix1_s, Int& ix2_s, Int& n_l, 
    Int& ic1_l, Int& ic2_l, Int& ix_p, Int& ixx, Int& ix_e, Int& ix_photon, 
    Int& mat6_calc, Int& tracking, Int& f_calc, Int& num_s, Int& int_ord, 
    Int& ptc, Int& t_ord, Int& aperture_at, Int& coupler_at, Int& symp, 
    Int& flip, Int& multi, Int& rad, Int& f_master, Int& is_on, 
    Int& internal, Int& logic, Int& girder, Int& csr_calc, Int& offset_moves_ap) {

  c.name                  = name;
  c.type                  = type;
  c.alias                 = alias; 
  c.attribute_name        = attrib;
  c.value                 << val;
  c.gen0                  << gen0; 
  c.vec0                  << vec0; 
  c.gamma_c               = gamma_c; 
  c.s                     = s; 
  c.descrip               = descrip;
  c.key                   = key;
  c.sub_key               = sub_key;
  c.control_type          = control;
  c.ix_value              = ix_v;
  c.n_slave               = n_s;
  c.ix1_slave             = ix1_s;
  c.ix2_slave             = ix2_s;
  c.n_lord                = n_l;
  c.ic1_lord              = ic1_l;
  c.ic2_lord              = ic2_l;
  c.ix_pointer            = ix_p;
  c.ixx                   = ixx;
  c.ix_ele                = ix_e;
  c.ix_photon_line        = ix_photon;
  c.mat6_calc_method      = mat6_calc;
  c.tracking_method       = tracking;
  c.field_calc            = f_calc;
  c.num_steps             = num_s;
  c.integrator_order      = int_ord;
  c.ref_orbit              = ptc;
  c.taylor_order          = t_ord;
  c.aperture_at           = aperture_at;
  c.coupler_at            = coupler_at;
  c.symplectify           = symp;
  c.mode_flip             = flip;
  c.multipoles_on         = multi;
  c.map_with_offsets      = rad;
  c.field_master          = f_master;
  c.is_on                 = is_on;
  c.old_is_on             = internal;
  c.logic                 = logic;
  c.on_a_girder           = girder;
  c.csr_calc_on           = csr_calc;
  c.offset_moves_aperture = offset_moves_ap;

  if (c.const_arr.size() != n_const) c.const_arr.resize(n_const);
  c.const_arr  << const_arr;

  if (c.wig_term.size() != n_wig) c.wig_term.resize(n_wig);

  if (nr1*nr2 == 0) {
    if (!c.r.size()) c.r.resize(0);
  } else {
    if (!(c.r.size() == nr1)) c.r.resize(nr1);
    if (!(c.r[0].size() == nr2)) {
      for (int i = 0; i < nr1; i++) c.r[i].resize(nr2);
    }
    c.r << r_arr;
  }

  x >> c.x;  y >> c.y;
  a >> c.a;  b >> c.b;  z >> c.z;
  floor >> c.floor;

  tlr0 >> c.taylor[0];   tlr1 >> c.taylor[1];   tlr2 >> c.taylor[2]; 
  tlr3 >> c.taylor[3];   tlr4 >> c.taylor[4];   tlr5 >> c.taylor[5]; 
  c.mat6 << mat6;
  c.c_mat << c_mat;
  if (wake_here) wake >> c.wake;  
  c.gen_field = gen;
  if (n_ab > 0) {
    c.a_pole.resize(Bmad::N_POLE_MAXX+1);
    c.b_pole.resize(Bmad::N_POLE_MAXX+1);
    c.a_pole = Real_Array(a_pole, Bmad::N_POLE_MAXX+1);
    c.b_pole = Real_Array(b_pole, Bmad::N_POLE_MAXX+1);
  }
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
  name                 = c.name;
  type                 = c.type;
  alias                = c.alias;
  attribute_name       = c.attribute_name;
  x                    = c.x;
  y                    = c.y;
  a                    = c.a;
  b                    = c.b;
  z                    = c.z;
  floor                = c.floor;
  value                << c.value;
  gen0                 << c.gen0;
  vec0                 << c.vec0;
  mat6                 << c.mat6;
  c_mat                << c.c_mat;
  gamma_c              = c.gamma_c;
  s                    = c.s;
  r                    << c.r;
  a_pole               << c.a_pole;
  b_pole               << c.b_pole;
  const_arr            << c.const_arr;
  descrip              = c.descrip;
  gen_field            = c.gen_field;
  taylor               << c.taylor;
  wig_term             << c.wig_term;
  wake                 = c.wake;
  key                  = c.key;
  sub_key              = c.sub_key;
  control_type         = c.control_type;
  ix_value             = c.ix_value;
  n_slave              = c.n_slave;
  ix1_slave            = c.ix1_slave;
  ix2_slave            = c.ix2_slave;
  n_lord               = c.n_lord;
  ic1_lord             = c.ic1_lord;
  ic2_lord             = c.ic2_lord;
  ix_pointer           = c.ix_pointer;
  ixx                  = c.ixx;
  ix_ele               = c.ix_ele;
  ix_photon_line       = c.ix_photon_line;
  mat6_calc_method     = c.mat6_calc_method;
  tracking_method      = c.tracking_method;
  field_calc           = c.field_calc;
  num_steps            = c.num_steps;
  integrator_order      = c.integrator_order;
  ref_orbit             = c.ref_orbit;
  taylor_order         = c.taylor_order;
  aperture_at          = c.aperture_at;
  symplectify          = c.symplectify;
  mode_flip            = c.mode_flip;
  multipoles_on        = c.multipoles_on;
  map_with_offsets     = c.map_with_offsets;
  field_master         = c.field_master;
  is_on                = c.is_on;
  old_is_on            = c.old_is_on;
  logic                = c.logic;
  on_a_girder          = c.on_a_girder;
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
    Char, Int&, C_mode_info&, C_mode_info&, C_mode_info&, C_param&, 
    Int&, Int&, Int&, Int&, Int&, Int&, Int&,
    C_ele&, Int&, IntArr, Int&);

extern "C" void ele_from_lat_to_f2_(lat_struct*, Int&, C_ele&);
extern "C" void control_from_lat_to_f2_(lat_struct*, Int&, C_control&);


extern "C" void lat_to_f_(C_lat& c, lat_struct* f) {
  const char* name  = c.name.data();            int n_name = c.name.size();
  const char* lat   = c.lattice.data();         int n_lat = c.lattice.size();
  const char* file  = c.input_file_name.data(); int n_file = c.input_file_name.size();
  const char* title = c.title.data();           int n_title = c.title.size();
  int n_con = c.control.size();
  int n_ic  = c.ic.size();
  int n_ele_max = c.ele.size() - 1;
  lat_to_f2_(f, name, n_name, lat, n_lat, file, n_file, title, n_title,
      c.x, c.y, c.z, c.param, 
      c.version, c.n_ele_use, c.n_ele_max, n_ele_max, c.n_control_max, 
      c.n_ic_max, c.input_taylor_order, c.ele_init, n_con, &c.ic[0], n_ic);
  for (int i = 0; i < n_ele_max+1; i++) {
    ele_from_lat_to_f2_(f, i, c.ele[i]);
  }
  for (int i = 0; i < n_con; i++) {
    control_from_lat_to_f2_(f, i+1, c.control[i]);
  }
}

extern "C" void lat_to_c2_(C_lat& c, char* name, char* lat, char* file,
    char* title, mode_info_struct* x, mode_info_struct* y, mode_info_struct* z,
    lat_param_struct* param, Int& ver, Int& n_use, Int& n_max, Int& n_maxx, 
    Int& n_con_max, Int& n_ic_max, Int& n_taylor, 
    ele_struct* ele_init, Int& n_con_array, IntArr ic, Int& n_ic_array) {
  c.name                = name;
  c.lattice             = lat;
  c.input_file_name     = file;
  c.title               = title;
  c.version             = ver;
  c.n_ele_use           = n_use;
  c.n_ele_max           = n_max;
  c.n_control_max       = n_con_max;
  c.n_ic_max            = n_ic_max;
  c.input_taylor_order  = n_taylor;

  if (c.control.size() != n_con_array) c.control.resize(n_con_array);
  if (!(c.ele.size() == n_maxx)) c.ele.resize(n_maxx);
  if (!(c.ic.size() == n_ic_array)) c.ic.resize(n_ic_array);

  c.ic << ic;
  x >> c.x;  y >> c.y;  z >> c.z;
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
  name               = c.name;
  lattice            = c.lattice;
  input_file_name    = c.input_file_name;
  title              = c.title;
  x                  = c.x;
  y                  = c.y;
  z                  = c.z;
  param              = c.param;
  version            = c.version;
  n_ele_use          = c.n_ele_use;
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
