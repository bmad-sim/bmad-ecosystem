#include <iostream>
#include <stdlib.h>
#include "cpp_and_bmad.h"

using namespace std;

//---------------------------------------------------

template <class T> bool is_all_true (const valarray<T>& v1, const valarray<T>& v2) {
  bool is_true = true;
  T t1, t2;
  if (v1.size() != v2.size()) return false;
  for (int i = 0; i < v1.size(); i++) {
    t1 = v1[i];
    t2 = v2[i];
    is_true = is_true && (v1[i] == v2[i]);
  }
  return is_true; 
}

template bool is_all_true(const Real_Array&, const Real_Array&);
template bool is_all_true(const Int_Array&, const Int_Array&);
template bool is_all_true(const C_taylor_term_array&, const C_taylor_term_array&);
template bool is_all_true(const C_sr1_wake_array&, const C_sr1_wake_array&);
template bool is_all_true(const C_sr2_wake_array&, const C_sr2_wake_array&);
template bool is_all_true(const C_lr_wake_array&, const C_lr_wake_array&);
template bool is_all_true(const C_taylor_array&, const C_taylor_array&);
template bool is_all_true(const C_wig_term_array&, const C_wig_term_array&);
template bool is_all_true(const C_control_array&, const C_control_array&);



//---------------------------------------------------

bool is_all_true (const Bool_Array& v) {
  bool is_true = true;
  for (int i = 0; i < v.size(); i++) is_true = is_true && v[i];
  return is_true; 
};

bool is_all_true (const Real_Matrix& v1, const Real_Matrix& v2) {
  bool is_true = true;
  double a1, a2;
  if (v1.size() != v2.size()) return false;
  for (int i = 0; i < v1.size(); i++) {
    if (v1[i].size() != v2[i].size()) return false;
    for (int j = 0; j < v1[i].size(); j++) {
      a1 = v1[i][j]; a2 = v2[i][j];
      is_true = is_true && (v1[i][j] == v2[i][j]);
    }
  }
  return is_true; 
};

//---------------------------------------------------

bool operator== (const C_coord& x, const C_coord& y) {
  return is_all_true(x.vec, y.vec);
};

bool operator== (const C_twiss& x, const C_twiss& y) {
  return (x.beta == y.beta) && (x.alpha == y.alpha) && (x.gamma == y.gamma) &&
     (x.phi == y.phi) && (x.eta == y.eta) && (x.etap == y.etap) &&
     (x.eta_lab == y.eta_lab) && (x.etap_lab == y.etap_lab) && (x.sigma == y.sigma);
};

bool operator== (const C_floor_position& x, const C_floor_position& y) {
 return (x.x == y.x) && (x.y == y.y) && (x.z == y.z) && (x.theta == y.theta) &&
        (x.phi == y.phi) && (x.psi == y.psi);
};

bool operator== (const C_wig_term& x, const C_wig_term& y) {
 return (x.coef == y.coef) && (x.kx == y.kx) && (x.ky == y.ky) && 
        (x.kz == y.kz) && (x.phi_z == y.phi_z) && (x.type == y.type);
};

bool operator== (const C_taylor_term& x, const C_taylor_term& y) {
  bool is_true = (x.coef == y.coef) && (x.exp.size() == y.exp.size());
  if (!is_true) return is_true;
  return is_all_true(x.exp, y.exp);
};

bool operator== (const C_taylor& x, const C_taylor& y) {
  return (x.ref == y.ref) && is_all_true(x.term, y.term);
};

bool operator== (const C_sr1_wake& x, const C_sr1_wake& y) {
  return (x.z == y.z) && (x.longitudinal == y.longitudinal) && (x.transverse == y.transverse);
};

bool operator== (const C_sr2_wake& x, const C_sr2_wake& y) {
  return (x.amp == y.amp) && (x.damp == y.damp) && 
         (x.k == y.k) && (x.phi == y.phi) && 
         (x.norm_sin == y.norm_sin) && (x.norm_cos == y.norm_cos) && 
         (x.skew_sin == y.skew_sin) && (x.skew_cos == y.skew_cos);
};

bool operator== (const C_lr_wake& x, const C_lr_wake& y) {
  return (x.freq == y.freq) && (x.freq_in == y.freq_in) && 
         (x.R_over_Q == y.R_over_Q) && (x.Q == y.Q) && (x.angle == y.angle) &&
         (x.norm_sin == y.norm_sin) && (x.norm_cos == y.norm_cos) && 
         (x.skew_sin == y.skew_sin) && (x.skew_cos == y.skew_cos) && 
         (x.m == y.m) && (x.polarized == y.polarized);
};

bool operator== (const C_wake& x, const C_wake& y) {
  bool all_true;
  all_true = (x.sr_file == y.sr_file) && (x.lr_file == y.lr_file) &&
            (x.sr1.size() == y.sr1.size()) && (x.sr2_long.size() == y.sr2_long.size()) && 
            (x.sr2_trans.size() == y.sr2_trans.size()) && (x.lr.size() == y.lr.size()); 
  if (!all_true) return all_true;
  return is_all_true(x.sr1, y.sr1) && is_all_true(x.sr2_long, y.sr2_long) && 
         is_all_true(x.sr2_trans, y.sr2_trans) && is_all_true(x.lr, y.lr);
}

bool operator== (const C_control& x, const C_control& y) {
  return (x.coef == y.coef) && (x.ix_lord == y.ix_lord) && 
         (x.ix_slave == y.ix_slave) && (x.ix_attrib == y.ix_attrib);
}

bool operator== (const C_param& x, const C_param& y) {
  return (x.n_part == y.n_part) && (x.total_length == y.total_length) && 
         (x.growth_rate == y.growth_rate) &&
         is_all_true(x.t1_with_RF, y.t1_with_RF) && 
         is_all_true(x.t1_no_RF, y.t1_no_RF) && 
         (x.particle == y.particle) && (x.ix_lost == y.ix_lost) && 
         (x.end_lost_at == y.end_lost_at) && 
         (x.lattice_type == y.lattice_type) && 
         (x.ixx == y.ixx) && (x.stable == y.stable) && 
         (x.aperture_limit_on == y.aperture_limit_on) && (x.lost == y.lost);
}

bool operator== (const C_amode& x, const C_amode& y) {
  return (x.emittance == y.emittance) && (x.synch_int4 == y.synch_int4) && 
         (x.synch_int5 == y.synch_int5) && (x.j_damp == y.j_damp) && 
         (x.alpha_damp == y.alpha_damp) && (x.chrom == y.chrom) && 
         (x.tune == y.tune);
};

bool operator== (const C_linac_mode& x, const C_linac_mode& y) {
  return (x.i2_E4 == y.i2_E4) && (x.i3_E7 == y.i3_E7) && 
         (x.i5a_E6 == y.i5a_E6) && (x.i5b_E6 == y.i5b_E6) && 
         (x.sig_E1 == y.sig_E1) &&
         (x.emittance_a == y.emittance_a) && (x.emittance_b == y.emittance_b);
};

bool operator== (const C_modes& x, const C_modes& y) {
  return (x.synch_int1 == y.synch_int1) && (x.synch_int2 == y.synch_int2) && 
         (x.synch_int3 == y.synch_int3) && (x.sigE_E == y.sigE_E) && 
         (x.sig_z == y.sig_z) && (x.e_loss == y.e_loss) && (x.a == y.a) && 
         (x.b == y.b) && (x.z == y.z) && (x.lin == y.lin);
}

bool operator== (const C_bmad_com& x, const C_bmad_com& y) {
  return is_all_true(x.d_orb, y.d_orb) && 
      (x.max_aperture_limit == y.max_aperture_limit) && (x.grad_loss_sr_wake == y.grad_loss_sr_wake) && 
      (x.rel_tollerance == y.rel_tollerance) && 
      (x.abs_tollerance == y.abs_tollerance) && 
      (x.taylor_order == y.taylor_order) && 
      (x.default_integ_order == y.default_integ_order) && 
      (x.default_num_steps == y.default_num_steps) &&  
      (x.canonical_coords == y.canonical_coords) && 
      (x.use_liar_lcavity == y.use_liar_lcavity) && 
      (x.sr_wakes_on == y.sr_wakes_on) &&  (x.lr_wakes_on == y.lr_wakes_on) &&  
      (x.mat6_track_symmetric ==  y.mat6_track_symmetric);
}

bool operator== (const C_em_field& x, const C_em_field& y) {
  return is_all_true(x.E, y.E) && is_all_true(x.B, y.B) && 
    is_all_true(x.kick, y.kick) && is_all_true(x.dE, y.dE) && 
    is_all_true(x.dB, y.dB) && is_all_true(x.dkick, y.dkick) && 
    (x.type == y.type);
}

bool operator== (const C_ele& x, const C_ele& y) {
  return (x.name == y.name) && (x.type == y.type) && (x.alias == y.alias) && 
    (x.attribute_name == y.attribute_name) && (x.x == y.x) && (x.y == y.y) && 
    (x.z == y.z) && (x.floor == y.floor) && is_all_true(x.value, y.value) && 
    is_all_true(x.gen0, y.gen0) && is_all_true(x.vec0, y.vec0) && 
    is_all_true(x.mat6, y.mat6) && is_all_true(x.c_mat, y.c_mat) && 
    (x.gamma_c == y.gamma_c) && (x.s == y.s) && is_all_true(x.r, y.r) && 
    is_all_true(x.a, y.a) && is_all_true(x.b, y.b) && 
    is_all_true(x.const_arr, y.const_arr) && (x.descrip == y.descrip) && 
    is_all_true(x.taylor, y.taylor) && 
    is_all_true(x.wig_term, y.wig_term) && (x.wake == y.wake) && 
    (x.key == y.key) && (x.sub_key == y.sub_key) && 
    (x.control_type == y.control_type) && (x.ix_value == y.ix_value) && 
    (x.n_slave == y.n_slave) && (x.ix1_slave == y.ix1_slave) && 
    (x.ix2_slave == y.ix2_slave) && (x.n_lord == y.n_lord) && 
    (x.ic1_lord == y.ic1_lord) && (x.ic2_lord == y.ic2_lord) && 
    (x.ix_pointer == y.ix_pointer) && (x.ixx == y.ixx) && 
    (x.ix_ele == y.ix_ele) && (x.mat6_calc_method == y.mat6_calc_method) && 
    (x.tracking_method == y.tracking_method) && (x.field_calc == y.field_calc) && 
    (x.num_steps == y.num_steps) && (x.integrator_order == y.integrator_order) && 
    (x.ptc_kind == y.ptc_kind) && (x.taylor_order == y.taylor_order) && 
    (x.aperture_at == y.aperture_at) && (x.symplectify == y.symplectify) && 
    (x.mode_flip == y.mode_flip) && (x.multipoles_on == y.multipoles_on) && 
    (x.exact_rad_int_calc == y.exact_rad_int_calc) && 
    (x.field_master == y.field_master) && 
    (x.is_on == y.is_on) && (x.internal_logic == y.internal_logic) && 
    (x.logic == y.logic) && (x.on_an_i_beam == y.on_an_i_beam);
}

bool operator== (const C_mode_info& x, const C_mode_info& y) {
  return (x.tune == y.tune) && (x.emit == y.emit) && (x.chrom == y.chrom);
}

bool operator== (const C_ring& x, const C_ring& y) {
  bool is_true = true;
  is_true = is_true && (x.name == y.name); 
  is_true = is_true && (x.lattice == y.lattice);
  is_true = is_true && (x.input_file_name == y.input_file_name);
  is_true = is_true && (x.title == y.title);
  is_true = is_true && (x.x == y.x);
  is_true = is_true && (x.y == y.y);
  is_true = is_true && (x.z == y.z);
  is_true = is_true && (x.param == y.param);
  is_true = is_true && (x.version == y.version);
  is_true = is_true && (x.n_ele_use == y.n_ele_use);
  is_true = is_true && (x.n_ele_max == y.n_ele_max);
  is_true = is_true && (x.n_control_max == y.n_control_max);
  is_true = is_true && (x.n_ic_max == y.n_ic_max);
  is_true = is_true && (x.input_taylor_order == y.input_taylor_order);
  is_true = is_true && (x.ele_init == y.ele_init);
  is_true = is_true && is_all_true(x.control, y.control);
  is_true = is_all_true(x.ic, y.ic); 
  for (int i = 0; i < x.n_ele_max+1; i++) {
    is_true = is_true && (x.ele[i] == y.ele[i]);
  }
  return is_true; 
}

//---------------------------------------------------------------------------

void ele_comp (const C_ele& x, const C_ele& y) {

  cout << "name: " << ((x.name == y.name) && (x.type == y.type) && 
      (x.alias == y.alias) && (x.attribute_name == y.attribute_name)) << endl;
  cout << "int:  " << ((x.gamma_c == y.gamma_c) && (x.s == y.s) && 
      (x.descrip == y.descrip) && (x.wake == y.wake) && 
      (x.key == y.key) && (x.sub_key == y.sub_key) && 
      (x.control_type == y.control_type) && (x.ix_value == y.ix_value) && 
      (x.n_slave == y.n_slave) && (x.ix1_slave == y.ix1_slave) && 
      (x.ix2_slave == y.ix2_slave) && (x.n_lord == y.n_lord) && 
      (x.ic1_lord == y.ic1_lord) && (x.ic2_lord == y.ic2_lord) && 
      (x.ix_pointer == y.ix_pointer) && (x.ixx == y.ixx) && 
      (x.ix_ele == y.ix_ele)) << endl;
  cout << "logic: " << ((x.mat6_calc_method == y.mat6_calc_method) && 
      (x.tracking_method == y.tracking_method) && (x.field_calc == y.field_calc) && 
      (x.num_steps == y.num_steps) && (x.integrator_order == y.integrator_order) && 
      (x.ptc_kind == y.ptc_kind) && (x.taylor_order == y.taylor_order) && 
      (x.aperture_at == y.aperture_at) && (x.symplectify == y.symplectify) && 
      (x.mode_flip == y.mode_flip) && (x.multipoles_on == y.multipoles_on) && 
      (x.exact_rad_int_calc == y.exact_rad_int_calc) && 
      (x.field_master == y.field_master) && 
      (x.is_on == y.is_on) && (x.internal_logic == y.internal_logic) && 
      (x.logic == y.logic) && (x.on_an_i_beam == y.on_an_i_beam)) << endl;

  cout << "xyz:    " << ((x.x == y.x) && (x.y == y.y) && (x.z == y.z)) << endl;
  cout << "floor:  " << ((x.floor == y.floor)) << endl;
  cout << "value:  " << is_all_true(x.value, y.value) << endl;
  cout << "gen0:   " << is_all_true(x.gen0, y.gen0) << endl;
  cout << "vec0:   " << is_all_true(x.vec0, y.vec0) << endl;
  cout << "mat6:   " << is_all_true(x.mat6, y.mat6) << endl; 
  cout << "c_mat:  " << is_all_true(x.c_mat, y.c_mat) << endl;
  cout << "a:      " << is_all_true(x.a, y.a) << endl;
  cout << "b:      " << is_all_true(x.b, y.b) << endl;
  cout << "const:  " << is_all_true(x.const_arr, y.const_arr) << endl;
  cout << "taylor: " << is_all_true(x.taylor, y.taylor) << endl;
  cout << "wig:    " << is_all_true(x.wig_term, y.wig_term)  << endl;
  cout << "r:      " << is_all_true(x.r, y.r) << endl;
}
 
