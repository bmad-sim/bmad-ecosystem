#include <iostream>
#include <stdlib.h>
#include "cpp_and_bmad.h"

using namespace std;

//---------------------------------------------------

template <class T> bool is_all_equal (const valarray<T>& v1, const valarray<T>& v2) {
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

template bool is_all_equal(const Complx_Array&, const Complx_Array&);
template bool is_all_equal(const Real_Array&, const Real_Array&);
template bool is_all_equal(const Int_Array&, const Int_Array&);
template bool is_all_equal(const C_taylor_term_array&, const C_taylor_term_array&);
template bool is_all_equal(const C_wig_term_array&, const C_wig_term_array&);
template bool is_all_equal(const C_rf_wake_sr_table_array&, const C_rf_wake_sr_table_array&);
template bool is_all_equal(const C_rf_wake_sr_mode_array&, const C_rf_wake_sr_mode_array&);
template bool is_all_equal(const C_rf_wake_lr_array&, const C_rf_wake_lr_array&);
template bool is_all_equal(const C_control_array&, const C_control_array&);
template bool is_all_equal(const C_branch_array&, const C_branch_array&);
template bool is_all_equal(const C_ele_array&, const C_ele_array&);
template bool is_all_equal(const C_taylor_array&, const C_taylor_array&);
template bool is_all_equal(const C_rf_field_mode_array&, const C_rf_field_mode_array&);
template bool is_all_equal(const C_rf_field_mode_term_array&, const C_rf_field_mode_term_array&);
template bool is_all_equal(const C_wall3d_vertex_array&, const C_wall3d_vertex_array&);
template bool is_all_equal(const C_wall3d_section_array&, const C_wall3d_section_array&);

//---------------------------------------------------

bool is_all_equal (const Bool_Array& v) {
  bool is_true = true;
  for (int i = 0; i < v.size(); i++) is_true = is_true && v[i];
  return is_true; 
};

bool is_all_equal (const Real_Matrix& v1, const Real_Matrix& v2) {
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
  return is_all_equal(x.vec, y.vec) && (x.s == y.s) && (x.t == y.t) && (x.spin1 == y.spin1) && (x.spin2 == y.spin2) && 
    (x.e_field_x == y.e_field_x) && (x.e_field_y == y.e_field_y) && (x.phase_x == y.phase_x) && (x.phase_y == y.phase_y);
};

bool operator== (const C_twiss& x, const C_twiss& y) {
  return (x.beta == y.beta) && (x.alpha == y.alpha) && (x.gamma == y.gamma) &&
     (x.phi == y.phi) && (x.eta == y.eta) && (x.etap == y.etap) &&
     (x.sigma == y.sigma) && (x.sigma_p == y.sigma_p) && (x.emit == y.emit) &&
     (x.norm_emit == y.norm_emit);
};

bool operator== (const C_xy_disp& x, const C_xy_disp& y) {
  return (x.eta == y.eta) && (x.etap == y.etap);
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
  bool is_true = (x.coef == y.coef) && (x.expn.size() == y.expn.size());
  if (!is_true) return is_true;
  return is_all_equal(x.expn, y.expn);
};

bool operator== (const C_taylor& x, const C_taylor& y) {
  return (x.ref == y.ref) && is_all_equal(x.term, y.term);
};

bool operator== (const C_rf_wake_sr_table& x, const C_rf_wake_sr_table& y) {
  return (x.z == y.z) && (x.longitudinal == y.longitudinal) && (x.transverse == y.transverse);
};

bool operator== (const C_rf_wake_sr_mode& x, const C_rf_wake_sr_mode& y) {
  return (x.amp == y.amp) && (x.damp == y.damp) && 
         (x.k == y.k) && (x.phi == y.phi) && 
         (x.b_sin == y.b_sin) && (x.b_cos == y.b_cos) && 
         (x.a_sin == y.a_sin) && (x.a_cos == y.a_cos);
};

bool operator== (const C_rf_wake_lr& x, const C_rf_wake_lr& y) {
  return (x.freq == y.freq) && (x.freq_in == y.freq_in) && 
         (x.R_over_Q == y.R_over_Q) && (x.Q == y.Q) && (x.angle == y.angle) &&
         (x.b_sin == y.b_sin) && (x.b_cos == y.b_cos) && 
         (x.a_sin == y.a_sin) && (x.a_cos == y.a_cos) && 
         (x.t_ref == y.t_ref) && (x.m == y.m) && (x.polarized == y.polarized);
};

bool operator== (const C_rf_wake& x, const C_rf_wake& y) {
  bool all_true = true;
  all_true = all_true && (x.sr_file == y.sr_file);
  all_true = all_true && (x.lr_file == y.lr_file);
  all_true = all_true && (x.sr_table.size() == y.sr_table.size());
  all_true = all_true && (x.sr_mode_long.size() == y.sr_mode_long.size());
  all_true = all_true && (x.sr_mode_trans.size() == y.sr_mode_trans.size());
  all_true = all_true && (x.lr.size() == y.lr.size());;
  all_true = all_true && is_all_equal(x.sr_table, y.sr_table);
  all_true = all_true && is_all_equal(x.sr_mode_long, y.sr_mode_long);
  all_true = all_true && is_all_equal(x.sr_mode_trans, y.sr_mode_trans);
  all_true = all_true && is_all_equal(x.lr, y.lr);
  return all_true;
}

bool operator== (const C_rf_field_mode_term& x, const C_rf_field_mode_term& y) {
  return (x.e == y.e) && (x.b == y.b);
}

bool operator== (const C_rf_field_mode& x, const C_rf_field_mode& y) {
  return (x.m == y.m) && (x.freq == y.freq) && (x.f_damp == y.f_damp) && (x.theta_t0 == y.theta_t0) && 
        (x.stored_energy == y.stored_energy) && (x.phi_0 == y.phi_0) && (x.dz == y.dz) && is_all_equal(x.term, y.term);
}

bool operator== (const C_rf_field& x, const C_rf_field& y) {
  return is_all_equal(x.mode, y.mode);
}

bool operator== (const C_rf& x, const C_rf& y) {
  return (x.wake == y.wake) && (x.field == y.field);
}

bool operator== (const C_wall3d_vertex& x, const C_wall3d_vertex& y) {
  return (x.x == y.x) && (x.y == y.y) && (x.radius_x == y.radius_x) && (x.radius_y == y.radius_y) && 
         (x.tilt == y.tilt) && (x.angle == y.angle) && (x.x0 == y.x0) && (x.y0 == y.y0);
}

bool operator== (const C_wall3d_section& x, const C_wall3d_section& y) {
  return (x.type == y.type) && (x.s == y.s) && is_all_equal(x.s_spline, y.s_spline) && 
         (x.n_slice_spline == y.n_slice_spline) && is_all_equal (x.v, y.v) && 
         (x.n_vertex_input == y.n_vertex_input);
}

bool operator== (const C_wall3d& x, const C_wall3d& y) {
  return is_all_equal (x.section, y.section);
}

bool operator== (const C_control& x, const C_control& y) {
  return (x.coef == y.coef) && (x.ix_lord == y.ix_lord) && 
         (x.ix_slave == y.ix_slave) && (x.ix_attrib == y.ix_attrib);
}

bool operator== (const C_lat_param& x, const C_lat_param& y) {
  return (x.n_part == y.n_part) && (x.total_length == y.total_length) && 
         (x.unstable_factor == y.unstable_factor) &&
         is_all_equal(x.t1_with_RF, y.t1_with_RF) && 
         is_all_equal(x.t1_no_RF, y.t1_no_RF) && 
         (x.particle == y.particle) && (x.ix_lost == y.ix_lost) && 
         (x.end_lost_at == y.end_lost_at) && 
         (x.lattice_type == y.lattice_type) && 
         (x.ixx == y.ixx) && (x.stable == y.stable) && 
         (x.aperture_limit_on == y.aperture_limit_on) && (x.lost == y.lost);
}

bool operator== (const C_anormal_mode& x, const C_anormal_mode& y) {
  return (x.emittance == y.emittance) && (x.synch_int4 == y.synch_int4) && 
         (x.synch_int5 == y.synch_int5) && (x.j_damp == y.j_damp) && 
         (x.alpha_damp == y.alpha_damp) && (x.chrom == y.chrom) && 
         (x.tune == y.tune);
};

bool operator== (const C_linac_normal_mode& x, const C_linac_normal_mode& y) {
  return (x.i2_E4 == y.i2_E4) && (x.i3_E7 == y.i3_E7) && 
         (x.i5a_E6 == y.i5a_E6) && (x.i5b_E6 == y.i5b_E6) && 
         (x.sig_E1 == y.sig_E1) && (x.a_emittance_end == y.a_emittance_end) && 
         (x.b_emittance_end == y.b_emittance_end);
};

bool operator== (const C_normal_modes& x, const C_normal_modes& y) {
  return is_all_equal(x.synch_int == y.synch_int) && (x.sigE_E == y.sigE_E) && 
         (x.sig_z == y.sig_z) && (x.e_loss == y.e_loss) && 
         (x.pz_aperture == y.pz_aperture) && (x.a == y.a) && 
         (x.b == y.b) && (x.z == y.z) && (x.lin == y.lin);
}

bool operator== (const C_branch& x, const C_branch& y) {
  return (x.name == y.name) && (x.ix_branch == y.ix_branch) && 
         (x.ix_from_branch == y.ix_from_branch) && (x.ix_from_ele == y.ix_from_ele) && 
         (x.n_ele_track == y.n_ele_track) && (x.n_ele_max == y.n_ele_max) && 
         is_all_equal(x.ele == y.ele) && (x.wall3d == y.wall3d) && (x.param == y.param);
}

bool operator== (const C_bmad_com& x, const C_bmad_com& y) {
  return is_all_equal(x.d_orb, y.d_orb) && 
      (x.max_aperture_limit == y.max_aperture_limit) && 
      (x.grad_loss_sr_wake == y.grad_loss_sr_wake) && 
      (x.rel_tolerance == y.rel_tolerance) && 
      (x.abs_tolerance == y.abs_tolerance) && 
      (x.rel_tol_adaptive_tracking == y.rel_tol_adaptive_tracking) && 
      (x.abs_tol_adaptive_tracking == y.abs_tol_adaptive_tracking) && 
      (x.taylor_order == y.taylor_order) && 
      (x.default_integ_order == y.default_integ_order) && 
      (x.default_ds_step == y.default_ds_step) &&  
      (x.canonical_coords == y.canonical_coords) && 
      (x.significant_longitudinal_length == y.significant_longitudinal_length) && 
      (x.sr_wakes_on == y.sr_wakes_on) &&  (x.lr_wakes_on == y.lr_wakes_on) &&  
      (x.mat6_track_symmetric ==  y.mat6_track_symmetric) &&
      (x.auto_bookkeeper == y.auto_bookkeeper) &&
      (x.trans_space_charge_on == y.trans_space_charge_on) &&
      (x.coherent_synch_rad_on == y.coherent_synch_rad_on) &&
      (x.spin_tracking_on == y.spin_tracking_on) &&
      (x.radiation_damping_on == y.radiation_damping_on) &&
      (x.radiation_fluctuations_on == y.radiation_fluctuations_on) &&
      (x.conserve_taylor_maps == y.conserve_taylor_maps);
}

bool operator== (const C_em_field& x, const C_em_field& y) {
  return is_all_equal(x.E, y.E) && is_all_equal(x.B, y.B) && 
    is_all_equal(x.dE, y.dE) && is_all_equal(x.dB, y.dB);
}

bool operator== (const C_ele& x, const C_ele& y) {
  bool is_true(true);
  is_true = is_true && (x.name == y.name);
  is_true = is_true && (x.type == y.type);
  is_true = is_true && (x.alias == y.alias);
  is_true = is_true && (x.component_name == y.component_name);
  is_true = is_true && (x.descrip == y.descrip);
  is_true = is_true && (x.a == y.a);
  is_true = is_true && (x.b == y.b);
  is_true = is_true && (x.z == y.z);
  is_true = is_true && (x.x == y.x);
  is_true = is_true && (x.y == y.y);
  is_true = is_true && (x.floor == y.floor);
  // is_true = is_true && (x.mode3 == y.mode3);
  is_true = is_true && (x.map_ref_orb_in == y.map_ref_orb_in);
  is_true = is_true && (x.map_ref_orb_out == y.map_ref_orb_out);
  // is_true = is_true && (x.gen_field == y.gen_field);
  is_true = is_true && is_all_equal(x.taylor, y.taylor);
  // is_true = is_true && (x.rf == y.rf);
  // is_true = is_true && is_all_equal(x.wig_term, y.wig_term);
  // is_true = is_true && (x.space_charge == y.space_charge);
  // is_true = is_true && (x.wall3d == y.wall3d);
  is_true = is_true && is_all_equal(x.value, y.value);
  is_true = is_true && is_all_equal(x.old_value, y.old_value);
  is_true = is_true && is_all_equal(x.gen0, y.gen0);
  is_true = is_true && is_all_equal(x.vec0, y.vec0);
  is_true = is_true && is_all_equal(x.mat6, y.mat6);
  is_true = is_true && is_all_equal(x.c_mat, y.c_mat);
  is_true = is_true && (x.gamma_c == y.gamma_c);
  is_true = is_true && (x.s == y.s);
  is_true = is_true && (x.ref_time == y.ref_time);
  is_true = is_true && is_all_equal(x.r, y.r);
  is_true = is_true && is_all_equal(x.a_pole, y.a_pole);
  is_true = is_true && is_all_equal(x.b_pole, y.b_pole);
  is_true = is_true && is_all_equal(x.const_arr, y.const_arr);
  is_true = is_true && (x.key == y.key);
  is_true = is_true && (x.sub_key == y.sub_key);
  is_true = is_true && (x.ix_ele == y.ix_ele);
  is_true = is_true && (x.ix_branch == y.ix_branch);
  is_true = is_true && (x.ix_value == y.ix_value);
  is_true = is_true && (x.slave_status == y.slave_status);
  is_true = is_true && (x.n_slave == y.n_slave);
  is_true = is_true && (x.ix1_slave == y.ix1_slave);
  is_true = is_true && (x.ix2_slave == y.ix2_slave);
  is_true = is_true && (x.lord_status == y.lord_status);
  is_true = is_true && (x.n_lord == y.n_lord);
  is_true = is_true && (x.ic1_lord == y.ic1_lord);
  is_true = is_true && (x.ic2_lord == y.ic2_lord);
  is_true = is_true && (x.ix_pointer == y.ix_pointer);
  is_true = is_true && (x.ixx == y.ixx);
  is_true = is_true && (x.mat6_calc_method == y.mat6_calc_method);
  is_true = is_true && (x.tracking_method == y.tracking_method);
  is_true = is_true && (x.field_calc == y.field_calc);
  is_true = is_true && (x.ref_orbit == y.ref_orbit);
  is_true = is_true && (x.aperture_at == y.aperture_at);
  is_true = is_true && (x.aperture_type == y.aperture_type);
  is_true = is_true && (x.symplectify == y.symplectify);
  is_true = is_true && (x.mode_flip == y.mode_flip);
  is_true = is_true && (x.multipoles_on == y.multipoles_on);
  is_true = is_true && (x.scale_multipoles == y.scale_multipoles);
  is_true = is_true && (x.map_with_offsets == y.map_with_offsets);
  is_true = is_true && (x.reversed == y.reversed);
  is_true = is_true && (x.is_on == y.is_on);
  is_true = is_true && (x.old_is_on == y.old_is_on);
  is_true = is_true && (x.logic == y.logic);
  is_true = is_true && (x.bmad_logic == y.bmad_logic);
  is_true = is_true && (x.on_a_girder == y.on_a_girder);
  is_true = is_true && (x.csr_calc_on == y.csr_calc_on);
  is_true = is_true && (x.offset_moves_aperture == y.offset_moves_aperture);
  return is_true;
}

bool operator== (const C_mode_info& x, const C_mode_info& y) {
  return (x.tune == y.tune) && (x.emit == y.emit) && (x.chrom == y.chrom);
}

bool operator== (const C_lat& x, const C_lat& y) {
  bool is_true = true;
  is_true = is_true && (x.use_name == y.use_name); 
  is_true = is_true && (x.lattice == y.lattice);
  is_true = is_true && (x.input_file_name == y.input_file_name);
  is_true = is_true && (x.title == y.title);
  is_true = is_true && (x.a == y.a);
  is_true = is_true && (x.b == y.b);
  is_true = is_true && (x.z == y.z);
  is_true = is_true && (x.param == y.param);
  is_true = is_true && (x.version == y.version);
  is_true = is_true && (x.n_ele_track == y.n_ele_track);
  is_true = is_true && (x.n_ele_max == y.n_ele_max);
  is_true = is_true && (x.n_control_max == y.n_control_max);
  is_true = is_true && (x.n_ic_max == y.n_ic_max);
  is_true = is_true && (x.input_taylor_order == y.input_taylor_order);
  is_true = is_true && (x.ele_init == y.ele_init);
  is_true = is_true && is_all_equal(x.control, y.control);
  is_true = is_all_equal(x.ic, y.ic); 
  for (int i = 0; i < x.n_ele_max+1; i++) {
    is_true = is_true && (x.ele[i] == y.ele[i]);
  }
  return is_true; 
}

 
