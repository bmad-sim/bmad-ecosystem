
program cpp_bmad_interface_test

use bmad_cpp_test_mod

logical ok, all_ok

!

all_ok = .true.
call test1_f_spline(ok); if (.not. ok) all_ok = .false.
call test1_f_spin_polar(ok); if (.not. ok) all_ok = .false.
call test1_f_surface_orientation(ok); if (.not. ok) all_ok = .false.
call test1_f_ac_kicker_time(ok); if (.not. ok) all_ok = .false.
call test1_f_ac_kicker_freq(ok); if (.not. ok) all_ok = .false.
call test1_f_ac_kicker(ok); if (.not. ok) all_ok = .false.
call test1_f_interval1_coef(ok); if (.not. ok) all_ok = .false.
call test1_f_photon_reflect_table(ok); if (.not. ok) all_ok = .false.
call test1_f_photon_reflect_surface(ok); if (.not. ok) all_ok = .false.
call test1_f_coord(ok); if (.not. ok) all_ok = .false.
call test1_f_coord_array(ok); if (.not. ok) all_ok = .false.
call test1_f_bpm_phase_coupling(ok); if (.not. ok) all_ok = .false.
call test1_f_expression_atom(ok); if (.not. ok) all_ok = .false.
call test1_f_wake_sr_mode(ok); if (.not. ok) all_ok = .false.
call test1_f_wake_sr(ok); if (.not. ok) all_ok = .false.
call test1_f_wake_lr_mode(ok); if (.not. ok) all_ok = .false.
call test1_f_wake_lr(ok); if (.not. ok) all_ok = .false.
call test1_f_lat_ele_loc(ok); if (.not. ok) all_ok = .false.
call test1_f_wake(ok); if (.not. ok) all_ok = .false.
call test1_f_taylor_term(ok); if (.not. ok) all_ok = .false.
call test1_f_taylor(ok); if (.not. ok) all_ok = .false.
call test1_f_em_taylor_term(ok); if (.not. ok) all_ok = .false.
call test1_f_em_taylor(ok); if (.not. ok) all_ok = .false.
call test1_f_cartesian_map_term1(ok); if (.not. ok) all_ok = .false.
call test1_f_cartesian_map_term(ok); if (.not. ok) all_ok = .false.
call test1_f_cartesian_map(ok); if (.not. ok) all_ok = .false.
call test1_f_cylindrical_map_term1(ok); if (.not. ok) all_ok = .false.
call test1_f_cylindrical_map_term(ok); if (.not. ok) all_ok = .false.
call test1_f_cylindrical_map(ok); if (.not. ok) all_ok = .false.
call test1_f_grid_field_pt1(ok); if (.not. ok) all_ok = .false.
call test1_f_grid_field_pt(ok); if (.not. ok) all_ok = .false.
call test1_f_grid_field(ok); if (.not. ok) all_ok = .false.
call test1_f_floor_position(ok); if (.not. ok) all_ok = .false.
call test1_f_high_energy_space_charge(ok); if (.not. ok) all_ok = .false.
call test1_f_xy_disp(ok); if (.not. ok) all_ok = .false.
call test1_f_twiss(ok); if (.not. ok) all_ok = .false.
call test1_f_mode3(ok); if (.not. ok) all_ok = .false.
call test1_f_bookkeeping_state(ok); if (.not. ok) all_ok = .false.
call test1_f_rad_map(ok); if (.not. ok) all_ok = .false.
call test1_f_rad_map_ele(ok); if (.not. ok) all_ok = .false.
call test1_f_gen_grad1(ok); if (.not. ok) all_ok = .false.
call test1_f_gen_grad_map(ok); if (.not. ok) all_ok = .false.
call test1_f_surface_grid_pt(ok); if (.not. ok) all_ok = .false.
call test1_f_surface_grid(ok); if (.not. ok) all_ok = .false.
call test1_f_target_point(ok); if (.not. ok) all_ok = .false.
call test1_f_surface_curvature(ok); if (.not. ok) all_ok = .false.
call test1_f_photon_target(ok); if (.not. ok) all_ok = .false.
call test1_f_photon_material(ok); if (.not. ok) all_ok = .false.
call test1_f_pixel_pt(ok); if (.not. ok) all_ok = .false.
call test1_f_pixel_detec(ok); if (.not. ok) all_ok = .false.
call test1_f_photon_element(ok); if (.not. ok) all_ok = .false.
call test1_f_wall3d_vertex(ok); if (.not. ok) all_ok = .false.
call test1_f_wall3d_section(ok); if (.not. ok) all_ok = .false.
call test1_f_wall3d(ok); if (.not. ok) all_ok = .false.
call test1_f_ramper_lord(ok); if (.not. ok) all_ok = .false.
call test1_f_control(ok); if (.not. ok) all_ok = .false.
call test1_f_control_var1(ok); if (.not. ok) all_ok = .false.
call test1_f_control_ramp1(ok); if (.not. ok) all_ok = .false.
call test1_f_controller(ok); if (.not. ok) all_ok = .false.
call test1_f_ellipse_beam_init(ok); if (.not. ok) all_ok = .false.
call test1_f_kv_beam_init(ok); if (.not. ok) all_ok = .false.
call test1_f_grid_beam_init(ok); if (.not. ok) all_ok = .false.
call test1_f_beam_init(ok); if (.not. ok) all_ok = .false.
call test1_f_lat_param(ok); if (.not. ok) all_ok = .false.
call test1_f_mode_info(ok); if (.not. ok) all_ok = .false.
call test1_f_pre_tracker(ok); if (.not. ok) all_ok = .false.
call test1_f_anormal_mode(ok); if (.not. ok) all_ok = .false.
call test1_f_linac_normal_mode(ok); if (.not. ok) all_ok = .false.
call test1_f_normal_modes(ok); if (.not. ok) all_ok = .false.
call test1_f_em_field(ok); if (.not. ok) all_ok = .false.
call test1_f_strong_beam(ok); if (.not. ok) all_ok = .false.
call test1_f_track_point(ok); if (.not. ok) all_ok = .false.
call test1_f_track(ok); if (.not. ok) all_ok = .false.
call test1_f_space_charge_common(ok); if (.not. ok) all_ok = .false.
call test1_f_bmad_common(ok); if (.not. ok) all_ok = .false.
call test1_f_rad_int1(ok); if (.not. ok) all_ok = .false.
call test1_f_rad_int_branch(ok); if (.not. ok) all_ok = .false.
call test1_f_rad_int_all_ele(ok); if (.not. ok) all_ok = .false.
call test1_f_ele(ok); if (.not. ok) all_ok = .false.
call test1_f_complex_taylor_term(ok); if (.not. ok) all_ok = .false.
call test1_f_complex_taylor(ok); if (.not. ok) all_ok = .false.
call test1_f_branch(ok); if (.not. ok) all_ok = .false.
call test1_f_lat(ok); if (.not. ok) all_ok = .false.
call test1_f_bunch(ok); if (.not. ok) all_ok = .false.
call test1_f_bunch_params(ok); if (.not. ok) all_ok = .false.
call test1_f_beam(ok); if (.not. ok) all_ok = .false.
call test1_f_aperture_point(ok); if (.not. ok) all_ok = .false.
call test1_f_aperture_param(ok); if (.not. ok) all_ok = .false.
call test1_f_aperture_scan(ok); if (.not. ok) all_ok = .false.

print *
if (all_ok) then
  print *, 'Bottom Line: Everything OK!'
else
  print *, 'BOTTOM LINE: PROBLEMS FOUND!'
endif

end program
