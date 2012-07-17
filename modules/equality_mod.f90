module equality_mod

use bmad_struct

interface operator (==)
  module procedure eq_coord, eq_twiss, eq_xy_disp, eq_floor_position
  module procedure eq_taylor_term, eq_taylor, eq_wig_term, eq_mode3
  module procedure eq_sr_table_wake, eq_sr_mode_wake, eq_lr_wake, eq_branch
  module procedure eq_rf_wake, eq_em_field_mode, eq_em_fields, eq_space_charge
  module procedure eq_modes, eq_bmad_com, eq_em_field, eq_ele, eq_mode_info
  module procedure eq_lat, eq_control, eq_param, eq_amode, eq_linac_mode
  module procedure eq_wall3d_vertex, eq_wall3d_section, eq_wall3d
end interface

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_coord (f1, f2) result (is_eq)

implicit none

type (coord_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (all(f1%vec == f2%vec) .and. all(f1%spin == f2%spin)) .and. &
        (f1%s == f2%s) .and. (f1%t == f2%t) .and. &
        (f1%e_field_x == f2%e_field_x) .and. (f1%e_field_y == f2%e_field_y) .and. &
        (f1%phase_x == f2%phase_x) .and. (f1%phase_y == f2%phase_y) .and. & 
        (f1%charge == f2%charge) .and. (f1%p0c == f2%p0c) .and. &
        (f1%beta == f2%beta) .and. (f1%species == f2%species) .and. &
        (f1%ix_ele == f2%ix_ele) .and. (f1%state == f2%state)
        

end function eq_coord

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_twiss (f1, f2) result (is_eq)

implicit none

type (twiss_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%beta == f2%beta) .and. (f1%alpha == f2%alpha) .and. &
          (f1%gamma == f2%gamma) .and. (f1%phi == f2%phi) .and. &
          (f1%eta == f2%eta) .and. (f1%etap == f2%etap) .and. &
          (f1%sigma == f2%sigma) .and. (f1%sigma_p == f2%sigma_p) .and. &
          (f1%emit == f2%emit) .and. (f1%norm_emit == f2%norm_emit)

end function eq_twiss

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_xy_disp (f1, f2) result (is_eq)

implicit none

type (xy_disp_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%eta == f2%eta) .and. (f1%etap == f2%etap)

end function eq_xy_disp

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_floor_position (f1, f2) result (is_eq)

implicit none

type (floor_position_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%x == f2%x) .and. (f1%y == f2%y) .and. &
          (f1%z == f2%z) .and. (f1%theta == f2%theta) .and. &
          (f1%phi == f2%phi) .and. (f1%psi == f2%psi)

end function eq_floor_position

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_wig_term (f1, f2) result (is_eq)

implicit none

type (wig_term_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%coef == f2%coef) .and. (f1%kx == f2%kx) .and. &
          (f1%ky == f2%ky) .and. (f1%kz == f2%kz) .and. &
          (f1%phi_z == f2%phi_z) .and. (f1%type == f2%type)

end function eq_wig_term

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_taylor_term (f1, f2) result (is_eq)

implicit none

type (taylor_term_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%coef == f2%coef) .and. all(f1%expn == f2%expn)

end function eq_taylor_term

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_taylor (f1, f2) result (is_eq)

implicit none

type (taylor_struct), intent(in) :: f1, f2
logical is_eq
integer i

!

is_eq = (f1%ref == f2%ref) .and. (associated(f1%term) .eqv. associated(f1%term))

if (associated(f1%term)) then
  is_eq = is_eq .and. (size(f1%term) == size(f2%term))
  if (.not. is_eq) return
  do i = 1, size(f1%term)
    is_eq = is_eq .and. (f1%term(i) == f2%term(i))
  enddo
endif

end function eq_taylor

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_sr_table_wake (f1, f2) result (is_eq)

implicit none

type (rf_wake_sr_table_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%z == f2%z) .and. (f1%long == f2%long) .and. (f1%trans == f2%trans)

end function eq_sr_table_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_sr_mode_wake (f1, f2) result (is_eq)

implicit none

type (rf_wake_sr_mode_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%amp == f2%amp) .and. (f1%damp == f2%damp) .and. &
        (f1%k == f2%k) .and. (f1%phi == f2%phi) .and. &
        (f1%b_sin == f2%b_sin) .and. (f1%b_cos == f2%b_cos) .and. &
        (f1%a_sin == f2%a_sin) .and. (f1%a_cos == f2%a_cos)

end function eq_sr_mode_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_lr_wake (f1, f2) result (is_eq)

implicit none

type (rf_wake_lr_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%freq == f2%freq) .and. (f1%r_over_q == f2%r_over_q) .and. &
        (f1%freq_in == f2%freq_in) .and. (f1%Q == f2%Q) .and. &
        (f1%m == f2%m) .and. (f1%b_sin == f2%b_sin) .and. &
        (f1%b_cos == f2%b_cos) .and. (f1%a_sin == f2%a_sin) .and. &
        (f1%a_cos == f2%a_cos) .and. (f1%t_ref == f2%t_ref) .and. &
        (f1%m == f2%m) .and. (f1%angle == f2%angle) .and. &
        (f1%polarized .eqv. f2%polarized)

end function eq_lr_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_em_field_mode (f1, f2) result (is_eq)

implicit none

type (em_field_mode_struct), intent(in) :: f1, f2
integer i
logical is_eq

!

is_eq = (f1%harmonic == f2%harmonic) .and. (f1%f_damp == f2%f_damp) .and. (f1%dphi0_ref == f2%dphi0_ref) .and. &
        (f1%stored_energy == f2%stored_energy) .and. (f1%m == f2%m) .and. (f1%phi0_azimuth == f2%phi0_azimuth) .and. &
        (f1%field_scale == f2%field_scale)
if (.not. is_eq) return

is_eq  = (associated(f1%map) .eqv. associated(f2%map))
is_eq  = (associated(f1%grid) .eqv. associated(f2%grid))
if (.not. is_eq) return

if (associated(f1%map)) is_eq = is_eq .and. (size(f1%map%term) == size(f2%map%term))
if (associated(f1%grid)) is_eq = is_eq .and. (size(f1%grid%pt) == size(f2%grid%pt))
if (.not. is_eq) return

! Notice that file names do not have to be the same

if (associated(f1%map)) then
  is_eq = (f1%map%dz == f2%map%dz) .and. all(f1%map%term%e_coef == f2%map%term%e_coef) .and. all(f1%map%term%b_coef == f2%map%term%b_coef) 
endif

if (associated(f1%grid)) then
  is_eq = is_eq .and. (f1%grid%type == f2%grid%type) .and. all(f1%grid%dr == f2%grid%dr)
  is_eq = is_eq .and. all(f1%grid%r0 == f2%grid%r0)
  if (.not. is_eq) return
  is_eq = is_eq .and. all(f1%grid%pt%e(1) == f2%grid%pt%e(1)) .and. all(f1%grid%pt%b(1) == f2%grid%pt%b(1)) 
  is_eq = is_eq .and. all(f1%grid%pt%e(2) == f2%grid%pt%e(2)) .and. all(f1%grid%pt%b(2) == f2%grid%pt%b(2)) 
  is_eq = is_eq .and. all(f1%grid%pt%e(3) == f2%grid%pt%e(3)) .and. all(f1%grid%pt%b(3) == f2%grid%pt%b(3)) 
  if (.not. is_eq) return
endif

end function eq_em_field_mode

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_em_fields (f1, f2) result (is_eq)

implicit none

type (em_fields_struct), intent(in) :: f1, f2
integer i
logical is_eq

!

is_eq  = (allocated(f1%mode) .eqv. allocated(f2%mode)) 
if (.not. is_eq) return

if (allocated(f1%mode)) then
  is_eq = (size(f1%mode) == size(f2%mode))
  if (.not. is_eq) return

  do i = 1, size(f1%mode)
    is_eq = (f1%mode(i) == f2%mode(i))
    if (.not. is_eq) return
  enddo
endif

end function eq_em_fields

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_rf_wake (f1, f2) result (is_eq)

implicit none

type (rf_wake_struct), intent(in) :: f1, f2
integer i
logical is_eq

!

is_eq = (f1%sr_file == f2%sr_file) .and. (f1%lr_file == f2%lr_file) .and. &
     (size(f1%sr_table) == size(f2%sr_table)) .and. (size(f2%sr_mode_long) == size(f2%sr_mode_long)) .and. &
     (size(f1%sr_mode_trans) == size(f2%sr_mode_trans)) .and. (size(f1%lr) == size(f2%lr)) .and. &
     (f1%z_sr_mode_max == f2%z_sr_mode_max)
if (.not. is_eq) return

do i = lbound(f1%sr_table, 1), ubound(f1%sr_table, 1)
  is_eq = is_eq .and. (f1%sr_table(i) == f2%sr_table(i)) 
enddo

do i = lbound(f1%sr_mode_long, 1), ubound(f1%sr_mode_long, 1)
  is_eq = is_eq .and. (f1%sr_mode_long(i) == f2%sr_mode_long(i)) 
enddo

do i = lbound(f1%sr_mode_trans, 1), ubound(f1%sr_mode_trans, 1)
  is_eq = is_eq .and. (f1%sr_mode_trans(i) == f2%sr_mode_trans(i)) 
enddo

do i = lbound(f1%lr, 1), ubound(f1%lr, 1)
  is_eq = is_eq .and. (f1%lr(i) == f2%lr(i)) 
enddo

end function eq_rf_wake

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_control (f1, f2) result (is_eq)

implicit none

type (control_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%coef == f2%coef) .and. (f1%ix_lord == f2%ix_lord) .and. &
        (f1%ix_slave == f2%ix_slave) .and. (f1%ix_attrib == f2%ix_attrib) .and. &
        (f1%ix_branch == f2%ix_branch) 

end function eq_control

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_param (f1, f2) result (is_eq)

implicit none

type (lat_param_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%n_part == f2%n_part) .and. (f1%total_length == f2%total_length) .and. &
     (f1%unstable_factor == f2%unstable_factor) .and. &
     all(f1%t1_with_RF == f2%t1_with_RF) .and. &
     all(f1%t1_no_RF == f2%t1_no_RF) .and. &
     (f1%particle == f2%particle) .and. &
     (f1%lattice_type == f2%lattice_type) .and. & 
     (f1%ixx == f2%ixx) .and. (f1%stable .eqv. f2%stable) .and. &
     (f1%aperture_limit_on .eqv. f2%aperture_limit_on) 

end function eq_param

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_amode (f1, f2) result (is_eq)

implicit none

type (anormal_mode_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%emittance == f2%emittance) .and. &
         all(f1%synch_int == f2%synch_int) .and. &
         (f1%j_damp == f2%j_damp) .and. (f1%alpha_damp == f2%alpha_damp) .and. &
         (f1%chrom == f2%chrom) .and. (f1%tune == f2%tune)

end function eq_amode

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_linac_mode (f1, f2) result (is_eq)

implicit none

type (linac_normal_mode_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%i2_E4 == f2%i2_E4) .and. (f1%i3_E7 == f2%i3_E7) .and. &
         (f1%i5a_E6 == f2%i5a_E6) .and. (f1%i5b_E6 == f2%i5b_E6) .and. &
         (f1%a_emittance_end == f2%a_emittance_end) .and. &
         (f1%b_emittance_end == f2%b_emittance_end)

end function eq_linac_mode

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_modes (f1, f2) result (is_eq)

implicit none

type (normal_modes_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = all(f1%synch_int == f2%synch_int) .and. (f1%sige_e == f2%sige_e) .and. &
 (f1%sig_z == f2%sig_z) .and. (f1%e_loss == f2%e_loss) .and. &
 (f1%rf_voltage == f2%rf_voltage) .and. (f1%pz_aperture == f2%pz_aperture) .and. &
 (f1%a == f2%a) .and. &
 (f1%b == f2%b) .and. (f1%z == f2%z) .and. (f1%lin == f2%lin)

end function eq_modes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_bmad_com (f1, f2) result (is_eq)

implicit none

type (bmad_common_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = all(f1%d_orb == f2%d_orb) .and. &
      (f1%max_aperture_limit == f2%max_aperture_limit) .and. &
      (f1%rel_tol_tracking == f2%rel_tol_tracking) .and. &
      (f1%abs_tol_tracking == f2%abs_tol_tracking) .and. &
      (f1%rel_tol_adaptive_tracking == f2%rel_tol_adaptive_tracking) .and. &
      (f1%abs_tol_adaptive_tracking == f2%abs_tol_adaptive_tracking) .and. &
      (f1%taylor_order == f2%taylor_order) .and. &
      (f1%default_integ_order == f2%default_integ_order) .and. &
      (f1%default_ds_step == f2%default_ds_step) .and.  &
      (f1%canonical_coords .eqv. f2%canonical_coords) .and. &
      (f1%significant_length  == f2%significant_length) .and. &
      (f1%sr_wakes_on .eqv. f2%sr_wakes_on) .and.  &
      (f1%lr_wakes_on .eqv. f2%lr_wakes_on) .and.  &
      (f1%mat6_track_symmetric .eqv.  f2%mat6_track_symmetric) .and. &
      (f1%auto_bookkeeper .eqv. f2%auto_bookkeeper) .and. &
      (f1%spin_tracking_on .eqv. f2%spin_tracking_on) .and. &
      (f1%space_charge_on .eqv. f2%space_charge_on) .and. &
      (f1%coherent_synch_rad_on .eqv. f2%coherent_synch_rad_on) .and. &
      (f1%radiation_damping_on .eqv. f2%radiation_damping_on) .and. &
      (f1%radiation_fluctuations_on .eqv. f2%radiation_fluctuations_on) .and. &
      (f1%conserve_taylor_maps .eqv. f2%conserve_taylor_maps)

end function eq_bmad_com

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_em_field (f1, f2) result (is_eq)

implicit none

type (em_field_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = all(f1%E == f2%E) .and. all(f1%b == f2%b) .and. &
    all(f1%dE == f2%dE) .and. all(f1%dB == f2%dB) 

end function eq_em_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_mode3 (f1, f2) result (is_eq)

implicit none

type (mode3_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = all(f1%v == f2%v) .and. (f1%a == f2%a) .and. (f1%b == f2%b) .and. &
        (f1%c == f2%c) .and. (f1%x == f2%x) .and. (f1%y == f2%y) 

end function eq_mode3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_space_charge (f1, f2) result (is_eq)

implicit none

type (space_charge_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%closed_orb == f2%closed_orb) .and. (f1%kick_const == f2%kick_const) .and. &
        (f1%sig_x == f2%sig_x) .and. (f1%sig_y == f2%sig_y) .and. (f1%phi == f2%phi) .and. &
        (f1%sin_phi == f2%sin_phi) .and. (f1%cos_phi == f2%cos_phi) .and. (f1%sig_z == f2%sig_z)

end function eq_space_charge

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_wall3d (f1, f2) result (is_eq)

implicit none

type (wall3d_struct), intent(in) :: f1, f2
logical is_eq
integer i

!

is_eq = (associated(f1%section) .eqv. associated(f2%section))
if (.not. is_eq) return

is_eq = .false.

if (associated (f1%section)) then
  if (size(f1%section) /= size(f2%section)) return
  do i = 1, size(f1%section)
    if (.not. (f1%section(i) == f2%section(i))) return
  enddo
endif

is_eq = .true.

end function eq_wall3d

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_wall3d_section (f1, f2) result (is_eq)

implicit none

type (wall3d_section_struct), intent(in) :: f1, f2
logical is_eq
integer i

!

is_eq = (allocated(f1%v) .eqv. allocated(f2%v)) .and. (f1%type == f2%type) .and. &
        (f1%s == f2%s) .and. (f1%x0 == f2%x0) .and. (f1%dx0_ds == f2%dx0_ds) .and. &
        (f1%dr_ds == f2%dr_ds) .and. &
        all(f1%p1_coef == f2%p1_coef) .and. all(f1%p2_coef == f2%p2_coef) .and. &
        all(f1%x0_coef == f2%x0_coef) .and. all(f1%y0_coef== f2%y0_coef) 

if (.not. is_eq) return

is_eq = .false.

if (allocated(f1%v)) then
  if (size(f1%v) /= size(f2%v)) return
  do i = 1, size(f1%v)
    if (.not. (f1%v(i) == f2%v(i))) return
  enddo
endif

is_eq = .true.

end function eq_wall3d_section

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_wall3d_vertex (f1, f2) result (is_eq)

implicit none

type (wall3d_vertex_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%x == f2%x) .and. (f1%y == f2%y) .and. (f1%radius_x == f2%radius_x) .and. &
        (f1%radius_y == f2%radius_y) .and. (f1%tilt == f2%tilt) .and. (f1%angle == f2%angle) .and. &
        (f1%x0 == f2%x0) .and. (f1%y0 == f2%y0)

end function eq_wall3d_vertex

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_ele (f1, f2) result (is_eq)

implicit none

type (ele_struct), intent(in) :: f1, f2
logical is_eq
integer i, j

!

is_eq = .true.
is_eq = is_eq .and. (f1%name == f2%name)
is_eq = is_eq .and. (f1%type == f2%type) 
is_eq = is_eq .and. (f1%alias == f2%alias) 
is_eq = is_eq .and. (f1%component_name == f2%component_name)
is_eq = is_eq .and. (associated(f1%descrip) .eqv. associated(f2%descrip));  if (.not. is_eq) return 
if (associated(f1%descrip)) then
  if (f1%descrip /= f2%descrip) then; is_eq = .false.; return; endif
endif
is_eq = is_eq .and. (f1%a == f2%a) 
is_eq = is_eq .and. (f1%b == f2%b)
is_eq = is_eq .and. (f1%z == f2%z)
is_eq = is_eq .and. (f1%x == f2%x)
is_eq = is_eq .and. (f1%y == f2%y) 
is_eq = is_eq .and. (f1%floor == f2%floor) 
is_eq = is_eq .and. (associated(f1%mode3) .eqv. associated(f2%mode3));  if (.not. is_eq) return 
if (associated(f1%mode3)) then
  if (.not. (f1%mode3 == f2%mode3)) then; is_eq = .false.; return; endif
endif
is_eq = is_eq .and. all(f1%map_ref_orb_in == f2%map_ref_orb_in)
is_eq = is_eq .and. all(f1%map_ref_orb_out == f2%map_ref_orb_out) 
is_eq = is_eq .and. (associated(f1%ptc_genfield) .eqv. associated(f2%ptc_genfield));  if (.not. is_eq) return 
is_eq = is_eq .and. eq_bookkeeping_state(f1%bookkeeping_state, f2%bookkeeping_state) 
do i = 1, size(f1%taylor)
  if (.not. (f1%taylor(i) == f2%taylor(i))) then; is_eq = .false.; return; endif
enddo
is_eq = is_eq .and. (associated(f1%em_field) .eqv. associated(f2%em_field));  if (.not. is_eq) return 
!! check rf field here...
!! is_eq = is_eq .and. (associated(f1%rf%wake) .eqv. associated(f2%rf%wake));  if (.not. is_eq) return 
!! if (associated(f1%rf%wake)) then
!!   if (.not. (f1%rf%wake == f2%rf%wake)) then; is_eq = .false.; return; endif
!! endif
!! is_eq = is_eq .and. (associated(f1%wig%term) .eqv. associated(f2%wig%term));  if (.not. is_eq) return 
!! if (associated(f1%wig%term)) then
!!   if (size(f1%wig%term) /= size(f2%wig%term)) then; is_eq = .false.; return; endif
!!   if (.not. all(f1%wig%term == f2%wig%term)) then; is_eq = .false.; return; endif
!! endif
!! is_eq = is_eq .and. (associated(f1%space_charge) .eqv. associated(f2%space_charge));  if (.not. is_eq) return
!! if (associated(f1%space_charge)) then
!!   if (.not. (f1%space_charge == f2%space_charge)) then; is_eq = .false.; return; endif
!! endif
!! is_eq = is_eq .and. (f1%wall3d == f2%wall3d)
is_eq = is_eq .and. all(f1%value == f2%value) 
is_eq = is_eq .and. all(f1%old_value == f2%old_value) 
is_eq = is_eq .and. all(f1%gen0 == f2%gen0)
is_eq = is_eq .and. all(f1%vec0 == f2%vec0) 
is_eq = is_eq .and. all(f1%mat6 == f2%mat6)
is_eq = is_eq .and. all(f1%c_mat == f2%c_mat) 
is_eq = is_eq .and. (f1%gamma_c == f2%gamma_c)
is_eq = is_eq .and. (f1%s == f2%s) 
is_eq = is_eq .and. (f1%ref_time == f2%ref_time)
is_eq = is_eq .and. (associated(f1%r) .eqv. associated(f2%r));  if (.not. is_eq) return 
if (associated(f1%r)) then
  if (any(lbound(f1%r) /= lbound(f2%r)) .or. any(ubound(f1%r) /= ubound(f2%r))) then; is_eq = .false.; return; endif
  if (any(f1%r /= f2%r)) then; is_eq = .false.; return; endif
endif
is_eq = is_eq .and. (associated(f1%a_pole) .eqv. associated(f2%a_pole));  if (.not. is_eq) return 
is_eq = is_eq .and. (associated(f1%b_pole) .eqv. associated(f2%b_pole));  if (.not. is_eq) return 
if (associated(f1%a_pole)) then
  if (any(f1%a_pole /= f2%a_pole) .or. any(f1%b_pole /= f2%b_pole)) then; is_eq = .false.; return; endif
endif
is_eq = is_eq .and. (f1%key == f2%key)
is_eq = is_eq .and. (f1%sub_key == f2%sub_key) 
is_eq = is_eq .and. (f1%ix_ele == f2%ix_ele)
is_eq = is_eq .and. (f1%ix_branch == f2%ix_branch)
is_eq = is_eq .and. (f1%ix_value == f2%ix_value)
is_eq = is_eq .and. (f1%slave_status == f2%slave_status) 
is_eq = is_eq .and. (f1%n_slave == f2%n_slave)
is_eq = is_eq .and. (f1%ix1_slave == f2%ix1_slave) 
is_eq = is_eq .and. (f1%ix2_slave == f2%ix2_slave)
is_eq = is_eq .and. (f1%lord_status == f2%lord_status)
is_eq = is_eq .and. (f1%n_lord == f2%n_lord) 
is_eq = is_eq .and. (f1%ic1_lord == f2%ic1_lord)
is_eq = is_eq .and. (f1%ic2_lord == f2%ic2_lord) 
is_eq = is_eq .and. (f1%ix_pointer == f2%ix_pointer)
is_eq = is_eq .and. (f1%ixx == f2%ixx) 
is_eq = is_eq .and. (f1%mat6_calc_method == f2%mat6_calc_method) 
is_eq = is_eq .and. (f1%tracking_method == f2%tracking_method) 
is_eq = is_eq .and. (f1%spin_tracking_method == f2%spin_tracking_method) 
is_eq = is_eq .and. (f1%field_calc == f2%field_calc)
is_eq = is_eq .and. (f1%ref_orbit == f2%ref_orbit)
is_eq = is_eq .and. (f1%aperture_at == f2%aperture_at)
is_eq = is_eq .and. (f1%aperture_type == f2%aperture_type) 
is_eq = is_eq .and. (f1%symplectify .eqv. f2%symplectify) 
is_eq = is_eq .and. (f1%mode_flip .eqv. f2%mode_flip)
is_eq = is_eq .and. (f1%multipoles_on .eqv. f2%multipoles_on) 
is_eq = is_eq .and. (f1%scale_multipoles .eqv. f2%scale_multipoles) 
is_eq = is_eq .and. (f1%map_with_offsets .eqv. f2%map_with_offsets)
is_eq = is_eq .and. (f1%field_master .eqv. f2%field_master) 
is_eq = is_eq .and. (f1%reversed .eqv. f2%reversed)
is_eq = is_eq .and. (f1%is_on .eqv. f2%is_on)
is_eq = is_eq .and. (f1%old_is_on .eqv. f2%old_is_on) 
is_eq = is_eq .and. (f1%logic .eqv. f2%logic)
is_eq = is_eq .and. (f1%bmad_logic .eqv. f2%bmad_logic)
is_eq = is_eq .and. (f1%on_a_girder .eqv. f2%on_a_girder) 
is_eq = is_eq .and. (f1%csr_calc_on .eqv. f2%csr_calc_on)
is_eq = is_eq .and. (f1%offset_moves_aperture .eqv. f2%offset_moves_aperture) 

end function eq_ele

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_bookkeeping_state (f1, f2) result (is_eq)

implicit none

type (bookkeeping_state_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = .true.
is_eq = is_eq .and. (f1%mat6 == f2%mat6) 
is_eq = is_eq .and. (f1%control == f2%control) 
is_eq = is_eq .and. (f1%floor_position == f2%floor_position) 
is_eq = is_eq .and. (f1%ref_energy == f2%ref_energy) 
is_eq = is_eq .and. (f1%attributes == f2%attributes) 
is_eq = is_eq .and. (f1%rad_int == f2%rad_int) 

end function eq_bookkeeping_state

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_mode_info (f1, f2) result (is_eq)

implicit none

type (mode_info_struct), intent(in) :: f1, f2
logical is_eq

!

is_eq = (f1%tune == f2%tune) .and. (f1%emit == f2%emit) .and. (f1%chrom == f2%chrom)

end function eq_mode_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_branch (f1, f2) result (is_eq)

implicit none

type (branch_struct), intent(in) :: f1, f2
logical is_eq
integer i

!

is_eq = (f1%name == f2%name) .and. (f1%ix_branch == f2%ix_branch) .and. &
        (f1%ix_from_branch == f2%ix_from_branch) .and. (f1%ix_from_ele == f2%ix_from_ele)
if (.not. is_eq) return

is_eq = (associated(f1%n_ele_track) .eqv. associated(f2%n_ele_track)) .and. &
        (associated(f1%n_ele_max) .eqv. associated(f2%n_ele_max)) .and. &
        (associated(f1%ele) .eqv. associated(f2%ele)) .and. &
        (associated(f1%param) .eqv. associated(f2%param)) .and. &
        (associated(f1%wall3d) .eqv. associated(f2%wall3d))
if (.not. is_eq) return

is_eq = .false.

if (associated(f1%n_ele_track)) then
  if (f1%n_ele_track /= f2%n_ele_track) return
endif

if (associated(f1%n_ele_max)) then
  if (f1%n_ele_max /= f2%n_ele_max) return
endif

if (associated(f1%param)) then
  if (.not. (f1%param == f2%param)) return
endif

if (associated(f1%wall3d)) then
  if (.not. (f1%wall3d == f2%wall3d)) return
endif

if (size(f1%ele) /= size(f2%ele)) return
do i = 0, f1%n_ele_max
  if (.not. (f1%ele(i) == f2%ele(i))) return
enddo

is_eq = .true.


end function eq_branch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function eq_lat (f1, f2) result (is_eq)

implicit none

type (lat_struct), intent(in) :: f1, f2
logical is_eq
integer i

!

is_eq = .true.
is_eq = is_eq .and. (f1%lattice == f2%lattice) 
is_eq = is_eq .and. (f1%input_file_name == f2%input_file_name) 
is_eq = is_eq .and. (f1%title == f2%title) 
is_eq = is_eq .and. (f1%a == f2%a) 
is_eq = is_eq .and. (f1%b == f2%b) 
is_eq = is_eq .and. (f1%z == f2%z) 
is_eq = is_eq .and. (f1%param == f2%param) 
is_eq = is_eq .and. (f1%version == f2%version) 
is_eq = is_eq .and. (f1%n_ele_track == f2%n_ele_track) 
is_eq = is_eq .and. (f1%n_ele_max == f2%n_ele_max) 
is_eq = is_eq .and. (f1%n_control_max == f2%n_control_max) 
is_eq = is_eq .and. (f1%n_ic_max == f2%n_ic_max) 
is_eq = is_eq .and. (f1%input_taylor_order == f2%input_taylor_order) 
is_eq = is_eq .and. (f1%ele_init == f2%ele_init) 
is_eq = is_eq .and. (size(f1%ele) == size(f2%ele)) 
is_eq = is_eq .and. (size(f1%control) == size(f2%control)) 
is_eq = is_eq .and. (size(f1%ic) == size(f2%ic)) 
is_eq = is_eq .and. (f1%wall3d == f2%wall3d)
is_eq = is_eq .and. (allocated(f1%branch) .eqv. allocated(f2%branch))

if (.not. is_eq) return

do i = 0, f1%n_ele_max
  is_eq = is_eq .and. (f1%ele(i) == f2%ele(i))
enddo

do i = 1, size(f1%control)
  is_eq = is_eq .and. (f1%control(i) == f2%control(i))
enddo

do i = 1, size(f1%ic)
  is_eq = is_eq .and. (f1%ic(i) == f2%ic(i))
enddo

if (.not. is_eq) return

is_eq = .false.

if (allocated(f1%branch)) then
  if (size(f1%branch) /= size(f2%branch)) return
  do i = 1, size(f1%branch)
    if (.not. (f1%branch(i) == f2%branch(i))) return
  enddo
endif

is_eq = .true.

end function eq_lat

end module

