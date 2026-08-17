!+
! Module tao_attrib_resolve_mod
!
! *** THIS FILE IS GENERATED. DO NOT EDIT. ***
!
! Regenerate with, from the root of the bmad-ecosystem repo:
!     python3 tao/scripts/generate_attrib_tables.py
!
! Each routine here resolves a structure component name such as
!     "floor_plan%ele_shape(3)%ele_id"
! onto a pointer to that component, replacing the use of a Fortran namelist read on a
! scratch file for the purpose of setting a structure component from a string.
!
! IMPORTANT: the structure passed to a resolver must have the TARGET attribute.
! The obj dummy argument below is declared with TARGET. If the actual argument does not
! also have TARGET then the returned pointer becomes undefined as soon as the resolver
! returns, since the compiler is free to pass a copy. That failure mode is nasty: it can
! appear to work for one structure and segfault for another, and the crash surfaces at
! whatever code runs next rather than at the call. So declare the variable as, eg:
!     type (bmad_common_struct), target :: this_bmad_com
! For a module variable that lacks TARGET, copy it into a local TARGET variable, resolve
! and set through that, then copy back. See tao_set_opti_de_param_cmd for an example.
!
! Structures covered: 49.  Components: 619.
!-

module tao_attrib_resolve_mod

use tao_attrib_ptr_mod
use bmad_struct
use tao_struct
use tao_input_struct
use quick_plot_struct
use opti_de_mod
use geodesic_lm
use spline_mod
use csr_and_space_charge_mod

implicit none

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_aperture_param_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a aperture_param_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_aperture_param_struct (obj, name, ptr, err, why)

type (aperture_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('min_angle')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%min_angle

case ('max_angle')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%max_angle

case ('n_angle')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_angle

case ('n_turn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_turn

case ('x_init')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_init

case ('y_init')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y_init

case ('rel_accuracy')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%rel_accuracy

case ('abs_accuracy')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%abs_accuracy

case ('start_ele')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%start_ele

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN APERTURE_PARAM_STRUCT)'
end select

end subroutine tao_res_aperture_param_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_aperture_param_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a aperture_param_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_aperture_param_struct_slot (obj, name, i_slot, ptr, err, why)

type (aperture_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%min_angle
  case (2);  ptr%r => obj%max_angle
  case (3);  ptr%i => obj%n_angle
  case (4);  ptr%i => obj%n_turn
  case (5);  ptr%r => obj%x_init
  case (6);  ptr%r => obj%y_init
  case (7);  ptr%r => obj%rel_accuracy
  case (8);  ptr%r => obj%abs_accuracy
  case (9);  ptr%str => obj%start_ele
  case default
    err = .true.
    why = 'TOO MANY VALUES. APERTURE_PARAM_STRUCT HAS 9 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('min_angle')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%min_angle

case ('max_angle')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%max_angle

case ('n_angle')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_angle

case ('n_turn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_turn

case ('x_init')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_init

case ('y_init')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y_init

case ('rel_accuracy')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%rel_accuracy

case ('abs_accuracy')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%abs_accuracy

case ('start_ele')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%start_ele

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN APERTURE_PARAM_STRUCT)'
end select

end subroutine tao_res_aperture_param_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_beam_init_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a beam_init_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_beam_init_struct (obj, name, ptr, err, why)

type (beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('position_file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%position_file

case ('distribution_type')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%str => obj%distribution_type(isub)
  else
    ptr%str1 => obj%distribution_type
  endif

case ('spin')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%r => obj%spin(isub)
  else
    ptr%r1 => obj%spin
  endif

case ('ellipse')
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  call tao_res_ellipse_beam_init_struct (obj%ellipse(isub), rest, ptr, err, why)

case ('kv')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_kv_beam_init_struct (obj%kv, rest, ptr, err, why)

case ('grid')
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  call tao_res_grid_beam_init_struct (obj%grid(isub), rest, ptr, err, why)

case ('center_jitter')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 6, err, why)) return
  if (has_sub) then
    ptr%r => obj%center_jitter(isub)
  else
    ptr%r1 => obj%center_jitter
  endif

case ('emit_jitter')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 2, err, why)) return
  if (has_sub) then
    ptr%r => obj%emit_jitter(isub)
  else
    ptr%r1 => obj%emit_jitter
  endif

case ('sig_z_jitter')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sig_z_jitter

case ('sig_pz_jitter')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sig_pz_jitter

case ('n_particle')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_particle

case ('renorm_center')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%renorm_center

case ('renorm_sigma')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%renorm_sigma

case ('random_engine')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%random_engine

case ('random_gauss_converter')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%random_gauss_converter

case ('random_sigma_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%random_sigma_cutoff

case ('a_norm_emit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%a_norm_emit

case ('b_norm_emit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%b_norm_emit

case ('a_emit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%a_emit

case ('b_emit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%b_emit

case ('dpz_dz')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%dpz_dz

case ('center')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 6, err, why)) return
  if (has_sub) then
    ptr%r => obj%center(isub)
  else
    ptr%r1 => obj%center
  endif

case ('t_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%t_offset

case ('dt_bunch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%dt_bunch

case ('sig_z')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sig_z

case ('sig_pz')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sig_pz

case ('bunch_charge')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%bunch_charge

case ('n_bunch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_bunch

case ('ix_turn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_turn

case ('species')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%species

case ('full_6d_coupling_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%full_6d_coupling_calc

case ('use_particle_start')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_particle_start

case ('use_t_coords')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_t_coords

case ('use_z_as_t')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_z_as_t

case ('file_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%file_name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_beam_init_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a beam_init_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)

type (beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%position_file
  case (2);  ptr%str => obj%distribution_type(1)
  case (3);  ptr%str => obj%distribution_type(2)
  case (4);  ptr%str => obj%distribution_type(3)
  case (5);  ptr%r => obj%spin(1)
  case (6);  ptr%r => obj%spin(2)
  case (7);  ptr%r => obj%spin(3)
  case (8);  ptr%i => obj%ellipse(1)%part_per_ellipse
  case (9);  ptr%i => obj%ellipse(1)%n_ellipse
  case (10);  ptr%r => obj%ellipse(1)%sigma_cutoff
  case (11);  ptr%i => obj%ellipse(2)%part_per_ellipse
  case (12);  ptr%i => obj%ellipse(2)%n_ellipse
  case (13);  ptr%r => obj%ellipse(2)%sigma_cutoff
  case (14);  ptr%i => obj%ellipse(3)%part_per_ellipse
  case (15);  ptr%i => obj%ellipse(3)%n_ellipse
  case (16);  ptr%r => obj%ellipse(3)%sigma_cutoff
  case (17);  ptr%i => obj%kv%part_per_phi(1)
  case (18);  ptr%i => obj%kv%part_per_phi(2)
  case (19);  ptr%i => obj%kv%n_i2
  case (20);  ptr%r => obj%kv%a
  case (21);  ptr%i => obj%grid(1)%n_x
  case (22);  ptr%i => obj%grid(1)%n_px
  case (23);  ptr%r => obj%grid(1)%x_min
  case (24);  ptr%r => obj%grid(1)%x_max
  case (25);  ptr%r => obj%grid(1)%px_min
  case (26);  ptr%r => obj%grid(1)%px_max
  case (27);  ptr%i => obj%grid(2)%n_x
  case (28);  ptr%i => obj%grid(2)%n_px
  case (29);  ptr%r => obj%grid(2)%x_min
  case (30);  ptr%r => obj%grid(2)%x_max
  case (31);  ptr%r => obj%grid(2)%px_min
  case (32);  ptr%r => obj%grid(2)%px_max
  case (33);  ptr%i => obj%grid(3)%n_x
  case (34);  ptr%i => obj%grid(3)%n_px
  case (35);  ptr%r => obj%grid(3)%x_min
  case (36);  ptr%r => obj%grid(3)%x_max
  case (37);  ptr%r => obj%grid(3)%px_min
  case (38);  ptr%r => obj%grid(3)%px_max
  case (39);  ptr%r => obj%center_jitter(1)
  case (40);  ptr%r => obj%center_jitter(2)
  case (41);  ptr%r => obj%center_jitter(3)
  case (42);  ptr%r => obj%center_jitter(4)
  case (43);  ptr%r => obj%center_jitter(5)
  case (44);  ptr%r => obj%center_jitter(6)
  case (45);  ptr%r => obj%emit_jitter(1)
  case (46);  ptr%r => obj%emit_jitter(2)
  case (47);  ptr%r => obj%sig_z_jitter
  case (48);  ptr%r => obj%sig_pz_jitter
  case (49);  ptr%i => obj%n_particle
  case (50);  ptr%l => obj%renorm_center
  case (51);  ptr%l => obj%renorm_sigma
  case (52);  ptr%str => obj%random_engine
  case (53);  ptr%str => obj%random_gauss_converter
  case (54);  ptr%r => obj%random_sigma_cutoff
  case (55);  ptr%r => obj%a_norm_emit
  case (56);  ptr%r => obj%b_norm_emit
  case (57);  ptr%r => obj%a_emit
  case (58);  ptr%r => obj%b_emit
  case (59);  ptr%r => obj%dpz_dz
  case (60);  ptr%r => obj%center(1)
  case (61);  ptr%r => obj%center(2)
  case (62);  ptr%r => obj%center(3)
  case (63);  ptr%r => obj%center(4)
  case (64);  ptr%r => obj%center(5)
  case (65);  ptr%r => obj%center(6)
  case (66);  ptr%r => obj%t_offset
  case (67);  ptr%r => obj%dt_bunch
  case (68);  ptr%r => obj%sig_z
  case (69);  ptr%r => obj%sig_pz
  case (70);  ptr%r => obj%bunch_charge
  case (71);  ptr%i => obj%n_bunch
  case (72);  ptr%i => obj%ix_turn
  case (73);  ptr%str => obj%species
  case (74);  ptr%l => obj%full_6d_coupling_calc
  case (75);  ptr%l => obj%use_particle_start
  case (76);  ptr%l => obj%use_t_coords
  case (77);  ptr%l => obj%use_z_as_t
  case (78);  ptr%str => obj%file_name
  case default
    err = .true.
    why = 'TOO MANY VALUES. BEAM_INIT_STRUCT HAS 78 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('position_file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%position_file

case ('distribution_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%distribution_type(1 + i_slot - 1)

case ('spin')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%spin(1 + i_slot - 1)

case ('ellipse')
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, 1, 3, err, why)) return
  call tao_res_ellipse_beam_init_struct_slot (obj%ellipse(isub), rest, i_slot, ptr, err, why)

case ('kv')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_kv_beam_init_struct_slot (obj%kv, rest, i_slot, ptr, err, why)

case ('grid')
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, 1, 3, err, why)) return
  call tao_res_grid_beam_init_struct_slot (obj%grid(isub), rest, i_slot, ptr, err, why)

case ('center_jitter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 6 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%center_jitter(1 + i_slot - 1)

case ('emit_jitter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 2 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%emit_jitter(1 + i_slot - 1)

case ('sig_z_jitter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sig_z_jitter

case ('sig_pz_jitter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sig_pz_jitter

case ('n_particle')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_particle

case ('renorm_center')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%renorm_center

case ('renorm_sigma')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%renorm_sigma

case ('random_engine')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%random_engine

case ('random_gauss_converter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%random_gauss_converter

case ('random_sigma_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%random_sigma_cutoff

case ('a_norm_emit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%a_norm_emit

case ('b_norm_emit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%b_norm_emit

case ('a_emit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%a_emit

case ('b_emit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%b_emit

case ('dpz_dz')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%dpz_dz

case ('center')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 6 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%center(1 + i_slot - 1)

case ('t_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%t_offset

case ('dt_bunch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%dt_bunch

case ('sig_z')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sig_z

case ('sig_pz')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sig_pz

case ('bunch_charge')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%bunch_charge

case ('n_bunch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_bunch

case ('ix_turn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_turn

case ('species')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%species

case ('full_6d_coupling_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%full_6d_coupling_calc

case ('use_particle_start')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_particle_start

case ('use_t_coords')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_t_coords

case ('use_z_as_t')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_z_as_t

case ('file_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%file_name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_beam_init_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_bmad_common_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a bmad_common_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_bmad_common_struct (obj, name, ptr, err, why)

type (bmad_common_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('max_aperture_limit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%max_aperture_limit

case ('d_orb')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 6, err, why)) return
  if (has_sub) then
    ptr%r => obj%d_orb(isub)
  else
    ptr%r1 => obj%d_orb
  endif

case ('default_ds_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%default_ds_step

case ('significant_length')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%significant_length

case ('rel_tol_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%rel_tol_tracking

case ('abs_tol_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%abs_tol_tracking

case ('rel_tol_adaptive_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%rel_tol_adaptive_tracking

case ('abs_tol_adaptive_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%abs_tol_adaptive_tracking

case ('init_ds_adaptive_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%init_ds_adaptive_tracking

case ('min_ds_adaptive_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%min_ds_adaptive_tracking

case ('fatal_ds_adaptive_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%fatal_ds_adaptive_tracking

case ('autoscale_amp_abs_tol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%autoscale_amp_abs_tol

case ('autoscale_amp_rel_tol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%autoscale_amp_rel_tol

case ('autoscale_phase_tol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%autoscale_phase_tol

case ('electric_dipole_moment')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%electric_dipole_moment

case ('synch_rad_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%synch_rad_scale

case ('sad_eps_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sad_eps_scale

case ('sad_amp_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sad_amp_max

case ('sad_n_div_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%sad_n_div_max

case ('taylor_order')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%taylor_order

case ('runge_kutta_order')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%runge_kutta_order

case ('default_integ_order')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%default_integ_order

case ('max_num_runge_kutta_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%max_num_runge_kutta_step

case ('rf_phase_below_transition_ref')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%rf_phase_below_transition_ref

case ('sr_wakes_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%sr_wakes_on

case ('lr_wakes_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%lr_wakes_on

case ('auto_bookkeeper')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%auto_bookkeeper

case ('high_energy_space_charge_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%high_energy_space_charge_on

case ('high_energy_space_charge_linear')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%high_energy_space_charge_linear

case ('csr_and_space_charge_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%csr_and_space_charge_on

case ('spin_tracking_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%spin_tracking_on

case ('spin_sokolov_ternov_flipping_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%spin_sokolov_ternov_flipping_on

case ('radiation_damping_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%radiation_damping_on

case ('radiation_zero_average')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%radiation_zero_average

case ('radiation_fluctuations_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%radiation_fluctuations_on

case ('conserve_taylor_maps')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%conserve_taylor_maps

case ('absolute_time_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%absolute_time_tracking

case ('absolute_time_ref_shift')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%absolute_time_ref_shift

case ('convert_to_kinetic_momentum')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%convert_to_kinetic_momentum

case ('normalize_twiss')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%normalize_twiss

case ('aperture_limit_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%aperture_limit_on

case ('spin_n0_direction_user_set')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%spin_n0_direction_user_set

case ('debug')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%debug

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN BMAD_COMMON_STRUCT)'
end select

end subroutine tao_res_bmad_common_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_bmad_common_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a bmad_common_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_bmad_common_struct_slot (obj, name, i_slot, ptr, err, why)

type (bmad_common_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%max_aperture_limit
  case (2);  ptr%r => obj%d_orb(1)
  case (3);  ptr%r => obj%d_orb(2)
  case (4);  ptr%r => obj%d_orb(3)
  case (5);  ptr%r => obj%d_orb(4)
  case (6);  ptr%r => obj%d_orb(5)
  case (7);  ptr%r => obj%d_orb(6)
  case (8);  ptr%r => obj%default_ds_step
  case (9);  ptr%r => obj%significant_length
  case (10);  ptr%r => obj%rel_tol_tracking
  case (11);  ptr%r => obj%abs_tol_tracking
  case (12);  ptr%r => obj%rel_tol_adaptive_tracking
  case (13);  ptr%r => obj%abs_tol_adaptive_tracking
  case (14);  ptr%r => obj%init_ds_adaptive_tracking
  case (15);  ptr%r => obj%min_ds_adaptive_tracking
  case (16);  ptr%r => obj%fatal_ds_adaptive_tracking
  case (17);  ptr%r => obj%autoscale_amp_abs_tol
  case (18);  ptr%r => obj%autoscale_amp_rel_tol
  case (19);  ptr%r => obj%autoscale_phase_tol
  case (20);  ptr%r => obj%electric_dipole_moment
  case (21);  ptr%r => obj%synch_rad_scale
  case (22);  ptr%r => obj%sad_eps_scale
  case (23);  ptr%r => obj%sad_amp_max
  case (24);  ptr%i => obj%sad_n_div_max
  case (25);  ptr%i => obj%taylor_order
  case (26);  ptr%i => obj%runge_kutta_order
  case (27);  ptr%i => obj%default_integ_order
  case (28);  ptr%i => obj%max_num_runge_kutta_step
  case (29);  ptr%l => obj%rf_phase_below_transition_ref
  case (30);  ptr%l => obj%sr_wakes_on
  case (31);  ptr%l => obj%lr_wakes_on
  case (32);  ptr%l => obj%auto_bookkeeper
  case (33);  ptr%l => obj%high_energy_space_charge_on
  case (34);  ptr%l => obj%high_energy_space_charge_linear
  case (35);  ptr%l => obj%csr_and_space_charge_on
  case (36);  ptr%l => obj%spin_tracking_on
  case (37);  ptr%l => obj%spin_sokolov_ternov_flipping_on
  case (38);  ptr%l => obj%radiation_damping_on
  case (39);  ptr%l => obj%radiation_zero_average
  case (40);  ptr%l => obj%radiation_fluctuations_on
  case (41);  ptr%l => obj%conserve_taylor_maps
  case (42);  ptr%l => obj%absolute_time_tracking
  case (43);  ptr%l => obj%absolute_time_ref_shift
  case (44);  ptr%l => obj%convert_to_kinetic_momentum
  case (45);  ptr%l => obj%normalize_twiss
  case (46);  ptr%l => obj%aperture_limit_on
  case (47);  ptr%l => obj%spin_n0_direction_user_set
  case (48);  ptr%l => obj%debug
  case default
    err = .true.
    why = 'TOO MANY VALUES. BMAD_COMMON_STRUCT HAS 48 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('max_aperture_limit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%max_aperture_limit

case ('d_orb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 6 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%d_orb(1 + i_slot - 1)

case ('default_ds_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%default_ds_step

case ('significant_length')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%significant_length

case ('rel_tol_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%rel_tol_tracking

case ('abs_tol_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%abs_tol_tracking

case ('rel_tol_adaptive_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%rel_tol_adaptive_tracking

case ('abs_tol_adaptive_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%abs_tol_adaptive_tracking

case ('init_ds_adaptive_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%init_ds_adaptive_tracking

case ('min_ds_adaptive_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%min_ds_adaptive_tracking

case ('fatal_ds_adaptive_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%fatal_ds_adaptive_tracking

case ('autoscale_amp_abs_tol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%autoscale_amp_abs_tol

case ('autoscale_amp_rel_tol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%autoscale_amp_rel_tol

case ('autoscale_phase_tol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%autoscale_phase_tol

case ('electric_dipole_moment')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%electric_dipole_moment

case ('synch_rad_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%synch_rad_scale

case ('sad_eps_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sad_eps_scale

case ('sad_amp_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sad_amp_max

case ('sad_n_div_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%sad_n_div_max

case ('taylor_order')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%taylor_order

case ('runge_kutta_order')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%runge_kutta_order

case ('default_integ_order')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%default_integ_order

case ('max_num_runge_kutta_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%max_num_runge_kutta_step

case ('rf_phase_below_transition_ref')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%rf_phase_below_transition_ref

case ('sr_wakes_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%sr_wakes_on

case ('lr_wakes_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%lr_wakes_on

case ('auto_bookkeeper')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%auto_bookkeeper

case ('high_energy_space_charge_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%high_energy_space_charge_on

case ('high_energy_space_charge_linear')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%high_energy_space_charge_linear

case ('csr_and_space_charge_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%csr_and_space_charge_on

case ('spin_tracking_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%spin_tracking_on

case ('spin_sokolov_ternov_flipping_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%spin_sokolov_ternov_flipping_on

case ('radiation_damping_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%radiation_damping_on

case ('radiation_zero_average')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%radiation_zero_average

case ('radiation_fluctuations_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%radiation_fluctuations_on

case ('conserve_taylor_maps')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%conserve_taylor_maps

case ('absolute_time_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%absolute_time_tracking

case ('absolute_time_ref_shift')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%absolute_time_ref_shift

case ('convert_to_kinetic_momentum')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%convert_to_kinetic_momentum

case ('normalize_twiss')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%normalize_twiss

case ('aperture_limit_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%aperture_limit_on

case ('spin_n0_direction_user_set')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%spin_n0_direction_user_set

case ('debug')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%debug

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN BMAD_COMMON_STRUCT)'
end select

end subroutine tao_res_bmad_common_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_ele_pointer_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a ele_pointer_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_ele_pointer_struct (obj, name, ptr, err, why)

type (ele_pointer_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele')
  err = .true.; why = 'COMPONENT ELE ' // &
          'IS A POINTER COMPONENT AND CANNOT BE SET'

case ('loc')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_lat_ele_loc_struct (obj%loc, rest, ptr, err, why)

case ('id')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%id

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN ELE_POINTER_STRUCT)'
end select

end subroutine tao_res_ele_pointer_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_ele_pointer_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a ele_pointer_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_ele_pointer_struct_slot (obj, name, i_slot, ptr, err, why)

type (ele_pointer_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%loc%ix_ele
  case (2);  ptr%i => obj%loc%ix_branch
  case (3);  ptr%i => obj%id
  case default
    err = .true.
    why = 'TOO MANY VALUES. ELE_POINTER_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('loc')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_lat_ele_loc_struct_slot (obj%loc, rest, i_slot, ptr, err, why)

case ('id')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%id

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN ELE_POINTER_STRUCT)'
end select

end subroutine tao_res_ele_pointer_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_ellipse_beam_init_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a ellipse_beam_init_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_ellipse_beam_init_struct (obj, name, ptr, err, why)

type (ellipse_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('part_per_ellipse')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%part_per_ellipse

case ('n_ellipse')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_ellipse

case ('sigma_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%sigma_cutoff

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN ELLIPSE_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_ellipse_beam_init_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_ellipse_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a ellipse_beam_init_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_ellipse_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)

type (ellipse_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%part_per_ellipse
  case (2);  ptr%i => obj%n_ellipse
  case (3);  ptr%r => obj%sigma_cutoff
  case default
    err = .true.
    why = 'TOO MANY VALUES. ELLIPSE_BEAM_INIT_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('part_per_ellipse')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%part_per_ellipse

case ('n_ellipse')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_ellipse

case ('sigma_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%sigma_cutoff

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN ELLIPSE_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_ellipse_beam_init_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_geodesic_lm_param_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a geodesic_lm_param_struct instance to a pointer.
! Structure defined in: sim_utils/geodesic_lm/geodesic_lm.f90
!-

recursive subroutine tao_res_geodesic_lm_param_struct (obj, name, ptr, err, why)

type (geodesic_lm_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('mode')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%mode

case ('maxiter')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%maxiter

case ('maxfev')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%maxfev

case ('maxjev')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%maxjev

case ('maxaev')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%maxaev

case ('print_level')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%print_level

case ('print_unit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%print_unit

case ('imethod')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%imethod

case ('iaccel')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%iaccel

case ('ibold')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ibold

case ('ibroyden')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ibroyden

case ('h1')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%h1

case ('h2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%h2

case ('maxlam')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%maxlam

case ('minlam')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%minlam

case ('artol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%artol

case ('cgoal')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%cgoal

case ('gtol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%gtol

case ('xtol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%xtol

case ('xrtol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%xrtol

case ('ftol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%ftol

case ('frtol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%frtol

case ('initialfactor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%initialfactor

case ('factoraccept')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%factoraccept

case ('factorreject')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%factorreject

case ('avmax')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%avmax

case ('analytic_jac')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%analytic_jac

case ('analytic_avv')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%analytic_avv

case ('center_diff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%center_diff

case ('geo_hit_limit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%geo_hit_limit

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN GEODESIC_LM_PARAM_STRUCT)'
end select

end subroutine tao_res_geodesic_lm_param_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_geodesic_lm_param_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a geodesic_lm_param_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_geodesic_lm_param_struct_slot (obj, name, i_slot, ptr, err, why)

type (geodesic_lm_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%mode
  case (2);  ptr%i => obj%maxiter
  case (3);  ptr%i => obj%maxfev
  case (4);  ptr%i => obj%maxjev
  case (5);  ptr%i => obj%maxaev
  case (6);  ptr%i => obj%print_level
  case (7);  ptr%i => obj%print_unit
  case (8);  ptr%i => obj%imethod
  case (9);  ptr%i => obj%iaccel
  case (10);  ptr%i => obj%ibold
  case (11);  ptr%i => obj%ibroyden
  case (12);  ptr%r => obj%h1
  case (13);  ptr%r => obj%h2
  case (14);  ptr%r => obj%maxlam
  case (15);  ptr%r => obj%minlam
  case (16);  ptr%r => obj%artol
  case (17);  ptr%r => obj%cgoal
  case (18);  ptr%r => obj%gtol
  case (19);  ptr%r => obj%xtol
  case (20);  ptr%r => obj%xrtol
  case (21);  ptr%r => obj%ftol
  case (22);  ptr%r => obj%frtol
  case (23);  ptr%r => obj%initialfactor
  case (24);  ptr%r => obj%factoraccept
  case (25);  ptr%r => obj%factorreject
  case (26);  ptr%r => obj%avmax
  case (27);  ptr%l => obj%analytic_jac
  case (28);  ptr%l => obj%analytic_avv
  case (29);  ptr%l => obj%center_diff
  case (30);  ptr%l => obj%geo_hit_limit
  case default
    err = .true.
    why = 'TOO MANY VALUES. GEODESIC_LM_PARAM_STRUCT HAS 30 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('mode')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%mode

case ('maxiter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%maxiter

case ('maxfev')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%maxfev

case ('maxjev')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%maxjev

case ('maxaev')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%maxaev

case ('print_level')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%print_level

case ('print_unit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%print_unit

case ('imethod')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%imethod

case ('iaccel')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%iaccel

case ('ibold')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ibold

case ('ibroyden')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ibroyden

case ('h1')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%h1

case ('h2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%h2

case ('maxlam')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%maxlam

case ('minlam')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%minlam

case ('artol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%artol

case ('cgoal')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%cgoal

case ('gtol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%gtol

case ('xtol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%xtol

case ('xrtol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%xrtol

case ('ftol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%ftol

case ('frtol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%frtol

case ('initialfactor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%initialfactor

case ('factoraccept')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%factoraccept

case ('factorreject')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%factorreject

case ('avmax')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%avmax

case ('analytic_jac')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%analytic_jac

case ('analytic_avv')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%analytic_avv

case ('center_diff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%center_diff

case ('geo_hit_limit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%geo_hit_limit

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN GEODESIC_LM_PARAM_STRUCT)'
end select

end subroutine tao_res_geodesic_lm_param_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_grid_beam_init_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a grid_beam_init_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_grid_beam_init_struct (obj, name, ptr, err, why)

type (grid_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('n_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_x

case ('n_px')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_px

case ('x_min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_min

case ('x_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_max

case ('px_min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%px_min

case ('px_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%px_max

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN GRID_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_grid_beam_init_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_grid_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a grid_beam_init_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_grid_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)

type (grid_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%n_x
  case (2);  ptr%i => obj%n_px
  case (3);  ptr%r => obj%x_min
  case (4);  ptr%r => obj%x_max
  case (5);  ptr%r => obj%px_min
  case (6);  ptr%r => obj%px_max
  case default
    err = .true.
    why = 'TOO MANY VALUES. GRID_BEAM_INIT_STRUCT HAS 6 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('n_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_x

case ('n_px')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_px

case ('x_min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_min

case ('x_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_max

case ('px_min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%px_min

case ('px_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%px_max

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN GRID_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_grid_beam_init_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_kv_beam_init_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a kv_beam_init_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_kv_beam_init_struct (obj, name, ptr, err, why)

type (kv_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('part_per_phi')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 2, err, why)) return
  if (has_sub) then
    ptr%i => obj%part_per_phi(isub)
  else
    ptr%i1 => obj%part_per_phi
  endif

case ('n_i2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_i2

case ('a')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%a

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN KV_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_kv_beam_init_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_kv_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a kv_beam_init_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_kv_beam_init_struct_slot (obj, name, i_slot, ptr, err, why)

type (kv_beam_init_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%part_per_phi(1)
  case (2);  ptr%i => obj%part_per_phi(2)
  case (3);  ptr%i => obj%n_i2
  case (4);  ptr%r => obj%a
  case default
    err = .true.
    why = 'TOO MANY VALUES. KV_BEAM_INIT_STRUCT HAS 4 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('part_per_phi')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 2 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%part_per_phi(1 + i_slot - 1)

case ('n_i2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_i2

case ('a')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%a

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN KV_BEAM_INIT_STRUCT)'
end select

end subroutine tao_res_kv_beam_init_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_lat_ele_loc_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a lat_ele_loc_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_lat_ele_loc_struct (obj, name, ptr, err, why)

type (lat_ele_loc_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ix_ele')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_ele

case ('ix_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_branch

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN LAT_ELE_LOC_STRUCT)'
end select

end subroutine tao_res_lat_ele_loc_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_lat_ele_loc_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a lat_ele_loc_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_lat_ele_loc_struct_slot (obj, name, i_slot, ptr, err, why)

type (lat_ele_loc_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%ix_ele
  case (2);  ptr%i => obj%ix_branch
  case default
    err = .true.
    why = 'TOO MANY VALUES. LAT_ELE_LOC_STRUCT HAS 2 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ix_ele')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_ele

case ('ix_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_branch

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN LAT_ELE_LOC_STRUCT)'
end select

end subroutine tao_res_lat_ele_loc_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_opti_de_param_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a opti_de_param_struct instance to a pointer.
! Structure defined in: sim_utils/optimizers/opti_de_mod.f90
!-

recursive subroutine tao_res_opti_de_param_struct (obj, name, ptr, err, why)

type (opti_de_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('cr')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%cr

case ('f')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%f

case ('l_best')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%l_best

case ('use_2nd_diff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_2nd_diff

case ('binomial_cross')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%binomial_cross

case ('randomize_f')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%randomize_f

case ('minimize_merit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%minimize_merit

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN OPTI_DE_PARAM_STRUCT)'
end select

end subroutine tao_res_opti_de_param_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_opti_de_param_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a opti_de_param_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_opti_de_param_struct_slot (obj, name, i_slot, ptr, err, why)

type (opti_de_param_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%cr
  case (2);  ptr%r => obj%f
  case (3);  ptr%r => obj%l_best
  case (4);  ptr%l => obj%use_2nd_diff
  case (5);  ptr%l => obj%binomial_cross
  case (6);  ptr%l => obj%randomize_f
  case (7);  ptr%l => obj%minimize_merit
  case default
    err = .true.
    why = 'TOO MANY VALUES. OPTI_DE_PARAM_STRUCT HAS 7 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('cr')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%cr

case ('f')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%f

case ('l_best')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%l_best

case ('use_2nd_diff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_2nd_diff

case ('binomial_cross')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%binomial_cross

case ('randomize_f')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%randomize_f

case ('minimize_merit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%minimize_merit

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN OPTI_DE_PARAM_STRUCT)'
end select

end subroutine tao_res_opti_de_param_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_axis_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_axis_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_axis_struct (obj, name, ptr, err, why)

type (qp_axis_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('label')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%label

case ('min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%min

case ('max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%max

case ('tick_min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%tick_min

case ('tick_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%tick_max

case ('eval_min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%eval_min

case ('eval_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%eval_max

case ('dtick')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%dtick

case ('number_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%number_offset

case ('label_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%label_offset

case ('major_tick_len')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%major_tick_len

case ('minor_tick_len')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%minor_tick_len

case ('label_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%label_color

case ('major_div')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%major_div

case ('major_div_nominal')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%major_div_nominal

case ('minor_div')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%minor_div

case ('minor_div_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%minor_div_max

case ('places')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%places

case ('type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%type

case ('bounds')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%bounds

case ('tick_side')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%tick_side

case ('number_side')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%number_side

case ('draw_label')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_label

case ('draw_numbers')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_numbers

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_AXIS_STRUCT)'
end select

end subroutine tao_res_qp_axis_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_axis_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_axis_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_axis_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_axis_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%label
  case (2);  ptr%r => obj%min
  case (3);  ptr%r => obj%max
  case (4);  ptr%r => obj%tick_min
  case (5);  ptr%r => obj%tick_max
  case (6);  ptr%r => obj%eval_min
  case (7);  ptr%r => obj%eval_max
  case (8);  ptr%r => obj%dtick
  case (9);  ptr%r => obj%number_offset
  case (10);  ptr%r => obj%label_offset
  case (11);  ptr%r => obj%major_tick_len
  case (12);  ptr%r => obj%minor_tick_len
  case (13);  ptr%str => obj%label_color
  case (14);  ptr%i => obj%major_div
  case (15);  ptr%i => obj%major_div_nominal
  case (16);  ptr%i => obj%minor_div
  case (17);  ptr%i => obj%minor_div_max
  case (18);  ptr%i => obj%places
  case (19);  ptr%str => obj%type
  case (20);  ptr%str => obj%bounds
  case (21);  ptr%i => obj%tick_side
  case (22);  ptr%i => obj%number_side
  case (23);  ptr%l => obj%draw_label
  case (24);  ptr%l => obj%draw_numbers
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_AXIS_STRUCT HAS 24 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('label')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%label

case ('min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%min

case ('max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%max

case ('tick_min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%tick_min

case ('tick_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%tick_max

case ('eval_min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%eval_min

case ('eval_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%eval_max

case ('dtick')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%dtick

case ('number_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%number_offset

case ('label_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%label_offset

case ('major_tick_len')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%major_tick_len

case ('minor_tick_len')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%minor_tick_len

case ('label_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%label_color

case ('major_div')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%major_div

case ('major_div_nominal')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%major_div_nominal

case ('minor_div')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%minor_div

case ('minor_div_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%minor_div_max

case ('places')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%places

case ('type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%type

case ('bounds')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%bounds

case ('tick_side')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%tick_side

case ('number_side')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%number_side

case ('draw_label')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_label

case ('draw_numbers')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_numbers

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_AXIS_STRUCT)'
end select

end subroutine tao_res_qp_axis_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_legend_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_legend_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_legend_struct (obj, name, ptr, err, why)

type (qp_legend_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('row_spacing')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%row_spacing

case ('line_length')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%line_length

case ('text_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%text_offset

case ('draw_line')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_line

case ('draw_symbol')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_symbol

case ('draw_text')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_text

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_LEGEND_STRUCT)'
end select

end subroutine tao_res_qp_legend_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_legend_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_legend_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_legend_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_legend_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%row_spacing
  case (2);  ptr%r => obj%line_length
  case (3);  ptr%r => obj%text_offset
  case (4);  ptr%l => obj%draw_line
  case (5);  ptr%l => obj%draw_symbol
  case (6);  ptr%l => obj%draw_text
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_LEGEND_STRUCT HAS 6 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('row_spacing')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%row_spacing

case ('line_length')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%line_length

case ('text_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%text_offset

case ('draw_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_line

case ('draw_symbol')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_symbol

case ('draw_text')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_text

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_LEGEND_STRUCT)'
end select

end subroutine tao_res_qp_legend_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_line_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_line_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_line_struct (obj, name, ptr, err, why)

type (qp_line_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%width

case ('color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%color

case ('pattern')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%pattern

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_LINE_STRUCT)'
end select

end subroutine tao_res_qp_line_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_line_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_line_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_line_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_line_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%i => obj%width
  case (2);  ptr%str => obj%color
  case (3);  ptr%str => obj%pattern
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_LINE_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%width

case ('color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%color

case ('pattern')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%pattern

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_LINE_STRUCT)'
end select

end subroutine tao_res_qp_line_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_point_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_point_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_point_struct (obj, name, ptr, err, why)

type (qp_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x

case ('y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y

case ('units')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%units

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_POINT_STRUCT)'
end select

end subroutine tao_res_qp_point_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_point_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_point_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_point_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%x
  case (2);  ptr%r => obj%y
  case (3);  ptr%str => obj%units
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_POINT_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x

case ('y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y

case ('units')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%units

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_POINT_STRUCT)'
end select

end subroutine tao_res_qp_point_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_rect_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_rect_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_rect_struct (obj, name, ptr, err, why)

type (qp_rect_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x1')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x1

case ('x2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x2

case ('y1')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y1

case ('y2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y2

case ('units')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%units

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_RECT_STRUCT)'
end select

end subroutine tao_res_qp_rect_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_rect_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_rect_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_rect_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_rect_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%x1
  case (2);  ptr%r => obj%x2
  case (3);  ptr%r => obj%y1
  case (4);  ptr%r => obj%y2
  case (5);  ptr%str => obj%units
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_RECT_STRUCT HAS 5 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x1')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x1

case ('x2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x2

case ('y1')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y1

case ('y2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y2

case ('units')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%units

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_RECT_STRUCT)'
end select

end subroutine tao_res_qp_rect_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_symbol_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a qp_symbol_struct instance to a pointer.
! Structure defined in: sim_utils/plot/quick_plot_struct.f90
!-

recursive subroutine tao_res_qp_symbol_struct (obj, name, ptr, err, why)

type (qp_symbol_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%type

case ('height')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%height

case ('color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%color

case ('fill_pattern')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%fill_pattern

case ('line_width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%line_width

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_SYMBOL_STRUCT)'
end select

end subroutine tao_res_qp_symbol_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_qp_symbol_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a qp_symbol_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_qp_symbol_struct_slot (obj, name, i_slot, ptr, err, why)

type (qp_symbol_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%type
  case (2);  ptr%r => obj%height
  case (3);  ptr%str => obj%color
  case (4);  ptr%str => obj%fill_pattern
  case (5);  ptr%i => obj%line_width
  case default
    err = .true.
    why = 'TOO MANY VALUES. QP_SYMBOL_STRUCT HAS 5 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%type

case ('height')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%height

case ('color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%color

case ('fill_pattern')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%fill_pattern

case ('line_width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%line_width

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN QP_SYMBOL_STRUCT)'
end select

end subroutine tao_res_qp_symbol_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_space_charge_common_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a space_charge_common_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_space_charge_common_struct (obj, name, ptr, err, why)

type (space_charge_common_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ds_track_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%ds_track_step

case ('dt_track_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%dt_track_step

case ('cathode_strength_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%cathode_strength_cutoff

case ('rel_tol_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%rel_tol_tracking

case ('abs_tol_tracking')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%abs_tol_tracking

case ('beam_chamber_height')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%beam_chamber_height

case ('lsc_sigma_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lsc_sigma_cutoff

case ('particle_sigma_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%particle_sigma_cutoff

case ('mesh_growth_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%mesh_growth_factor

case ('mesh_shrink_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%mesh_shrink_factor

case ('space_charge_mesh_size')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%i => obj%space_charge_mesh_size(isub)
  else
    ptr%i1 => obj%space_charge_mesh_size
  endif

case ('csr3d_mesh_size')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%i => obj%csr3d_mesh_size(isub)
  else
    ptr%i1 => obj%csr3d_mesh_size
  endif

case ('n_bin')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_bin

case ('particle_bin_span')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%particle_bin_span

case ('n_shield_images')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_shield_images

case ('sc_min_in_bin')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%sc_min_in_bin

case ('lsc_kick_transverse_dependence')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%lsc_kick_transverse_dependence

case ('debug')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%debug

case ('diagnostic_output_file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%diagnostic_output_file

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN SPACE_CHARGE_COMMON_STRUCT)'
end select

end subroutine tao_res_space_charge_common_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_space_charge_common_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a space_charge_common_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_space_charge_common_struct_slot (obj, name, i_slot, ptr, err, why)

type (space_charge_common_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%ds_track_step
  case (2);  ptr%r => obj%dt_track_step
  case (3);  ptr%r => obj%cathode_strength_cutoff
  case (4);  ptr%r => obj%rel_tol_tracking
  case (5);  ptr%r => obj%abs_tol_tracking
  case (6);  ptr%r => obj%beam_chamber_height
  case (7);  ptr%r => obj%lsc_sigma_cutoff
  case (8);  ptr%r => obj%particle_sigma_cutoff
  case (9);  ptr%r => obj%mesh_growth_factor
  case (10);  ptr%r => obj%mesh_shrink_factor
  case (11);  ptr%i => obj%space_charge_mesh_size(1)
  case (12);  ptr%i => obj%space_charge_mesh_size(2)
  case (13);  ptr%i => obj%space_charge_mesh_size(3)
  case (14);  ptr%i => obj%csr3d_mesh_size(1)
  case (15);  ptr%i => obj%csr3d_mesh_size(2)
  case (16);  ptr%i => obj%csr3d_mesh_size(3)
  case (17);  ptr%i => obj%n_bin
  case (18);  ptr%i => obj%particle_bin_span
  case (19);  ptr%i => obj%n_shield_images
  case (20);  ptr%i => obj%sc_min_in_bin
  case (21);  ptr%l => obj%lsc_kick_transverse_dependence
  case (22);  ptr%l => obj%debug
  case (23);  ptr%str => obj%diagnostic_output_file
  case default
    err = .true.
    why = 'TOO MANY VALUES. SPACE_CHARGE_COMMON_STRUCT HAS 23 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ds_track_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%ds_track_step

case ('dt_track_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%dt_track_step

case ('cathode_strength_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%cathode_strength_cutoff

case ('rel_tol_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%rel_tol_tracking

case ('abs_tol_tracking')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%abs_tol_tracking

case ('beam_chamber_height')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%beam_chamber_height

case ('lsc_sigma_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lsc_sigma_cutoff

case ('particle_sigma_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%particle_sigma_cutoff

case ('mesh_growth_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%mesh_growth_factor

case ('mesh_shrink_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%mesh_shrink_factor

case ('space_charge_mesh_size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%space_charge_mesh_size(1 + i_slot - 1)

case ('csr3d_mesh_size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%csr3d_mesh_size(1 + i_slot - 1)

case ('n_bin')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_bin

case ('particle_bin_span')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%particle_bin_span

case ('n_shield_images')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_shield_images

case ('sc_min_in_bin')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%sc_min_in_bin

case ('lsc_kick_transverse_dependence')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%lsc_kick_transverse_dependence

case ('debug')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%debug

case ('diagnostic_output_file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%diagnostic_output_file

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN SPACE_CHARGE_COMMON_STRUCT)'
end select

end subroutine tao_res_space_charge_common_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_spin_axis_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a spin_axis_struct instance to a pointer.
! Structure defined in: bmad/modules/bmad_struct.f90
!-

recursive subroutine tao_res_spin_axis_struct (obj, name, ptr, err, why)

type (spin_axis_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('l')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%r => obj%l(isub)
  else
    ptr%r1 => obj%l
  endif

case ('n0')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%r => obj%n0(isub)
  else
    ptr%r1 => obj%n0
  endif

case ('m')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 3, err, why)) return
  if (has_sub) then
    ptr%r => obj%m(isub)
  else
    ptr%r1 => obj%m
  endif

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN SPIN_AXIS_STRUCT)'
end select

end subroutine tao_res_spin_axis_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_spin_axis_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a spin_axis_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_spin_axis_struct_slot (obj, name, i_slot, ptr, err, why)

type (spin_axis_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%l(1)
  case (2);  ptr%r => obj%l(2)
  case (3);  ptr%r => obj%l(3)
  case (4);  ptr%r => obj%n0(1)
  case (5);  ptr%r => obj%n0(2)
  case (6);  ptr%r => obj%n0(3)
  case (7);  ptr%r => obj%m(1)
  case (8);  ptr%r => obj%m(2)
  case (9);  ptr%r => obj%m(3)
  case default
    err = .true.
    why = 'TOO MANY VALUES. SPIN_AXIS_STRUCT HAS 9 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('l')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%l(1 + i_slot - 1)

case ('n0')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%n0(1 + i_slot - 1)

case ('m')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 3 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%m(1 + i_slot - 1)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN SPIN_AXIS_STRUCT)'
end select

end subroutine tao_res_spin_axis_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_building_wall_point_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_building_wall_point_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_building_wall_point_struct (obj, name, ptr, err, why)

type (tao_building_wall_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('z')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%z

case ('x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x

case ('radius')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%radius

case ('z_center')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%z_center

case ('x_center')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_center

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_BUILDING_WALL_POINT_STRUCT)'
end select

end subroutine tao_res_tao_building_wall_point_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_building_wall_point_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_building_wall_point_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_building_wall_point_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_building_wall_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%z
  case (2);  ptr%r => obj%x
  case (3);  ptr%r => obj%radius
  case (4);  ptr%r => obj%z_center
  case (5);  ptr%r => obj%x_center
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_BUILDING_WALL_POINT_STRUCT HAS 5 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('z')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%z

case ('x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x

case ('radius')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%radius

case ('z_center')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%z_center

case ('x_center')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_center

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_BUILDING_WALL_POINT_STRUCT)'
end select

end subroutine tao_res_tao_building_wall_point_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_color_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_curve_color_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_curve_color_struct (obj, name, ptr, err, why)

type (tao_curve_color_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('data_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type

case ('is_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%is_on

case ('min')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%min

case ('max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%max

case ('autoscale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_COLOR_STRUCT)'
end select

end subroutine tao_res_tao_curve_color_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_color_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_curve_color_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_curve_color_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_curve_color_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%data_type
  case (2);  ptr%l => obj%is_on
  case (3);  ptr%r => obj%min
  case (4);  ptr%r => obj%max
  case (5);  ptr%l => obj%autoscale
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_CURVE_COLOR_STRUCT HAS 5 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('data_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type

case ('is_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%is_on

case ('min')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%min

case ('max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%max

case ('autoscale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_COLOR_STRUCT)'
end select

end subroutine tao_res_tao_curve_color_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_curve_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_curve_input (obj, name, ptr, err, why)

type (tao_curve_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('data_source')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_source

case ('data_type_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type_x

case ('data_type_z')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type_z

case ('data_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type

case ('data_index')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_index

case ('legend_text')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%legend_text

case ('units')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%units

case ('component')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%component

case ('y_axis_scale_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y_axis_scale_factor

case ('z_color0')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%z_color0

case ('z_color1')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%z_color1

case ('symbol_every')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%symbol_every

case ('ix_universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_universe

case ('n_turn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_turn

case ('draw_line')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_line

case ('draw_symbols')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_symbols

case ('draw_symbol_index')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_symbol_index

case ('draw_error_bars')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_error_bars

case ('use_y2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_y2

case ('use_z_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_z_color

case ('autoscale_z_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_z_color

case ('smooth_line_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%smooth_line_calc

case ('ele_ref_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_ref_name

case ('ix_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_branch

case ('ix_bunch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_bunch

case ('line')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_line_struct (obj%line, rest, ptr, err, why)

case ('symbol')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_symbol_struct (obj%symbol, rest, ptr, err, why)

case ('hist')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_histogram_struct (obj%hist, rest, ptr, err, why)

case ('orbit')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_curve_orbit_struct (obj%orbit, rest, ptr, err, why)

case ('z_color')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_curve_color_struct (obj%z_color, rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_INPUT)'
end select

end subroutine tao_res_tao_curve_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_curve_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_curve_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_curve_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case (2);  ptr%str => obj%data_source
  case (3);  ptr%str => obj%data_type_x
  case (4);  ptr%str => obj%data_type_z
  case (5);  ptr%str => obj%data_type
  case (6);  ptr%str => obj%data_index
  case (7);  ptr%str => obj%legend_text
  case (8);  ptr%str => obj%units
  case (9);  ptr%str => obj%component
  case (10);  ptr%r => obj%y_axis_scale_factor
  case (11);  ptr%r => obj%z_color0
  case (12);  ptr%r => obj%z_color1
  case (13);  ptr%i => obj%symbol_every
  case (14);  ptr%i => obj%ix_universe
  case (15);  ptr%i => obj%n_turn
  case (16);  ptr%l => obj%draw_line
  case (17);  ptr%l => obj%draw_symbols
  case (18);  ptr%l => obj%draw_symbol_index
  case (19);  ptr%l => obj%draw_error_bars
  case (20);  ptr%l => obj%use_y2
  case (21);  ptr%l => obj%use_z_color
  case (22);  ptr%l => obj%autoscale_z_color
  case (23);  ptr%l => obj%smooth_line_calc
  case (24);  ptr%str => obj%ele_ref_name
  case (25);  ptr%i => obj%ix_branch
  case (26);  ptr%i => obj%ix_bunch
  case (27);  ptr%i => obj%line%width
  case (28);  ptr%str => obj%line%color
  case (29);  ptr%str => obj%line%pattern
  case (30);  ptr%str => obj%symbol%type
  case (31);  ptr%r => obj%symbol%height
  case (32);  ptr%str => obj%symbol%color
  case (33);  ptr%str => obj%symbol%fill_pattern
  case (34);  ptr%i => obj%symbol%line_width
  case (35);  ptr%l => obj%hist%density_normalized
  case (36);  ptr%l => obj%hist%weight_by_charge
  case (37);  ptr%r => obj%hist%minimum
  case (38);  ptr%r => obj%hist%maximum
  case (39);  ptr%r => obj%hist%width
  case (40);  ptr%r => obj%hist%center
  case (41);  ptr%i => obj%hist%number
  case (42);  ptr%r => obj%orbit%x
  case (43);  ptr%r => obj%orbit%y
  case (44);  ptr%r => obj%orbit%t
  case (45);  ptr%str => obj%z_color%data_type
  case (46);  ptr%l => obj%z_color%is_on
  case (47);  ptr%r => obj%z_color%min
  case (48);  ptr%r => obj%z_color%max
  case (49);  ptr%l => obj%z_color%autoscale
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_CURVE_INPUT HAS 49 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('data_source')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_source

case ('data_type_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type_x

case ('data_type_z')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type_z

case ('data_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type

case ('data_index')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_index

case ('legend_text')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%legend_text

case ('units')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%units

case ('component')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%component

case ('y_axis_scale_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y_axis_scale_factor

case ('z_color0')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%z_color0

case ('z_color1')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%z_color1

case ('symbol_every')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%symbol_every

case ('ix_universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_universe

case ('n_turn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_turn

case ('draw_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_line

case ('draw_symbols')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_symbols

case ('draw_symbol_index')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_symbol_index

case ('draw_error_bars')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_error_bars

case ('use_y2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_y2

case ('use_z_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_z_color

case ('autoscale_z_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_z_color

case ('smooth_line_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%smooth_line_calc

case ('ele_ref_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_ref_name

case ('ix_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_branch

case ('ix_bunch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_bunch

case ('line')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_line_struct_slot (obj%line, rest, i_slot, ptr, err, why)

case ('symbol')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_symbol_struct_slot (obj%symbol, rest, i_slot, ptr, err, why)

case ('hist')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_histogram_struct_slot (obj%hist, rest, i_slot, ptr, err, why)

case ('orbit')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_curve_orbit_struct_slot (obj%orbit, rest, i_slot, ptr, err, why)

case ('z_color')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_curve_color_struct_slot (obj%z_color, rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_INPUT)'
end select

end subroutine tao_res_tao_curve_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_orbit_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_curve_orbit_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_curve_orbit_struct (obj, name, ptr, err, why)

type (tao_curve_orbit_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x

case ('y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y

case ('t')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%t

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_ORBIT_STRUCT)'
end select

end subroutine tao_res_tao_curve_orbit_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_orbit_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_curve_orbit_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_curve_orbit_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_curve_orbit_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%x
  case (2);  ptr%r => obj%y
  case (3);  ptr%r => obj%t
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_CURVE_ORBIT_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x

case ('y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y

case ('t')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%t

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_ORBIT_STRUCT)'
end select

end subroutine tao_res_tao_curve_orbit_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_curve_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_curve_struct (obj, name, ptr, err, why)

type (tao_curve_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('data_source')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_source

case ('data_index')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_index

case ('data_type_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type_x

case ('data_type')
  if (tao_res_alloc_bad(head, allocated(obj%data_type), err, why)) return
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type

case ('ele_ref_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_ref_name

case ('legend_text')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%legend_text

case ('message_text')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%message_text

case ('component')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%component

case ('why_invalid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%why_invalid

case ('g')
  err = .true.; why = 'COMPONENT G ' // &
          'IS A POINTER COMPONENT AND CANNOT BE SET'

case ('hist')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_histogram_struct (obj%hist, rest, ptr, err, why)

case ('z_color')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_curve_color_struct (obj%z_color, rest, ptr, err, why)

case ('x_line')
  if (tao_res_alloc_bad(head, allocated(obj%x_line), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%x_line,1), ubound(obj%x_line,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%x_line(isub)
  else
    ptr%r1 => obj%x_line
  endif

case ('y_line')
  if (tao_res_alloc_bad(head, allocated(obj%y_line), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%y_line,1), ubound(obj%y_line,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%y_line(isub)
  else
    ptr%r1 => obj%y_line
  endif

case ('y2_line')
  if (tao_res_alloc_bad(head, allocated(obj%y2_line), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%y2_line,1), ubound(obj%y2_line,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%y2_line(isub)
  else
    ptr%r1 => obj%y2_line
  endif

case ('ix_line')
  if (tao_res_alloc_bad(head, allocated(obj%ix_line), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%ix_line,1), ubound(obj%ix_line,1), err, why)) return
  if (has_sub) then
    ptr%i => obj%ix_line(isub)
  else
    ptr%i1 => obj%ix_line
  endif

case ('x_symb')
  if (tao_res_alloc_bad(head, allocated(obj%x_symb), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%x_symb,1), ubound(obj%x_symb,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%x_symb(isub)
  else
    ptr%r1 => obj%x_symb
  endif

case ('y_symb')
  if (tao_res_alloc_bad(head, allocated(obj%y_symb), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%y_symb,1), ubound(obj%y_symb,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%y_symb(isub)
  else
    ptr%r1 => obj%y_symb
  endif

case ('z_symb')
  if (tao_res_alloc_bad(head, allocated(obj%z_symb), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%z_symb,1), ubound(obj%z_symb,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%z_symb(isub)
  else
    ptr%r1 => obj%z_symb
  endif

case ('err_symb')
  if (tao_res_alloc_bad(head, allocated(obj%err_symb), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%err_symb,1), ubound(obj%err_symb,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%err_symb(isub)
  else
    ptr%r1 => obj%err_symb
  endif

case ('symb_size')
  if (tao_res_alloc_bad(head, allocated(obj%symb_size), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%symb_size,1), ubound(obj%symb_size,1), err, why)) return
  if (has_sub) then
    ptr%r => obj%symb_size(isub)
  else
    ptr%r1 => obj%symb_size
  endif

case ('ix_symb')
  if (tao_res_alloc_bad(head, allocated(obj%ix_symb), err, why)) return
  if (tao_res_array_bad(head, rest, has_sub, isub, lbound(obj%ix_symb,1), ubound(obj%ix_symb,1), err, why)) return
  if (has_sub) then
    ptr%i => obj%ix_symb(isub)
  else
    ptr%i1 => obj%ix_symb
  endif

case ('y_axis_scale_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y_axis_scale_factor

case ('line')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_line_struct (obj%line, rest, ptr, err, why)

case ('symbol')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_symbol_struct (obj%symbol, rest, ptr, err, why)

case ('orbit')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_curve_orbit_struct (obj%orbit, rest, ptr, err, why)

case ('ix_universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_universe

case ('symbol_every')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%symbol_every

case ('ix_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_branch

case ('ix_bunch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_bunch

case ('n_turn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_turn

case ('use_y2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%use_y2

case ('draw_line')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_line

case ('draw_symbols')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_symbols

case ('draw_symbol_index')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_symbol_index

case ('draw_error_bars')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_error_bars

case ('smooth_line_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%smooth_line_calc

case ('valid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%valid

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_STRUCT)'
end select

end subroutine tao_res_tao_curve_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_curve_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_curve_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_curve_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_curve_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_CURVE_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('data_source')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_source

case ('data_index')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_index

case ('data_type_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type_x

case ('data_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type

case ('ele_ref_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_ref_name

case ('legend_text')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%legend_text

case ('message_text')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%message_text

case ('component')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%component

case ('why_invalid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%why_invalid

case ('hist')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_histogram_struct_slot (obj%hist, rest, i_slot, ptr, err, why)

case ('z_color')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_curve_color_struct_slot (obj%z_color, rest, i_slot, ptr, err, why)

case ('x_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%x_line), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%x_line,1) - lbound(obj%x_line,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_line(lbound(obj%x_line,1) + i_slot - 1)

case ('y_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%y_line), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%y_line,1) - lbound(obj%y_line,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y_line(lbound(obj%y_line,1) + i_slot - 1)

case ('y2_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%y2_line), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%y2_line,1) - lbound(obj%y2_line,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y2_line(lbound(obj%y2_line,1) + i_slot - 1)

case ('ix_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%ix_line), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%ix_line,1) - lbound(obj%ix_line,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_line(lbound(obj%ix_line,1) + i_slot - 1)

case ('x_symb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%x_symb), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%x_symb,1) - lbound(obj%x_symb,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_symb(lbound(obj%x_symb,1) + i_slot - 1)

case ('y_symb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%y_symb), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%y_symb,1) - lbound(obj%y_symb,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y_symb(lbound(obj%y_symb,1) + i_slot - 1)

case ('z_symb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%z_symb), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%z_symb,1) - lbound(obj%z_symb,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%z_symb(lbound(obj%z_symb,1) + i_slot - 1)

case ('err_symb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%err_symb), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%err_symb,1) - lbound(obj%err_symb,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%err_symb(lbound(obj%err_symb,1) + i_slot - 1)

case ('symb_size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%symb_size), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%symb_size,1) - lbound(obj%symb_size,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%symb_size(lbound(obj%symb_size,1) + i_slot - 1)

case ('ix_symb')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (tao_res_alloc_bad(head, allocated(obj%ix_symb), err, why)) return
  if (i_slot < 1 .or. i_slot > ubound(obj%ix_symb,1) - lbound(obj%ix_symb,1) + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_symb(lbound(obj%ix_symb,1) + i_slot - 1)

case ('y_axis_scale_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y_axis_scale_factor

case ('line')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_line_struct_slot (obj%line, rest, i_slot, ptr, err, why)

case ('symbol')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_symbol_struct_slot (obj%symbol, rest, i_slot, ptr, err, why)

case ('orbit')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_curve_orbit_struct_slot (obj%orbit, rest, i_slot, ptr, err, why)

case ('ix_universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_universe

case ('symbol_every')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%symbol_every

case ('ix_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_branch

case ('ix_bunch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_bunch

case ('n_turn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_turn

case ('use_y2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%use_y2

case ('draw_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_line

case ('draw_symbols')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_symbols

case ('draw_symbol_index')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_symbol_index

case ('draw_error_bars')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_error_bars

case ('smooth_line_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%smooth_line_calc

case ('valid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%valid

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_CURVE_STRUCT)'
end select

end subroutine tao_res_tao_curve_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_d1_data_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_d1_data_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_d1_data_input (obj, name, ptr, err, why)

type (tao_d1_data_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_D1_DATA_INPUT)'
end select

end subroutine tao_res_tao_d1_data_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_d1_data_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_d1_data_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_d1_data_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_d1_data_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_D1_DATA_INPUT HAS 1 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_D1_DATA_INPUT)'
end select

end subroutine tao_res_tao_d1_data_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_d2_data_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_d2_data_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_d2_data_input (obj, name, ptr, err, why)

type (tao_d2_data_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_D2_DATA_INPUT)'
end select

end subroutine tao_res_tao_d2_data_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_d2_data_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_d2_data_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_d2_data_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_d2_data_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_D2_DATA_INPUT HAS 1 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_D2_DATA_INPUT)'
end select

end subroutine tao_res_tao_d2_data_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_datum_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_datum_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_datum_input (obj, name, ptr, err, why)

type (tao_datum_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('data_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_type

case ('ele_ref_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_ref_name

case ('ele_start_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_start_name

case ('ele_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_name

case ('merit_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%merit_type

case ('meas')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%meas

case ('weight')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%weight

case ('good_user')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%good_user

case ('good_opt')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%good_opt

case ('data_source')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%data_source

case ('eval_point')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%eval_point

case ('s_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%s_offset

case ('ref_s_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%ref_s_offset

case ('ix_bunch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_bunch

case ('spin_axis')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_spin_axis_struct (obj%spin_axis, rest, ptr, err, why)

case ('invalid_value')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%invalid_value

case ('error_rms')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%error_rms

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DATUM_INPUT)'
end select

end subroutine tao_res_tao_datum_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_datum_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_datum_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_datum_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_datum_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%data_type
  case (2);  ptr%str => obj%ele_ref_name
  case (3);  ptr%str => obj%ele_start_name
  case (4);  ptr%str => obj%ele_name
  case (5);  ptr%str => obj%merit_type
  case (6);  ptr%r => obj%meas
  case (7);  ptr%r => obj%weight
  case (8);  ptr%l => obj%good_user
  case (9);  ptr%l => obj%good_opt
  case (10);  ptr%str => obj%data_source
  case (11);  ptr%str => obj%eval_point
  case (12);  ptr%r => obj%s_offset
  case (13);  ptr%r => obj%ref_s_offset
  case (14);  ptr%i => obj%ix_bunch
  case (15);  ptr%r => obj%spin_axis%l(1)
  case (16);  ptr%r => obj%spin_axis%l(2)
  case (17);  ptr%r => obj%spin_axis%l(3)
  case (18);  ptr%r => obj%spin_axis%n0(1)
  case (19);  ptr%r => obj%spin_axis%n0(2)
  case (20);  ptr%r => obj%spin_axis%n0(3)
  case (21);  ptr%r => obj%spin_axis%m(1)
  case (22);  ptr%r => obj%spin_axis%m(2)
  case (23);  ptr%r => obj%spin_axis%m(3)
  case (24);  ptr%r => obj%invalid_value
  case (25);  ptr%r => obj%error_rms
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_DATUM_INPUT HAS 25 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('data_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_type

case ('ele_ref_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_ref_name

case ('ele_start_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_start_name

case ('ele_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_name

case ('merit_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%merit_type

case ('meas')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%meas

case ('weight')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%weight

case ('good_user')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%good_user

case ('good_opt')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%good_opt

case ('data_source')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%data_source

case ('eval_point')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%eval_point

case ('s_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%s_offset

case ('ref_s_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%ref_s_offset

case ('ix_bunch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_bunch

case ('spin_axis')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_spin_axis_struct_slot (obj%spin_axis, rest, i_slot, ptr, err, why)

case ('invalid_value')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%invalid_value

case ('error_rms')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%error_rms

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DATUM_INPUT)'
end select

end subroutine tao_res_tao_datum_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_design_lat_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_design_lat_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_design_lat_input (obj, name, ptr, err, why)

type (tao_design_lat_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%file

case ('file2')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%file2

case ('language')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%language

case ('use_line')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%use_line

case ('one_turn_map_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%one_turn_map_calc

case ('dynamic_aperture_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%dynamic_aperture_calc

case ('reverse_lattice')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%reverse_lattice

case ('start_branch_at')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%start_branch_at

case ('slice_lattice')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%slice_lattice

case ('use_element_range')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 2, err, why)) return
  if (has_sub) then
    ptr%str => obj%use_element_range(isub)
  else
    ptr%str1 => obj%use_element_range
  endif

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DESIGN_LAT_INPUT)'
end select

end subroutine tao_res_tao_design_lat_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_design_lat_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_design_lat_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_design_lat_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_design_lat_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%file
  case (2);  ptr%str => obj%file2
  case (3);  ptr%str => obj%language
  case (4);  ptr%str => obj%use_line
  case (5);  ptr%l => obj%one_turn_map_calc
  case (6);  ptr%l => obj%dynamic_aperture_calc
  case (7);  ptr%l => obj%reverse_lattice
  case (8);  ptr%str => obj%start_branch_at
  case (9);  ptr%str => obj%slice_lattice
  case (10);  ptr%str => obj%use_element_range(1)
  case (11);  ptr%str => obj%use_element_range(2)
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_DESIGN_LAT_INPUT HAS 11 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%file

case ('file2')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%file2

case ('language')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%language

case ('use_line')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%use_line

case ('one_turn_map_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%one_turn_map_calc

case ('dynamic_aperture_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%dynamic_aperture_calc

case ('reverse_lattice')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%reverse_lattice

case ('start_branch_at')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%start_branch_at

case ('slice_lattice')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%slice_lattice

case ('use_element_range')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 2 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%use_element_range(1 + i_slot - 1)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DESIGN_LAT_INPUT)'
end select

end subroutine tao_res_tao_design_lat_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_drawing_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_drawing_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_drawing_struct (obj, name, ptr, err, why)

type (tao_drawing_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_shape')
  if (tao_res_alloc_bad(head, allocated(obj%ele_shape), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%ele_shape,1), ubound(obj%ele_shape,1), err, why)) return
  call tao_res_tao_ele_shape_struct (obj%ele_shape(isub), rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DRAWING_STRUCT)'
end select

end subroutine tao_res_tao_drawing_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_drawing_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_drawing_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_drawing_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_drawing_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_DRAWING_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_shape')
  if (tao_res_alloc_bad(head, allocated(obj%ele_shape), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%ele_shape,1), ubound(obj%ele_shape,1), err, why)) return
  call tao_res_tao_ele_shape_struct_slot (obj%ele_shape(isub), rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_DRAWING_STRUCT)'
end select

end subroutine tao_res_tao_drawing_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_pointer_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_ele_pointer_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_ele_pointer_struct (obj, name, ptr, err, why)

type (tao_ele_pointer_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('eles')
  if (tao_res_alloc_bad(head, allocated(obj%eles), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%eles,1), ubound(obj%eles,1), err, why)) return
  call tao_res_ele_pointer_struct (obj%eles(isub), rest, ptr, err, why)

case ('n_loc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_loc

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_POINTER_STRUCT)'
end select

end subroutine tao_res_tao_ele_pointer_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_pointer_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_ele_pointer_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_ele_pointer_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_ele_pointer_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_ELE_POINTER_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('eles')
  if (tao_res_alloc_bad(head, allocated(obj%eles), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%eles,1), ubound(obj%eles,1), err, why)) return
  call tao_res_ele_pointer_struct_slot (obj%eles(isub), rest, i_slot, ptr, err, why)

case ('n_loc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_loc

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_POINTER_STRUCT)'
end select

end subroutine tao_res_tao_ele_pointer_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_shape_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_ele_shape_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_ele_shape_input (obj, name, ptr, err, why)

type (tao_ele_shape_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_id')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_id

case ('shape')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%shape

case ('color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%color

case ('size')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%size

case ('label')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%label

case ('draw')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw

case ('multi')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%multi

case ('line_width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%line_width

case ('offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%offset

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_SHAPE_INPUT)'
end select

end subroutine tao_res_tao_ele_shape_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_shape_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_ele_shape_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_ele_shape_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_ele_shape_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%ele_id
  case (2);  ptr%str => obj%shape
  case (3);  ptr%str => obj%color
  case (4);  ptr%r => obj%size
  case (5);  ptr%str => obj%label
  case (6);  ptr%l => obj%draw
  case (7);  ptr%l => obj%multi
  case (8);  ptr%i => obj%line_width
  case (9);  ptr%r => obj%offset
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_ELE_SHAPE_INPUT HAS 9 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_id')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_id

case ('shape')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%shape

case ('color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%color

case ('size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%size

case ('label')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%label

case ('draw')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw

case ('multi')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%multi

case ('line_width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%line_width

case ('offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%offset

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_SHAPE_INPUT)'
end select

end subroutine tao_res_tao_ele_shape_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_shape_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_ele_shape_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_ele_shape_struct (obj, name, ptr, err, why)

type (tao_ele_shape_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_id')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_id

case ('shape')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%shape

case ('color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%color

case ('size')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%size

case ('label')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%label

case ('draw')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw

case ('multi')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%multi

case ('line_width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%line_width

case ('offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%offset

case ('ix_key')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_key

case ('name_ele')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name_ele

case ('uni')
  if (tao_res_alloc_bad(head, allocated(obj%uni), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%uni,1), ubound(obj%uni,1), err, why)) return
  call tao_res_tao_ele_pointer_struct (obj%uni(isub), rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_SHAPE_STRUCT)'
end select

end subroutine tao_res_tao_ele_shape_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_ele_shape_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_ele_shape_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_ele_shape_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_ele_shape_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_ELE_SHAPE_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_id')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_id

case ('shape')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%shape

case ('color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%color

case ('size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%size

case ('label')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%label

case ('draw')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw

case ('multi')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%multi

case ('line_width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%line_width

case ('offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%offset

case ('ix_key')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_key

case ('name_ele')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name_ele

case ('uni')
  if (tao_res_alloc_bad(head, allocated(obj%uni), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%uni,1), ubound(obj%uni,1), err, why)) return
  call tao_res_tao_ele_pointer_struct_slot (obj%uni(isub), rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_ELE_SHAPE_STRUCT)'
end select

end subroutine tao_res_tao_ele_shape_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_floor_plan_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_floor_plan_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_floor_plan_struct (obj, name, ptr, err, why)

type (tao_floor_plan_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('view')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%view

case ('rotation')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%rotation

case ('correct_distortion')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%correct_distortion

case ('flip_label_side')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%flip_label_side

case ('size_is_absolute')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%size_is_absolute

case ('draw_only_first_pass')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_only_first_pass

case ('draw_building_wall')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_building_wall

case ('orbit_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%orbit_scale

case ('orbit_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%orbit_color

case ('orbit_pattern')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%orbit_pattern

case ('orbit_lattice')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%orbit_lattice

case ('orbit_width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%orbit_width

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_FLOOR_PLAN_STRUCT)'
end select

end subroutine tao_res_tao_floor_plan_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_floor_plan_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_floor_plan_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_floor_plan_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_floor_plan_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%view
  case (2);  ptr%r => obj%rotation
  case (3);  ptr%l => obj%correct_distortion
  case (4);  ptr%l => obj%flip_label_side
  case (5);  ptr%l => obj%size_is_absolute
  case (6);  ptr%l => obj%draw_only_first_pass
  case (7);  ptr%l => obj%draw_building_wall
  case (8);  ptr%r => obj%orbit_scale
  case (9);  ptr%str => obj%orbit_color
  case (10);  ptr%str => obj%orbit_pattern
  case (11);  ptr%str => obj%orbit_lattice
  case (12);  ptr%i => obj%orbit_width
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_FLOOR_PLAN_STRUCT HAS 12 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('view')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%view

case ('rotation')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%rotation

case ('correct_distortion')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%correct_distortion

case ('flip_label_side')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%flip_label_side

case ('size_is_absolute')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%size_is_absolute

case ('draw_only_first_pass')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_only_first_pass

case ('draw_building_wall')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_building_wall

case ('orbit_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%orbit_scale

case ('orbit_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%orbit_color

case ('orbit_pattern')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%orbit_pattern

case ('orbit_lattice')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%orbit_lattice

case ('orbit_width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%orbit_width

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_FLOOR_PLAN_STRUCT)'
end select

end subroutine tao_res_tao_floor_plan_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_global_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_global_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_global_struct (obj, name, ptr, err, why)

type (tao_global_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('beam_dead_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%beam_dead_cutoff

case ('lm_opt_deriv_reinit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lm_opt_deriv_reinit

case ('de_lm_step_ratio')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%de_lm_step_ratio

case ('de_var_to_population_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%de_var_to_population_factor

case ('lmdif_eps')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lmdif_eps

case ('lmdif_negligible_merit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lmdif_negligible_merit

case ('svd_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%svd_cutoff

case ('unstable_penalty')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%unstable_penalty

case ('merit_stop_value')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%merit_stop_value

case ('dmerit_stop_value')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%dmerit_stop_value

case ('random_sigma_cutoff')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%random_sigma_cutoff

case ('delta_e_chrom')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%delta_e_chrom

case ('max_plot_time')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%max_plot_time

case ('default_universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%default_universe

case ('default_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%default_branch

case ('n_opti_cycles')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_opti_cycles

case ('n_opti_loops')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_opti_loops

case ('n_threads')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_threads

case ('phase_units')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%phase_units

case ('bunch_to_plot')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%bunch_to_plot

case ('random_seed')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%random_seed

case ('n_top10_merit')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_top10_merit

case ('srdt_gen_n_slices')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%srdt_gen_n_slices

case ('datum_err_messages_max')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%datum_err_messages_max

case ('srdt_sxt_n_slices')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%srdt_sxt_n_slices

case ('srdt_use_cache')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%srdt_use_cache

case ('quiet')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%quiet

case ('random_engine')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%random_engine

case ('random_gauss_converter')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%random_gauss_converter

case ('track_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%track_type

case ('lat_sigma_calc_uses_emit_from')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%lat_sigma_calc_uses_emit_from

case ('prompt_string')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%prompt_string

case ('prompt_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%prompt_color

case ('optimizer')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%optimizer

case ('print_command')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%print_command

case ('var_out_file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%var_out_file

case ('history_file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%history_file

case ('beam_timer_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%beam_timer_on

case ('box_plots')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%box_plots

case ('blank_line_between_commands')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%blank_line_between_commands

case ('cmd_file_abort_on_error')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%cmd_file_abort_on_error

case ('concatenate_maps')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%concatenate_maps

case ('derivative_recalc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%derivative_recalc

case ('derivative_uses_design')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%derivative_uses_design

case ('disable_smooth_line_calc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%disable_smooth_line_calc

case ('draw_curve_off_scale_warn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_curve_off_scale_warn

case ('external_plotting')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%external_plotting

case ('label_lattice_elements')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%label_lattice_elements

case ('label_keys')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%label_keys

case ('lattice_calc_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%lattice_calc_on

case ('only_limit_opt_vars')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%only_limit_opt_vars

case ('opt_with_ref')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%opt_with_ref

case ('opt_with_base')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%opt_with_base

case ('opt_match_auto_recalc')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%opt_match_auto_recalc

case ('opti_write_var_file')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%opti_write_var_file

case ('optimizer_allow_user_abort')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%optimizer_allow_user_abort

case ('optimizer_var_limit_warn')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%optimizer_var_limit_warn

case ('plot_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%plot_on

case ('rad_int_user_calc_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%rad_int_user_calc_on

case ('rf_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%rf_on

case ('single_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%single_step

case ('stop_on_error')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%stop_on_error

case ('svd_retreat_on_merit_increase')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%svd_retreat_on_merit_increase

case ('var_limits_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%var_limits_on

case ('wait_for_cr_in_single_mode')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%wait_for_cr_in_single_mode

case ('symbol_import')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%symbol_import

case ('debug_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%debug_on

case ('expression_tree_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%expression_tree_on

case ('verbose_on')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%verbose_on

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GLOBAL_STRUCT)'
end select

end subroutine tao_res_tao_global_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_global_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_global_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_global_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_global_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%beam_dead_cutoff
  case (2);  ptr%r => obj%lm_opt_deriv_reinit
  case (3);  ptr%r => obj%de_lm_step_ratio
  case (4);  ptr%r => obj%de_var_to_population_factor
  case (5);  ptr%r => obj%lmdif_eps
  case (6);  ptr%r => obj%lmdif_negligible_merit
  case (7);  ptr%r => obj%svd_cutoff
  case (8);  ptr%r => obj%unstable_penalty
  case (9);  ptr%r => obj%merit_stop_value
  case (10);  ptr%r => obj%dmerit_stop_value
  case (11);  ptr%r => obj%random_sigma_cutoff
  case (12);  ptr%r => obj%delta_e_chrom
  case (13);  ptr%r => obj%max_plot_time
  case (14);  ptr%i => obj%default_universe
  case (15);  ptr%i => obj%default_branch
  case (16);  ptr%i => obj%n_opti_cycles
  case (17);  ptr%i => obj%n_opti_loops
  case (18);  ptr%i => obj%n_threads
  case (19);  ptr%i => obj%phase_units
  case (20);  ptr%i => obj%bunch_to_plot
  case (21);  ptr%i => obj%random_seed
  case (22);  ptr%i => obj%n_top10_merit
  case (23);  ptr%i => obj%srdt_gen_n_slices
  case (24);  ptr%i => obj%datum_err_messages_max
  case (25);  ptr%i => obj%srdt_sxt_n_slices
  case (26);  ptr%l => obj%srdt_use_cache
  case (27);  ptr%str => obj%quiet
  case (28);  ptr%str => obj%random_engine
  case (29);  ptr%str => obj%random_gauss_converter
  case (30);  ptr%str => obj%track_type
  case (31);  ptr%str => obj%lat_sigma_calc_uses_emit_from
  case (32);  ptr%str => obj%prompt_string
  case (33);  ptr%str => obj%prompt_color
  case (34);  ptr%str => obj%optimizer
  case (35);  ptr%str => obj%print_command
  case (36);  ptr%str => obj%var_out_file
  case (37);  ptr%str => obj%history_file
  case (38);  ptr%l => obj%beam_timer_on
  case (39);  ptr%l => obj%box_plots
  case (40);  ptr%l => obj%blank_line_between_commands
  case (41);  ptr%l => obj%cmd_file_abort_on_error
  case (42);  ptr%l => obj%concatenate_maps
  case (43);  ptr%l => obj%derivative_recalc
  case (44);  ptr%l => obj%derivative_uses_design
  case (45);  ptr%l => obj%disable_smooth_line_calc
  case (46);  ptr%l => obj%draw_curve_off_scale_warn
  case (47);  ptr%l => obj%external_plotting
  case (48);  ptr%l => obj%label_lattice_elements
  case (49);  ptr%l => obj%label_keys
  case (50);  ptr%l => obj%lattice_calc_on
  case (51);  ptr%l => obj%only_limit_opt_vars
  case (52);  ptr%l => obj%opt_with_ref
  case (53);  ptr%l => obj%opt_with_base
  case (54);  ptr%l => obj%opt_match_auto_recalc
  case (55);  ptr%l => obj%opti_write_var_file
  case (56);  ptr%l => obj%optimizer_allow_user_abort
  case (57);  ptr%l => obj%optimizer_var_limit_warn
  case (58);  ptr%l => obj%plot_on
  case (59);  ptr%l => obj%rad_int_user_calc_on
  case (60);  ptr%l => obj%rf_on
  case (61);  ptr%l => obj%single_step
  case (62);  ptr%l => obj%stop_on_error
  case (63);  ptr%l => obj%svd_retreat_on_merit_increase
  case (64);  ptr%l => obj%var_limits_on
  case (65);  ptr%l => obj%wait_for_cr_in_single_mode
  case (66);  ptr%l => obj%symbol_import
  case (67);  ptr%l => obj%debug_on
  case (68);  ptr%l => obj%expression_tree_on
  case (69);  ptr%l => obj%verbose_on
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_GLOBAL_STRUCT HAS 69 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('beam_dead_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%beam_dead_cutoff

case ('lm_opt_deriv_reinit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lm_opt_deriv_reinit

case ('de_lm_step_ratio')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%de_lm_step_ratio

case ('de_var_to_population_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%de_var_to_population_factor

case ('lmdif_eps')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lmdif_eps

case ('lmdif_negligible_merit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lmdif_negligible_merit

case ('svd_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%svd_cutoff

case ('unstable_penalty')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%unstable_penalty

case ('merit_stop_value')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%merit_stop_value

case ('dmerit_stop_value')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%dmerit_stop_value

case ('random_sigma_cutoff')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%random_sigma_cutoff

case ('delta_e_chrom')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%delta_e_chrom

case ('max_plot_time')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%max_plot_time

case ('default_universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%default_universe

case ('default_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%default_branch

case ('n_opti_cycles')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_opti_cycles

case ('n_opti_loops')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_opti_loops

case ('n_threads')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_threads

case ('phase_units')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%phase_units

case ('bunch_to_plot')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%bunch_to_plot

case ('random_seed')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%random_seed

case ('n_top10_merit')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_top10_merit

case ('srdt_gen_n_slices')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%srdt_gen_n_slices

case ('datum_err_messages_max')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%datum_err_messages_max

case ('srdt_sxt_n_slices')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%srdt_sxt_n_slices

case ('srdt_use_cache')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%srdt_use_cache

case ('quiet')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%quiet

case ('random_engine')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%random_engine

case ('random_gauss_converter')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%random_gauss_converter

case ('track_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%track_type

case ('lat_sigma_calc_uses_emit_from')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%lat_sigma_calc_uses_emit_from

case ('prompt_string')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%prompt_string

case ('prompt_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%prompt_color

case ('optimizer')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%optimizer

case ('print_command')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%print_command

case ('var_out_file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%var_out_file

case ('history_file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%history_file

case ('beam_timer_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%beam_timer_on

case ('box_plots')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%box_plots

case ('blank_line_between_commands')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%blank_line_between_commands

case ('cmd_file_abort_on_error')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%cmd_file_abort_on_error

case ('concatenate_maps')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%concatenate_maps

case ('derivative_recalc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%derivative_recalc

case ('derivative_uses_design')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%derivative_uses_design

case ('disable_smooth_line_calc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%disable_smooth_line_calc

case ('draw_curve_off_scale_warn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_curve_off_scale_warn

case ('external_plotting')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%external_plotting

case ('label_lattice_elements')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%label_lattice_elements

case ('label_keys')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%label_keys

case ('lattice_calc_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%lattice_calc_on

case ('only_limit_opt_vars')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%only_limit_opt_vars

case ('opt_with_ref')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%opt_with_ref

case ('opt_with_base')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%opt_with_base

case ('opt_match_auto_recalc')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%opt_match_auto_recalc

case ('opti_write_var_file')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%opti_write_var_file

case ('optimizer_allow_user_abort')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%optimizer_allow_user_abort

case ('optimizer_var_limit_warn')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%optimizer_var_limit_warn

case ('plot_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%plot_on

case ('rad_int_user_calc_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%rad_int_user_calc_on

case ('rf_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%rf_on

case ('single_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%single_step

case ('stop_on_error')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%stop_on_error

case ('svd_retreat_on_merit_increase')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%svd_retreat_on_merit_increase

case ('var_limits_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%var_limits_on

case ('wait_for_cr_in_single_mode')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%wait_for_cr_in_single_mode

case ('symbol_import')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%symbol_import

case ('debug_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%debug_on

case ('expression_tree_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%expression_tree_on

case ('verbose_on')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%verbose_on

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GLOBAL_STRUCT)'
end select

end subroutine tao_res_tao_global_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_graph_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_graph_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_graph_input (obj, name, ptr, err, why)

type (tao_graph_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%type

case ('title')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%title

case ('component')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%component

case ('text_legend')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 10, err, why)) return
  if (has_sub) then
    ptr%str => obj%text_legend(isub)
  else
    ptr%str1 => obj%text_legend
  endif

case ('floor_plan_view')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%floor_plan_view

case ('floor_plan_orbit_color')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%floor_plan_orbit_color

case ('box')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 4, err, why)) return
  if (has_sub) then
    ptr%i => obj%box(isub)
  else
    ptr%i1 => obj%box
  endif

case ('ix_universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_universe

case ('ix_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_branch

case ('n_curve')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_curve

case ('x_axis_scale_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_axis_scale_factor

case ('symbol_size_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%symbol_size_scale

case ('floor_plan_rotation')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_rotation

case ('floor_plan_orbit_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_orbit_scale

case ('floor_plan_flip_label_side')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%floor_plan_flip_label_side

case ('floor_plan_size_is_absolute')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%floor_plan_size_is_absolute

case ('floor_plan_draw_only_first_pass')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%floor_plan_draw_only_first_pass

case ('correct_xy_distortion')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%correct_xy_distortion

case ('clip')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%clip

case ('draw_title')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_title

case ('draw_axes')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_axes

case ('draw_grid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_grid

case ('draw_curve_legend')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_curve_legend

case ('draw_only_good_user_data_or_vars')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_only_good_user_data_or_vars

case ('allow_wrap_around')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%allow_wrap_around

case ('floor_plan')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_floor_plan_struct (obj%floor_plan, rest, ptr, err, why)

case ('text_legend_origin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_point_struct (obj%text_legend_origin, rest, ptr, err, why)

case ('curve_legend_origin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_point_struct (obj%curve_legend_origin, rest, ptr, err, why)

case ('margin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%margin, rest, ptr, err, why)

case ('scale_margin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%scale_margin, rest, ptr, err, why)

case ('x')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%x, rest, ptr, err, why)

case ('y')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%y, rest, ptr, err, why)

case ('x2')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%x2, rest, ptr, err, why)

case ('y2')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%y2, rest, ptr, err, why)

case ('curve_legend')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_legend_struct (obj%curve_legend, rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GRAPH_INPUT)'
end select

end subroutine tao_res_tao_graph_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_graph_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_graph_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_graph_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_graph_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case (2);  ptr%str => obj%type
  case (3);  ptr%str => obj%title
  case (4);  ptr%str => obj%component
  case (5);  ptr%str => obj%text_legend(1)
  case (6);  ptr%str => obj%text_legend(2)
  case (7);  ptr%str => obj%text_legend(3)
  case (8);  ptr%str => obj%text_legend(4)
  case (9);  ptr%str => obj%text_legend(5)
  case (10);  ptr%str => obj%text_legend(6)
  case (11);  ptr%str => obj%text_legend(7)
  case (12);  ptr%str => obj%text_legend(8)
  case (13);  ptr%str => obj%text_legend(9)
  case (14);  ptr%str => obj%text_legend(10)
  case (15);  ptr%str => obj%floor_plan_view
  case (16);  ptr%str => obj%floor_plan_orbit_color
  case (17);  ptr%i => obj%box(1)
  case (18);  ptr%i => obj%box(2)
  case (19);  ptr%i => obj%box(3)
  case (20);  ptr%i => obj%box(4)
  case (21);  ptr%i => obj%ix_universe
  case (22);  ptr%i => obj%ix_branch
  case (23);  ptr%i => obj%n_curve
  case (24);  ptr%r => obj%x_axis_scale_factor
  case (25);  ptr%r => obj%symbol_size_scale
  case (26);  ptr%r => obj%floor_plan_rotation
  case (27);  ptr%r => obj%floor_plan_orbit_scale
  case (28);  ptr%l => obj%floor_plan_flip_label_side
  case (29);  ptr%l => obj%floor_plan_size_is_absolute
  case (30);  ptr%l => obj%floor_plan_draw_only_first_pass
  case (31);  ptr%l => obj%correct_xy_distortion
  case (32);  ptr%l => obj%clip
  case (33);  ptr%l => obj%draw_title
  case (34);  ptr%l => obj%draw_axes
  case (35);  ptr%l => obj%draw_grid
  case (36);  ptr%l => obj%draw_curve_legend
  case (37);  ptr%l => obj%draw_only_good_user_data_or_vars
  case (38);  ptr%l => obj%allow_wrap_around
  case (39);  ptr%str => obj%floor_plan%view
  case (40);  ptr%r => obj%floor_plan%rotation
  case (41);  ptr%l => obj%floor_plan%correct_distortion
  case (42);  ptr%l => obj%floor_plan%flip_label_side
  case (43);  ptr%l => obj%floor_plan%size_is_absolute
  case (44);  ptr%l => obj%floor_plan%draw_only_first_pass
  case (45);  ptr%l => obj%floor_plan%draw_building_wall
  case (46);  ptr%r => obj%floor_plan%orbit_scale
  case (47);  ptr%str => obj%floor_plan%orbit_color
  case (48);  ptr%str => obj%floor_plan%orbit_pattern
  case (49);  ptr%str => obj%floor_plan%orbit_lattice
  case (50);  ptr%i => obj%floor_plan%orbit_width
  case (51);  ptr%r => obj%text_legend_origin%x
  case (52);  ptr%r => obj%text_legend_origin%y
  case (53);  ptr%str => obj%text_legend_origin%units
  case (54);  ptr%r => obj%curve_legend_origin%x
  case (55);  ptr%r => obj%curve_legend_origin%y
  case (56);  ptr%str => obj%curve_legend_origin%units
  case (57);  ptr%r => obj%margin%x1
  case (58);  ptr%r => obj%margin%x2
  case (59);  ptr%r => obj%margin%y1
  case (60);  ptr%r => obj%margin%y2
  case (61);  ptr%str => obj%margin%units
  case (62);  ptr%r => obj%scale_margin%x1
  case (63);  ptr%r => obj%scale_margin%x2
  case (64);  ptr%r => obj%scale_margin%y1
  case (65);  ptr%r => obj%scale_margin%y2
  case (66);  ptr%str => obj%scale_margin%units
  case (67);  ptr%str => obj%x%label
  case (68);  ptr%r => obj%x%min
  case (69);  ptr%r => obj%x%max
  case (70);  ptr%r => obj%x%tick_min
  case (71);  ptr%r => obj%x%tick_max
  case (72);  ptr%r => obj%x%eval_min
  case (73);  ptr%r => obj%x%eval_max
  case (74);  ptr%r => obj%x%dtick
  case (75);  ptr%r => obj%x%number_offset
  case (76);  ptr%r => obj%x%label_offset
  case (77);  ptr%r => obj%x%major_tick_len
  case (78);  ptr%r => obj%x%minor_tick_len
  case (79);  ptr%str => obj%x%label_color
  case (80);  ptr%i => obj%x%major_div
  case (81);  ptr%i => obj%x%major_div_nominal
  case (82);  ptr%i => obj%x%minor_div
  case (83);  ptr%i => obj%x%minor_div_max
  case (84);  ptr%i => obj%x%places
  case (85);  ptr%str => obj%x%type
  case (86);  ptr%str => obj%x%bounds
  case (87);  ptr%i => obj%x%tick_side
  case (88);  ptr%i => obj%x%number_side
  case (89);  ptr%l => obj%x%draw_label
  case (90);  ptr%l => obj%x%draw_numbers
  case (91);  ptr%str => obj%y%label
  case (92);  ptr%r => obj%y%min
  case (93);  ptr%r => obj%y%max
  case (94);  ptr%r => obj%y%tick_min
  case (95);  ptr%r => obj%y%tick_max
  case (96);  ptr%r => obj%y%eval_min
  case (97);  ptr%r => obj%y%eval_max
  case (98);  ptr%r => obj%y%dtick
  case (99);  ptr%r => obj%y%number_offset
  case (100);  ptr%r => obj%y%label_offset
  case (101);  ptr%r => obj%y%major_tick_len
  case (102);  ptr%r => obj%y%minor_tick_len
  case (103);  ptr%str => obj%y%label_color
  case (104);  ptr%i => obj%y%major_div
  case (105);  ptr%i => obj%y%major_div_nominal
  case (106);  ptr%i => obj%y%minor_div
  case (107);  ptr%i => obj%y%minor_div_max
  case (108);  ptr%i => obj%y%places
  case (109);  ptr%str => obj%y%type
  case (110);  ptr%str => obj%y%bounds
  case (111);  ptr%i => obj%y%tick_side
  case (112);  ptr%i => obj%y%number_side
  case (113);  ptr%l => obj%y%draw_label
  case (114);  ptr%l => obj%y%draw_numbers
  case (115);  ptr%str => obj%x2%label
  case (116);  ptr%r => obj%x2%min
  case (117);  ptr%r => obj%x2%max
  case (118);  ptr%r => obj%x2%tick_min
  case (119);  ptr%r => obj%x2%tick_max
  case (120);  ptr%r => obj%x2%eval_min
  case (121);  ptr%r => obj%x2%eval_max
  case (122);  ptr%r => obj%x2%dtick
  case (123);  ptr%r => obj%x2%number_offset
  case (124);  ptr%r => obj%x2%label_offset
  case (125);  ptr%r => obj%x2%major_tick_len
  case (126);  ptr%r => obj%x2%minor_tick_len
  case (127);  ptr%str => obj%x2%label_color
  case (128);  ptr%i => obj%x2%major_div
  case (129);  ptr%i => obj%x2%major_div_nominal
  case (130);  ptr%i => obj%x2%minor_div
  case (131);  ptr%i => obj%x2%minor_div_max
  case (132);  ptr%i => obj%x2%places
  case (133);  ptr%str => obj%x2%type
  case (134);  ptr%str => obj%x2%bounds
  case (135);  ptr%i => obj%x2%tick_side
  case (136);  ptr%i => obj%x2%number_side
  case (137);  ptr%l => obj%x2%draw_label
  case (138);  ptr%l => obj%x2%draw_numbers
  case (139);  ptr%str => obj%y2%label
  case (140);  ptr%r => obj%y2%min
  case (141);  ptr%r => obj%y2%max
  case (142);  ptr%r => obj%y2%tick_min
  case (143);  ptr%r => obj%y2%tick_max
  case (144);  ptr%r => obj%y2%eval_min
  case (145);  ptr%r => obj%y2%eval_max
  case (146);  ptr%r => obj%y2%dtick
  case (147);  ptr%r => obj%y2%number_offset
  case (148);  ptr%r => obj%y2%label_offset
  case (149);  ptr%r => obj%y2%major_tick_len
  case (150);  ptr%r => obj%y2%minor_tick_len
  case (151);  ptr%str => obj%y2%label_color
  case (152);  ptr%i => obj%y2%major_div
  case (153);  ptr%i => obj%y2%major_div_nominal
  case (154);  ptr%i => obj%y2%minor_div
  case (155);  ptr%i => obj%y2%minor_div_max
  case (156);  ptr%i => obj%y2%places
  case (157);  ptr%str => obj%y2%type
  case (158);  ptr%str => obj%y2%bounds
  case (159);  ptr%i => obj%y2%tick_side
  case (160);  ptr%i => obj%y2%number_side
  case (161);  ptr%l => obj%y2%draw_label
  case (162);  ptr%l => obj%y2%draw_numbers
  case (163);  ptr%r => obj%curve_legend%row_spacing
  case (164);  ptr%r => obj%curve_legend%line_length
  case (165);  ptr%r => obj%curve_legend%text_offset
  case (166);  ptr%l => obj%curve_legend%draw_line
  case (167);  ptr%l => obj%curve_legend%draw_symbol
  case (168);  ptr%l => obj%curve_legend%draw_text
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_GRAPH_INPUT HAS 168 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%type

case ('title')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%title

case ('component')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%component

case ('text_legend')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 10 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%text_legend(1 + i_slot - 1)

case ('floor_plan_view')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%floor_plan_view

case ('floor_plan_orbit_color')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%floor_plan_orbit_color

case ('box')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 4 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%box(1 + i_slot - 1)

case ('ix_universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_universe

case ('ix_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_branch

case ('n_curve')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_curve

case ('x_axis_scale_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_axis_scale_factor

case ('symbol_size_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%symbol_size_scale

case ('floor_plan_rotation')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_rotation

case ('floor_plan_orbit_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_orbit_scale

case ('floor_plan_flip_label_side')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%floor_plan_flip_label_side

case ('floor_plan_size_is_absolute')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%floor_plan_size_is_absolute

case ('floor_plan_draw_only_first_pass')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%floor_plan_draw_only_first_pass

case ('correct_xy_distortion')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%correct_xy_distortion

case ('clip')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%clip

case ('draw_title')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_title

case ('draw_axes')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_axes

case ('draw_grid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_grid

case ('draw_curve_legend')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_curve_legend

case ('draw_only_good_user_data_or_vars')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_only_good_user_data_or_vars

case ('allow_wrap_around')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%allow_wrap_around

case ('floor_plan')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_floor_plan_struct_slot (obj%floor_plan, rest, i_slot, ptr, err, why)

case ('text_legend_origin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_point_struct_slot (obj%text_legend_origin, rest, i_slot, ptr, err, why)

case ('curve_legend_origin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_point_struct_slot (obj%curve_legend_origin, rest, i_slot, ptr, err, why)

case ('margin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%margin, rest, i_slot, ptr, err, why)

case ('scale_margin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%scale_margin, rest, i_slot, ptr, err, why)

case ('x')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%x, rest, i_slot, ptr, err, why)

case ('y')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%y, rest, i_slot, ptr, err, why)

case ('x2')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%x2, rest, i_slot, ptr, err, why)

case ('y2')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%y2, rest, i_slot, ptr, err, why)

case ('curve_legend')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_legend_struct_slot (obj%curve_legend, rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GRAPH_INPUT)'
end select

end subroutine tao_res_tao_graph_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_graph_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_graph_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_graph_struct (obj, name, ptr, err, why)

type (tao_graph_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%type

case ('title')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%title

case ('title_suffix')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%title_suffix

case ('text_legend')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 10, err, why)) return
  if (has_sub) then
    ptr%str => obj%text_legend(isub)
  else
    ptr%str1 => obj%text_legend
  endif

case ('text_legend_out')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 10, err, why)) return
  if (has_sub) then
    ptr%str => obj%text_legend_out(isub)
  else
    ptr%str1 => obj%text_legend_out
  endif

case ('why_invalid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%why_invalid

case ('curve')
  if (tao_res_alloc_bad(head, allocated(obj%curve), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%curve,1), ubound(obj%curve,1), err, why)) return
  call tao_res_tao_curve_struct (obj%curve(isub), rest, ptr, err, why)

case ('p')
  err = .true.; why = 'COMPONENT P ' // &
          'IS A POINTER COMPONENT AND CANNOT BE SET'

case ('floor_plan')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_floor_plan_struct (obj%floor_plan, rest, ptr, err, why)

case ('text_legend_origin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_point_struct (obj%text_legend_origin, rest, ptr, err, why)

case ('curve_legend_origin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_point_struct (obj%curve_legend_origin, rest, ptr, err, why)

case ('curve_legend')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_legend_struct (obj%curve_legend, rest, ptr, err, why)

case ('x')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%x, rest, ptr, err, why)

case ('y')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%y, rest, ptr, err, why)

case ('x2')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%x2, rest, ptr, err, why)

case ('y2')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%y2, rest, ptr, err, why)

case ('margin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%margin, rest, ptr, err, why)

case ('scale_margin')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%scale_margin, rest, ptr, err, why)

case ('x_axis_scale_factor')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x_axis_scale_factor

case ('symbol_size_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%symbol_size_scale

case ('box')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 4, err, why)) return
  if (has_sub) then
    ptr%i => obj%box(isub)
  else
    ptr%i1 => obj%box
  endif

case ('ix_branch')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_branch

case ('ix_universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_universe

case ('clip')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%clip

case ('y2_mirrors_y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%y2_mirrors_y

case ('limited')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%limited

case ('draw_axes')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_axes

case ('draw_curve_legend')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_curve_legend

case ('draw_grid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_grid

case ('draw_title')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_title

case ('draw_only_good_user_data_or_vars')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_only_good_user_data_or_vars

case ('allow_wrap_around')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%allow_wrap_around

case ('is_valid')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%is_valid

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GRAPH_STRUCT)'
end select

end subroutine tao_res_tao_graph_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_graph_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_graph_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_graph_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_graph_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_GRAPH_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%type

case ('title')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%title

case ('title_suffix')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%title_suffix

case ('text_legend')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 10 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%text_legend(1 + i_slot - 1)

case ('text_legend_out')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 10 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%text_legend_out(1 + i_slot - 1)

case ('why_invalid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%why_invalid

case ('curve')
  if (tao_res_alloc_bad(head, allocated(obj%curve), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%curve,1), ubound(obj%curve,1), err, why)) return
  call tao_res_tao_curve_struct_slot (obj%curve(isub), rest, i_slot, ptr, err, why)

case ('floor_plan')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_floor_plan_struct_slot (obj%floor_plan, rest, i_slot, ptr, err, why)

case ('text_legend_origin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_point_struct_slot (obj%text_legend_origin, rest, i_slot, ptr, err, why)

case ('curve_legend_origin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_point_struct_slot (obj%curve_legend_origin, rest, i_slot, ptr, err, why)

case ('curve_legend')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_legend_struct_slot (obj%curve_legend, rest, i_slot, ptr, err, why)

case ('x')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%x, rest, i_slot, ptr, err, why)

case ('y')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%y, rest, i_slot, ptr, err, why)

case ('x2')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%x2, rest, i_slot, ptr, err, why)

case ('y2')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%y2, rest, i_slot, ptr, err, why)

case ('margin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%margin, rest, i_slot, ptr, err, why)

case ('scale_margin')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%scale_margin, rest, i_slot, ptr, err, why)

case ('x_axis_scale_factor')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x_axis_scale_factor

case ('symbol_size_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%symbol_size_scale

case ('box')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 4 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%box(1 + i_slot - 1)

case ('ix_branch')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_branch

case ('ix_universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_universe

case ('clip')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%clip

case ('y2_mirrors_y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%y2_mirrors_y

case ('limited')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%limited

case ('draw_axes')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_axes

case ('draw_curve_legend')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_curve_legend

case ('draw_grid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_grid

case ('draw_title')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_title

case ('draw_only_good_user_data_or_vars')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_only_good_user_data_or_vars

case ('allow_wrap_around')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%allow_wrap_around

case ('is_valid')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%is_valid

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_GRAPH_STRUCT)'
end select

end subroutine tao_res_tao_graph_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_histogram_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_histogram_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_histogram_struct (obj, name, ptr, err, why)

type (tao_histogram_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('density_normalized')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%density_normalized

case ('weight_by_charge')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%weight_by_charge

case ('minimum')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%minimum

case ('maximum')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%maximum

case ('width')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%width

case ('center')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%center

case ('number')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%number

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_HISTOGRAM_STRUCT)'
end select

end subroutine tao_res_tao_histogram_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_histogram_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_histogram_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_histogram_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_histogram_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%l => obj%density_normalized
  case (2);  ptr%l => obj%weight_by_charge
  case (3);  ptr%r => obj%minimum
  case (4);  ptr%r => obj%maximum
  case (5);  ptr%r => obj%width
  case (6);  ptr%r => obj%center
  case (7);  ptr%i => obj%number
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_HISTOGRAM_STRUCT HAS 7 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('density_normalized')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%density_normalized

case ('weight_by_charge')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%weight_by_charge

case ('minimum')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%minimum

case ('maximum')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%maximum

case ('width')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%width

case ('center')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%center

case ('number')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%number

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_HISTOGRAM_STRUCT)'
end select

end subroutine tao_res_tao_histogram_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_key_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_key_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_key_input (obj, name, ptr, err, why)

type (tao_key_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_name

case ('attrib_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%attrib_name

case ('delta')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%delta

case ('universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%universe

case ('small_step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%small_step

case ('low_lim')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%low_lim

case ('high_lim')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%high_lim

case ('weight')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%weight

case ('good_opt')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%good_opt

case ('merit_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%merit_type

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_KEY_INPUT)'
end select

end subroutine tao_res_tao_key_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_key_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_key_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_key_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_key_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%ele_name
  case (2);  ptr%str => obj%attrib_name
  case (3);  ptr%r => obj%delta
  case (4);  ptr%str => obj%universe
  case (5);  ptr%r => obj%small_step
  case (6);  ptr%r => obj%low_lim
  case (7);  ptr%r => obj%high_lim
  case (8);  ptr%r => obj%weight
  case (9);  ptr%l => obj%good_opt
  case (10);  ptr%str => obj%merit_type
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_KEY_INPUT HAS 10 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_name

case ('attrib_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%attrib_name

case ('delta')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%delta

case ('universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%universe

case ('small_step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%small_step

case ('low_lim')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%low_lim

case ('high_lim')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%high_lim

case ('weight')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%weight

case ('good_opt')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%good_opt

case ('merit_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%merit_type

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_KEY_INPUT)'
end select

end subroutine tao_res_tao_key_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_place_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_place_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_place_input (obj, name, ptr, err, why)

type (tao_place_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('region')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%region

case ('plot')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%plot

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLACE_INPUT)'
end select

end subroutine tao_res_tao_place_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_place_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_place_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_place_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_place_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%region
  case (2);  ptr%str => obj%plot
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_PLACE_INPUT HAS 2 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('region')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%region

case ('plot')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%plot

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLACE_INPUT)'
end select

end subroutine tao_res_tao_place_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_plot_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_plot_input (obj, name, ptr, err, why)

type (tao_plot_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('description')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%description

case ('x_axis_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%x_axis_type

case ('n_graph')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_graph

case ('autoscale_gang_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_gang_x

case ('autoscale_gang_y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_gang_y

case ('autoscale_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_x

case ('autoscale_y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_y

case ('n_curve_pts')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_curve_pts

case ('x')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_axis_struct (obj%x, rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_INPUT)'
end select

end subroutine tao_res_tao_plot_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_plot_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_plot_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_plot_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case (2);  ptr%str => obj%description
  case (3);  ptr%str => obj%x_axis_type
  case (4);  ptr%i => obj%n_graph
  case (5);  ptr%l => obj%autoscale_gang_x
  case (6);  ptr%l => obj%autoscale_gang_y
  case (7);  ptr%l => obj%autoscale_x
  case (8);  ptr%l => obj%autoscale_y
  case (9);  ptr%i => obj%n_curve_pts
  case (10);  ptr%str => obj%x%label
  case (11);  ptr%r => obj%x%min
  case (12);  ptr%r => obj%x%max
  case (13);  ptr%r => obj%x%tick_min
  case (14);  ptr%r => obj%x%tick_max
  case (15);  ptr%r => obj%x%eval_min
  case (16);  ptr%r => obj%x%eval_max
  case (17);  ptr%r => obj%x%dtick
  case (18);  ptr%r => obj%x%number_offset
  case (19);  ptr%r => obj%x%label_offset
  case (20);  ptr%r => obj%x%major_tick_len
  case (21);  ptr%r => obj%x%minor_tick_len
  case (22);  ptr%str => obj%x%label_color
  case (23);  ptr%i => obj%x%major_div
  case (24);  ptr%i => obj%x%major_div_nominal
  case (25);  ptr%i => obj%x%minor_div
  case (26);  ptr%i => obj%x%minor_div_max
  case (27);  ptr%i => obj%x%places
  case (28);  ptr%str => obj%x%type
  case (29);  ptr%str => obj%x%bounds
  case (30);  ptr%i => obj%x%tick_side
  case (31);  ptr%i => obj%x%number_side
  case (32);  ptr%l => obj%x%draw_label
  case (33);  ptr%l => obj%x%draw_numbers
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_PLOT_INPUT HAS 33 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('description')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%description

case ('x_axis_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%x_axis_type

case ('n_graph')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_graph

case ('autoscale_gang_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_gang_x

case ('autoscale_gang_y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_gang_y

case ('autoscale_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_x

case ('autoscale_y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_y

case ('n_curve_pts')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_curve_pts

case ('x')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_axis_struct_slot (obj%x, rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_INPUT)'
end select

end subroutine tao_res_tao_plot_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_page_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_plot_page_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_plot_page_input (obj, name, ptr, err, why)

type (tao_plot_page_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('title')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_title_struct (obj%title, rest, ptr, err, why)

case ('subtitle')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_title_struct (obj%subtitle, rest, ptr, err, why)

case ('border')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%border, rest, ptr, err, why)

case ('plot_display_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%plot_display_type

case ('size')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 2, err, why)) return
  if (has_sub) then
    ptr%r => obj%size(isub)
  else
    ptr%r1 => obj%size
  endif

case ('text_height')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%text_height

case ('main_title_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%main_title_text_scale

case ('graph_title_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%graph_title_text_scale

case ('axis_number_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%axis_number_text_scale

case ('axis_label_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%axis_label_text_scale

case ('legend_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%legend_text_scale

case ('key_table_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%key_table_text_scale

case ('floor_plan_shape_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_shape_scale

case ('floor_plan_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_text_scale

case ('lat_layout_shape_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lat_layout_shape_scale

case ('lat_layout_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lat_layout_text_scale

case ('curve_legend_line_len')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%curve_legend_line_len

case ('curve_legend_text_offset')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%curve_legend_text_offset

case ('n_curve_pts')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_curve_pts

case ('delete_overlapping_plots')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%delete_overlapping_plots

case ('draw_graph_title_suffix')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_graph_title_suffix

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_PAGE_INPUT)'
end select

end subroutine tao_res_tao_plot_page_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_page_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_plot_page_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_plot_page_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_plot_page_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%title%string
  case (2);  ptr%r => obj%title%x
  case (3);  ptr%r => obj%title%y
  case (4);  ptr%str => obj%title%units
  case (5);  ptr%str => obj%title%justify
  case (6);  ptr%l => obj%title%draw_it
  case (7);  ptr%str => obj%subtitle%string
  case (8);  ptr%r => obj%subtitle%x
  case (9);  ptr%r => obj%subtitle%y
  case (10);  ptr%str => obj%subtitle%units
  case (11);  ptr%str => obj%subtitle%justify
  case (12);  ptr%l => obj%subtitle%draw_it
  case (13);  ptr%r => obj%border%x1
  case (14);  ptr%r => obj%border%x2
  case (15);  ptr%r => obj%border%y1
  case (16);  ptr%r => obj%border%y2
  case (17);  ptr%str => obj%border%units
  case (18);  ptr%str => obj%plot_display_type
  case (19);  ptr%r => obj%size(1)
  case (20);  ptr%r => obj%size(2)
  case (21);  ptr%r => obj%text_height
  case (22);  ptr%r => obj%main_title_text_scale
  case (23);  ptr%r => obj%graph_title_text_scale
  case (24);  ptr%r => obj%axis_number_text_scale
  case (25);  ptr%r => obj%axis_label_text_scale
  case (26);  ptr%r => obj%legend_text_scale
  case (27);  ptr%r => obj%key_table_text_scale
  case (28);  ptr%r => obj%floor_plan_shape_scale
  case (29);  ptr%r => obj%floor_plan_text_scale
  case (30);  ptr%r => obj%lat_layout_shape_scale
  case (31);  ptr%r => obj%lat_layout_text_scale
  case (32);  ptr%r => obj%curve_legend_line_len
  case (33);  ptr%r => obj%curve_legend_text_offset
  case (34);  ptr%i => obj%n_curve_pts
  case (35);  ptr%l => obj%delete_overlapping_plots
  case (36);  ptr%l => obj%draw_graph_title_suffix
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_PLOT_PAGE_INPUT HAS 36 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('title')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_title_struct_slot (obj%title, rest, i_slot, ptr, err, why)

case ('subtitle')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_title_struct_slot (obj%subtitle, rest, i_slot, ptr, err, why)

case ('border')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%border, rest, i_slot, ptr, err, why)

case ('plot_display_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%plot_display_type

case ('size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 2 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%size(1 + i_slot - 1)

case ('text_height')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%text_height

case ('main_title_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%main_title_text_scale

case ('graph_title_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%graph_title_text_scale

case ('axis_number_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%axis_number_text_scale

case ('axis_label_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%axis_label_text_scale

case ('legend_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%legend_text_scale

case ('key_table_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%key_table_text_scale

case ('floor_plan_shape_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_shape_scale

case ('floor_plan_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_text_scale

case ('lat_layout_shape_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lat_layout_shape_scale

case ('lat_layout_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lat_layout_text_scale

case ('curve_legend_line_len')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%curve_legend_line_len

case ('curve_legend_text_offset')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%curve_legend_text_offset

case ('n_curve_pts')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_curve_pts

case ('delete_overlapping_plots')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%delete_overlapping_plots

case ('draw_graph_title_suffix')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_graph_title_suffix

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_PAGE_INPUT)'
end select

end subroutine tao_res_tao_plot_page_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_page_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_plot_page_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_plot_page_struct (obj, name, ptr, err, why)

type (tao_plot_page_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('title')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_title_struct (obj%title, rest, ptr, err, why)

case ('subtitle')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_title_struct (obj%subtitle, rest, ptr, err, why)

case ('border')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_rect_struct (obj%border, rest, ptr, err, why)

case ('floor_plan')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_drawing_struct (obj%floor_plan, rest, ptr, err, why)

case ('lat_layout')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_drawing_struct (obj%lat_layout, rest, ptr, err, why)

case ('pattern')
  if (tao_res_alloc_bad(head, allocated(obj%pattern), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%pattern,1), ubound(obj%pattern,1), err, why)) return
  call tao_res_tao_shape_pattern_struct (obj%pattern(isub), rest, ptr, err, why)

case ('template')
  if (tao_res_alloc_bad(head, allocated(obj%template), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%template,1), ubound(obj%template,1), err, why)) return
  call tao_res_tao_plot_struct (obj%template(isub), rest, ptr, err, why)

case ('region')
  if (tao_res_alloc_bad(head, allocated(obj%region), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%region,1), ubound(obj%region,1), err, why)) return
  call tao_res_tao_plot_region_struct (obj%region(isub), rest, ptr, err, why)

case ('plot_display_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%plot_display_type

case ('size')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 2, err, why)) return
  if (has_sub) then
    ptr%r => obj%size(isub)
  else
    ptr%r1 => obj%size
  endif

case ('text_height')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%text_height

case ('main_title_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%main_title_text_scale

case ('graph_title_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%graph_title_text_scale

case ('axis_number_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%axis_number_text_scale

case ('axis_label_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%axis_label_text_scale

case ('legend_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%legend_text_scale

case ('key_table_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%key_table_text_scale

case ('floor_plan_shape_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_shape_scale

case ('floor_plan_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%floor_plan_text_scale

case ('lat_layout_shape_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lat_layout_shape_scale

case ('lat_layout_text_scale')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%lat_layout_text_scale

case ('n_curve_pts')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_curve_pts

case ('id_window')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%id_window

case ('delete_overlapping_plots')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%delete_overlapping_plots

case ('draw_graph_title_suffix')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_graph_title_suffix

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_PAGE_STRUCT)'
end select

end subroutine tao_res_tao_plot_page_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_page_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_plot_page_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_plot_page_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_plot_page_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_PLOT_PAGE_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('title')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_title_struct_slot (obj%title, rest, i_slot, ptr, err, why)

case ('subtitle')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_title_struct_slot (obj%subtitle, rest, i_slot, ptr, err, why)

case ('border')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_rect_struct_slot (obj%border, rest, i_slot, ptr, err, why)

case ('floor_plan')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_drawing_struct_slot (obj%floor_plan, rest, i_slot, ptr, err, why)

case ('lat_layout')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_drawing_struct_slot (obj%lat_layout, rest, i_slot, ptr, err, why)

case ('pattern')
  if (tao_res_alloc_bad(head, allocated(obj%pattern), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%pattern,1), ubound(obj%pattern,1), err, why)) return
  call tao_res_tao_shape_pattern_struct_slot (obj%pattern(isub), rest, i_slot, ptr, err, why)

case ('template')
  if (tao_res_alloc_bad(head, allocated(obj%template), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%template,1), ubound(obj%template,1), err, why)) return
  call tao_res_tao_plot_struct_slot (obj%template(isub), rest, i_slot, ptr, err, why)

case ('region')
  if (tao_res_alloc_bad(head, allocated(obj%region), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%region,1), ubound(obj%region,1), err, why)) return
  call tao_res_tao_plot_region_struct_slot (obj%region(isub), rest, i_slot, ptr, err, why)

case ('plot_display_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%plot_display_type

case ('size')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 2 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%size(1 + i_slot - 1)

case ('text_height')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%text_height

case ('main_title_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%main_title_text_scale

case ('graph_title_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%graph_title_text_scale

case ('axis_number_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%axis_number_text_scale

case ('axis_label_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%axis_label_text_scale

case ('legend_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%legend_text_scale

case ('key_table_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%key_table_text_scale

case ('floor_plan_shape_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_shape_scale

case ('floor_plan_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%floor_plan_text_scale

case ('lat_layout_shape_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lat_layout_shape_scale

case ('lat_layout_text_scale')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%lat_layout_text_scale

case ('n_curve_pts')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_curve_pts

case ('id_window')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%id_window

case ('delete_overlapping_plots')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%delete_overlapping_plots

case ('draw_graph_title_suffix')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_graph_title_suffix

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_PAGE_STRUCT)'
end select

end subroutine tao_res_tao_plot_page_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_region_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_plot_region_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_plot_region_struct (obj, name, ptr, err, why)

type (tao_plot_region_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('plot')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_tao_plot_struct (obj%plot, rest, ptr, err, why)

case ('location')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 4, err, why)) return
  if (has_sub) then
    ptr%r => obj%location(isub)
  else
    ptr%r1 => obj%location
  endif

case ('visible')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%visible

case ('list_with_show_plot_command')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%list_with_show_plot_command

case ('setup_done')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%setup_done

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_REGION_STRUCT)'
end select

end subroutine tao_res_tao_plot_region_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_region_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_plot_region_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_plot_region_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_plot_region_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_PLOT_REGION_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('plot')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_tao_plot_struct_slot (obj%plot, rest, i_slot, ptr, err, why)

case ('location')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 4 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%location(1 + i_slot - 1)

case ('visible')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%visible

case ('list_with_show_plot_command')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%list_with_show_plot_command

case ('setup_done')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%setup_done

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_REGION_STRUCT)'
end select

end subroutine tao_res_tao_plot_region_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_plot_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_plot_struct (obj, name, ptr, err, why)

type (tao_plot_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('description')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%description

case ('graph')
  if (tao_res_alloc_bad(head, allocated(obj%graph), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%graph,1), ubound(obj%graph,1), err, why)) return
  call tao_res_tao_graph_struct (obj%graph(isub), rest, ptr, err, why)

case ('r')
  err = .true.; why = 'COMPONENT R ' // &
          'IS A POINTER COMPONENT AND CANNOT BE SET'

case ('ix_plot')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%ix_plot

case ('n_curve_pts')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%i => obj%n_curve_pts

case ('type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%type

case ('x_axis_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%x_axis_type

case ('autoscale_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_x

case ('autoscale_y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_y

case ('autoscale_gang_x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_gang_x

case ('autoscale_gang_y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%autoscale_gang_y

case ('list_with_show_plot_command')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%list_with_show_plot_command

case ('phantom')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%phantom

case ('default_plot')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%default_plot

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_STRUCT)'
end select

end subroutine tao_res_tao_plot_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_plot_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_plot_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_plot_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_plot_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_PLOT_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('description')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%description

case ('graph')
  if (tao_res_alloc_bad(head, allocated(obj%graph), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%graph,1), ubound(obj%graph,1), err, why)) return
  call tao_res_tao_graph_struct_slot (obj%graph(isub), rest, i_slot, ptr, err, why)

case ('ix_plot')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%ix_plot

case ('n_curve_pts')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%i => obj%n_curve_pts

case ('type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%type

case ('x_axis_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%x_axis_type

case ('autoscale_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_x

case ('autoscale_y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_y

case ('autoscale_gang_x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_gang_x

case ('autoscale_gang_y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%autoscale_gang_y

case ('list_with_show_plot_command')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%list_with_show_plot_command

case ('phantom')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%phantom

case ('default_plot')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%default_plot

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_PLOT_STRUCT)'
end select

end subroutine tao_res_tao_plot_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_region_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_region_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_region_input (obj, name, ptr, err, why)

type (tao_region_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('location')
  if (tao_res_array_bad(head, rest, has_sub, isub, 1, 4, err, why)) return
  if (has_sub) then
    ptr%r => obj%location(isub)
  else
    ptr%r1 => obj%location
  endif

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_REGION_INPUT)'
end select

end subroutine tao_res_tao_region_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_region_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_region_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_region_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_region_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case (2);  ptr%r => obj%location(1)
  case (3);  ptr%r => obj%location(2)
  case (4);  ptr%r => obj%location(3)
  case (5);  ptr%r => obj%location(4)
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_REGION_INPUT HAS 5 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('location')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot < 1 .or. i_slot > 4 - 1 + 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%location(1 + i_slot - 1)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_REGION_INPUT)'
end select

end subroutine tao_res_tao_region_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_shape_pattern_point_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_shape_pattern_point_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_shape_pattern_point_struct (obj, name, ptr, err, why)

type (tao_shape_pattern_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('s')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%s

case ('y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y

case ('radius')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%radius

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_SHAPE_PATTERN_POINT_STRUCT)'
end select

end subroutine tao_res_tao_shape_pattern_point_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_shape_pattern_point_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_shape_pattern_point_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_shape_pattern_point_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_shape_pattern_point_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%r => obj%s
  case (2);  ptr%r => obj%y
  case (3);  ptr%r => obj%radius
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_SHAPE_PATTERN_POINT_STRUCT HAS 3 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('s')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%s

case ('y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y

case ('radius')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%radius

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_SHAPE_PATTERN_POINT_STRUCT)'
end select

end subroutine tao_res_tao_shape_pattern_point_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_shape_pattern_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_shape_pattern_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_shape_pattern_struct (obj, name, ptr, err, why)

type (tao_shape_pattern_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case ('line')
  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return
  call tao_res_qp_line_struct (obj%line, rest, ptr, err, why)

case ('pt')
  if (tao_res_alloc_bad(head, allocated(obj%pt), err, why)) return
  if (tao_res_struct_array_bad(head, rest, has_sub, isub, lbound(obj%pt,1), ubound(obj%pt,1), err, why)) return
  call tao_res_tao_shape_pattern_point_struct (obj%pt(isub), rest, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_SHAPE_PATTERN_STRUCT)'
end select

end subroutine tao_res_tao_shape_pattern_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_shape_pattern_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_shape_pattern_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_shape_pattern_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_shape_pattern_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  ! This structure has an allocatable or deferred shape component so it has no
  ! fixed component ordering and cannot be assigned positionally.
  err = .true.
  why = 'STRUCTURE TAO_SHAPE_PATTERN_STRUCT CANNOT BE SET FROM A POSITIONAL VALUE LIST'
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case ('line')
  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return
  call tao_res_qp_line_struct_slot (obj%line, rest, i_slot, ptr, err, why)

case ('pt')
  if (tao_res_alloc_bad(head, allocated(obj%pt), err, why)) return
  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, lbound(obj%pt,1), ubound(obj%pt,1), err, why)) return
  call tao_res_tao_shape_pattern_point_struct_slot (obj%pt(isub), rest, i_slot, ptr, err, why)

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_SHAPE_PATTERN_STRUCT)'
end select

end subroutine tao_res_tao_shape_pattern_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_title_struct (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_title_struct instance to a pointer.
! Structure defined in: tao/code/tao_struct.f90
!-

recursive subroutine tao_res_tao_title_struct (obj, name, ptr, err, why)

type (tao_title_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('string')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%string

case ('x')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%x

case ('y')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%y

case ('units')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%units

case ('justify')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%justify

case ('draw_it')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%draw_it

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_TITLE_STRUCT)'
end select

end subroutine tao_res_tao_title_struct

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_title_struct_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_title_struct.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_title_struct_slot (obj, name, i_slot, ptr, err, why)

type (tao_title_struct), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%string
  case (2);  ptr%r => obj%x
  case (3);  ptr%r => obj%y
  case (4);  ptr%str => obj%units
  case (5);  ptr%str => obj%justify
  case (6);  ptr%l => obj%draw_it
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_TITLE_STRUCT HAS 6 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('string')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%string

case ('x')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%x

case ('y')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%y

case ('units')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%units

case ('justify')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%justify

case ('draw_it')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%draw_it

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_TITLE_STRUCT)'
end select

end subroutine tao_res_tao_title_struct_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_v1_var_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_v1_var_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_v1_var_input (obj, name, ptr, err, why)

type (tao_v1_var_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_V1_VAR_INPUT)'
end select

end subroutine tao_res_tao_v1_var_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_v1_var_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_v1_var_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_v1_var_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_v1_var_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%name
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_V1_VAR_INPUT HAS 1 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%name

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_V1_VAR_INPUT)'
end select

end subroutine tao_res_tao_v1_var_input_slot

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_var_input (obj, name, ptr, err, why)
!
! Resolve a component name of a tao_var_input instance to a pointer.
! Structure defined in: tao/code/tao_input_struct.f90
!-

recursive subroutine tao_res_tao_var_input (obj, name, ptr, err, why)

type (tao_var_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer isub

logical, intent(out) :: err
logical has_sub

!

! ptr is cleared on entry. Without this a component left associated by an
! earlier call would win over the one set here, since tao_set_ptr_value takes
! the first associated component it finds.

ptr = tao_ptr_struct()

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_name')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%ele_name

case ('attribute')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%attribute

case ('universe')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%universe

case ('weight')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%weight

case ('step')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%step

case ('low_lim')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%low_lim

case ('high_lim')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%high_lim

case ('merit_type')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%str => obj%merit_type

case ('good_user')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%good_user

case ('key_bound')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%l => obj%key_bound

case ('key_delta')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%key_delta

case ('meas')
  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return
  ptr%r => obj%meas

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_VAR_INPUT)'
end select

end subroutine tao_res_tao_var_input

!--------------------------------------------------------------------------
!+
! Subroutine tao_res_tao_var_input_slot (obj, name, i_slot, ptr, err, why)
!
! Pointer to the i_slot-th ultimate component, in declaration order, of the thing
! that name selects within a tao_var_input.
!
! This is the order a namelist read uses when something is given a positional value
! list. A blank name means the whole structure, as in
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! and a non blank name selects a component, which may itself be a structure, as in
!     plot_page%border = 0, 0, 0, 0, "%PAGE"
!
! err is set when i_slot is past the last component, which is how the caller
! detects that too many values were supplied.
!-

recursive subroutine tao_res_tao_var_input_slot (obj, name, i_slot, ptr, err, why)

type (tao_var_input), target :: obj
type (tao_ptr_struct) ptr

character(*), intent(in) :: name
character(*), intent(out) :: why
character(:), allocatable :: head, rest

integer, intent(in) :: i_slot
integer isub

logical, intent(out) :: err
logical has_sub

!

ptr = tao_ptr_struct()
err = .false.
why = ''

! A blank name means the whole structure.

if (name == '') then
  select case (i_slot)
  case (1);  ptr%str => obj%ele_name
  case (2);  ptr%str => obj%attribute
  case (3);  ptr%str => obj%universe
  case (4);  ptr%r => obj%weight
  case (5);  ptr%r => obj%step
  case (6);  ptr%r => obj%low_lim
  case (7);  ptr%r => obj%high_lim
  case (8);  ptr%str => obj%merit_type
  case (9);  ptr%l => obj%good_user
  case (10);  ptr%l => obj%key_bound
  case (11);  ptr%r => obj%key_delta
  case (12);  ptr%r => obj%meas
  case default
    err = .true.
    why = 'TOO MANY VALUES. TAO_VAR_INPUT HAS 12 COMPONENTS'
  end select
  return
endif

! A non blank name selects a component to descend into.

call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
if (err) return

select case (head)

case ('ele_name')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%ele_name

case ('attribute')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%attribute

case ('universe')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%universe

case ('weight')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%weight

case ('step')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%step

case ('low_lim')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%low_lim

case ('high_lim')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%high_lim

case ('merit_type')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%str => obj%merit_type

case ('good_user')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%good_user

case ('key_bound')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%l => obj%key_bound

case ('key_delta')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%key_delta

case ('meas')
  if (tao_res_slot_leaf_bad(head, rest, err, why)) return
  if (i_slot /= 1) then
    err = .true.
    why = 'TOO MANY VALUES FOR COMPONENT: ' // head
    return
  endif
  ptr%r => obj%meas

case default
  err = .true.
  why = 'NO SUCH COMPONENT: ' // head // '  (IN TAO_VAR_INPUT)'
end select

end subroutine tao_res_tao_var_input_slot

end module tao_attrib_resolve_mod
