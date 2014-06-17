module tao_show_mod

use tao_mod
use tao_data_and_eval_mod
use tao_top10_mod
use tao_command_mod, only: tao_next_switch
use tao_lattice_calc_mod

type show_common_struct
  type (ele_struct), pointer :: ele 
  type (coord_struct), pointer :: orbit 
  type (bunch_params_struct), pointer :: bunch_params
  type (tao_universe_struct), pointer :: u
  integer ix_ele
end type

private write_real

contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine tao_show_cmd (what, stuff)
!
! Show information on variable, parameters, elements, etc.
!
! Input:
!   what  -- Character(*): What to show.
!   stuff -- Character(*): Particular stuff to show.
!-

subroutine tao_show_cmd (what, stuff)

implicit none


integer iu, ix, n, nl
integer :: n_write_file = 0            ! used for indexing 'show write' files

character(*) what, stuff
character(100) file_name, result_id
character(len(what)) what2
character(len(stuff)) stuff2
character(n_char_show), allocatable, save :: lines(:)
character(16) :: r_name = 'tao_show_cmd'
character(20) switch

logical opened, err, doprint

! Init

what2 = what
stuff2 = stuff
opened = .false.
doprint = s%com%print_to_terminal

! See if the results need to be written to a file.

do
  call tao_next_switch (what2, ['-append ', '-write  ', '-noprint'], switch, err, ix)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-noprint')
    doprint = .false.

  case ('-append', '-write')
    call string_trim(stuff2, stuff2, ix)
    file_name = stuff2(:ix)
    call string_trim(stuff2(ix+1:), stuff2, ix)

    ix = index(file_name, '*')
    if (ix /= 0) then
      n_write_file = n_write_file + 1
      write (file_name, '(a, i3.3, a)') file_name(1:ix-1), n_write_file, trim(file_name(ix+1:))
    endif

    iu = lunget()
    if (switch == '-append') then
      open (iu, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = n_char_show)
    else
      open (iu, file = file_name, status = 'REPLACE', recl = n_char_show)
    endif

    opened = .true.
  end select

  call string_trim (stuff2, stuff2, ix)
  what2 = stuff2(1:ix)
  stuff2 = stuff2(ix+1:)

end do

! Get restults.
! Result_id is for tao_show_this to show exactly what it did.
! This info can be helpful in tao_hook_show_cmd.

if (opened) call output_direct (iu)  ! tell out_io to write to a file

call tao_show_this (what2, stuff2, result_id, lines, nl)  
call tao_hook_show_cmd (what2, stuff2, result_id, lines, nl)

if (nl > 0) then
  if (result_id == 'ERROR') then
    call out_io (s_error$, r_name, lines(1:nl))
  else
    call output_direct (do_print = doprint)
    call out_io (s_blank$, r_name, lines(1:nl))
  endif
endif

! Finish

call output_direct (0, do_print=s%com%print_to_terminal)  ! reset to not write to a file

if (opened) then
  close (iu)
  call out_io (s_blank$, r_name, '', 'Written to file: ' // file_name)
endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine tao_show_this (what, stuff, result_id, lines, nl)

use random_mod
use csr_mod, only: csr_param
use location_encode_mod
use transfer_map_mod
use opti_de_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_branch_struct), pointer :: lat_branch
type (tao_lattice_struct), pointer :: tao_lat
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_building_wall_section_struct), pointer :: section
type (tao_building_wall_point_struct), pointer :: pt
type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, save, target :: v_array(:)
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: c1
type (tao_plot_region_struct), pointer :: region
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_ele_shape_struct), pointer :: shape
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (ele_struct), target :: ele3, ele0
type (bunch_struct), pointer :: bunch
type (wake_lr_struct), pointer :: lr
type (coord_struct), target :: orb
type (bunch_params_struct) bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (taylor_struct) taylor(6)
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch
type (wall3d_struct), pointer :: wall
type (wall3d_section_struct), pointer :: wall_sec
type (wall3d_vertex_struct), pointer :: v
type (random_state_struct) ran_state
type (normal_form_struct), pointer :: normal_form
type (aperture_scan_struct), pointer :: aperture_scan

type show_lat_column_struct
  character(80) name
  character(16) format
  integer field_width
  character(32) label
  logical remove_line_if_zero
end type

type (show_lat_column_struct) column(50)

real(rp) f_phi, s_pos, l_lat, gam, s_ele, s1, s2, gamma2, val, z, angle, r
real(rp) mat6(6,6), vec0(6), vec_in(6)
real(rp), allocatable, save :: value(:)

character(*) :: what, stuff
character(3) undef_str
character(20) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name, what_to_print
character(24) :: var_name, blank_str = ''
character(24)  :: plane, imt, lmt, amt, iamt, ramt, f3mt, rmt, irmt, iimt
character(40) ele_name, sub_name, ele1_name, ele2_name, switch
character(40) replacement_for_blank
character(60) nam
character(80) :: word1, fmt, fmt2, fmt3
character(100) file_name, name, why_invalid
character(120) header, str
character(200), allocatable :: alloc_lines(:)

character(16) :: show_what, show_names(31) = [ &
   'data            ', 'variable        ', 'global          ', 'alias           ', 'top10           ', &
   'optimizer       ', 'element         ', 'lattice         ', 'constraints     ', 'plot            ', &
   'beam            ', 'tune            ', 'graph           ', 'curve           ', 'particle        ', &
   'hom             ', 'key_bindings    ', 'universe        ', 'orbit           ', 'derivative      ', &
   'branch          ', 'use             ', 'taylor_map      ', 'value           ', 'wave            ', &
   'twiss_and_orbit ', 'building_wall   ', 'wall            ', 'normal_form     ', 'dynamic_aperture', &
   'matrix          ']

character(*), allocatable :: lines(:)
character(*) result_id
character(n_char_show) line, line1, line2, line3
character(n_char_show) stuff2
character(9) angle_str

integer :: data_number, ix_plane, ix_class, n_live, n_order, i1, i2, ix_branch
integer nl, loc, ixl, iu, nc, n_size, ix_u, ios, ie, nb, id, iv, jd, jv, stat, lat_type
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju, ios1, ios2, i_uni, ix_remove
integer num_locations, ix_ele, n_name, n_start, n_ele, n_ref, n_tot, ix_p, ix_word
integer xfer_mat_print, twiss_out, ix_sec

logical bmad_format, good_opt_only, show_lords, print_wall, show_lost
logical err, found, at_ends, first_time, by_s, print_header_lines, all_lat, limited
logical show_sym, show_line, show_shape, print_data, ok, print_tail_lines, print_slaves
logical show_all, name_found, print_taylor, print_em_field, print_all, print_ran_state
logical print_global, print_optimization, print_bmad_com, print_csr_param, print_ptc
logical valid_value, print_floor, show_section, is_complex
logical, allocatable, save :: picked_uni(:)
logical, allocatable, save :: picked_ele(:)
logical, allocatable, save :: good(:)

namelist / custom_show_list / column

!

call re_allocate (lines, 200)

err = .false.

lines = " "
nl = 0

rmt  = '(a, 9es16.8)'
f3mt = '(a, 9f0.3)'
irmt = '(a, i0, a, es16.8)'
imt  = '(a, 9i8)'
iimt = '(a, i0, a, i8)'
lmt  = '(a, 9(l1, 2x))'
amt  = '(9a)'
iamt = '(a, i0, 2x, 9a)'
ramt = '(a, f0.3, 2x, 9a)'

u => tao_pointer_to_universe(-1)
lat => u%model%lat
branch => lat%branch(s%com%default_branch)

if (s%global%phase_units == radians$) f_phi = 1
if (s%global%phase_units == degrees$) f_phi = 180 / pi
if (s%global%phase_units == cycles$)  f_phi = 1 / twopi

! find what to show

result_id = 'ERROR'

if (what == ' ') then
  nl=1; lines(1) = 'SHOW WHAT?' 
  return
endif

call match_word (what, show_names, ix, matched_name = show_what)
if (ix == 0) then
  nl=1; lines(1) = 'SHOW WHAT? WORD NOT RECOGNIZED: ' // what
  return
endif

if (ix < 0) then
  nl=1; lines(1) = 'SHOW WHAT? AMBIGUOUS: ' // what
  return
endif

call string_trim (stuff, stuff2, ix_word)
word1 = stuff2(:ix_word)

select case (show_what)

!----------------------------------------------------------------------
! alias

case ('alias')

  call re_allocate (lines, s%com%n_alias+10, .false.)
  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, s%com%n_alias
    nl=nl+1; lines(nl) = trim(s%com%alias(i)%name) // ' = "' // &
                                    trim(s%com%alias(i)%string) // '"'
  enddo
  
  result_id = show_what

!----------------------------------------------------------------------
! beam

case ('beam')

  ! no element index

  if (word1 == '') then

    ix_branch = s%com%default_branch
    uni_branch => u%uni_branch(ix_branch)
    nl=nl+1; write(lines(nl), '(a, i0, a, i0)') 'Universe: ', u%ix_uni, '  of: ', ubound(s%u, 1)
    nl=nl+1; write(lines(nl), '(a, i3)') 'Branch:   ', ix_branch
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), amt) 's%com%beam_file           = ', s%com%beam_file
    nl=nl+1; write(lines(nl), amt) 'beam0_file                  = ', u%beam%beam0_file
    nl=nl+1; write(lines(nl), amt) 'beam_all_file               = ', u%beam%beam_all_file
    beam => uni_branch%ele(0)%beam
    if (allocated(beam%bunch)) then
      nl=nl+1; write(lines(nl), imt) 'n_particle                  = ', size(beam%bunch(1)%particle)
      nl=nl+1; write(lines(nl), imt) 'n_bunch                     = ', size(beam%bunch)
      nl=nl+1; write(lines(nl), rmt) 'bunch_charge_tot            = ', beam%bunch(1)%charge_tot
      nl=nl+1; write(lines(nl), amt) 'bunch_species               = ', trim(species_name(beam%bunch(1)%particle(1)%species))
    endif
    if (u%beam%beam_all_file == '' .and. u%beam%beam0_file == '') then
      beam_init => u%beam%beam_init
      nl=nl+1; lines(nl) = 'beam_init components:'
      nl=nl+1; write(lines(nl), amt) '  %distribution_type      = ', beam_init%distribution_type
      nl=nl+1; write(lines(nl), rmt) '  %center                 = ', beam_init%center
      nl=nl+1; write(lines(nl), rmt) '  %a_norm_emit            = ', beam_init%a_norm_emit
      nl=nl+1; write(lines(nl), rmt) '  %b_norm_emit            = ', beam_init%b_norm_emit
      nl=nl+1; write(lines(nl), rmt) '  %a_emit                 = ', beam_init%a_emit
      nl=nl+1; write(lines(nl), rmt) '  %b_emit                 = ', beam_init%b_emit
      nl=nl+1; write(lines(nl), rmt) '  %dPz_dz                 = ', beam_init%dPz_dz
      nl=nl+1; write(lines(nl), rmt) '  %dt_bunch               = ', beam_init%dt_bunch
      nl=nl+1; write(lines(nl), rmt) '  %sig_z                  = ', beam_init%sig_z
      nl=nl+1; write(lines(nl), rmt) '  %sig_e                  = ', beam_init%sig_e
      nl=nl+1; write(lines(nl), rmt) '  %center_jitter          = ', beam_init%center_jitter
      nl=nl+1; write(lines(nl), rmt) '  %emit_jitter            = ', beam_init%emit_jitter
      nl=nl+1; write(lines(nl), rmt) '  %sig_z_jitter           = ', beam_init%sig_z_jitter
      nl=nl+1; write(lines(nl), rmt) '  %sig_e_jitter           = ', beam_init%sig_e_jitter
      nl=nl+1; write(lines(nl), rmt) '  %spin%polarization      = ', beam_init%spin%polarization
      nl=nl+1; write(lines(nl), rmt) '  %spin%theta             = ', beam_init%spin%theta
      nl=nl+1; write(lines(nl), rmt) '  %spin%phi               = ', beam_init%spin%phi
      nl=nl+1; write(lines(nl), lmt) '  %renorm_center          = ', beam_init%renorm_center
      nl=nl+1; write(lines(nl), lmt) '  %renorm_sigma           = ', beam_init%renorm_sigma
      nl=nl+1; write(lines(nl), amt) '  %random_engine          = ', beam_init%random_engine
      nl=nl+1; write(lines(nl), amt) '  %random_gauss_converter = ', beam_init%random_gauss_converter
      nl=nl+1; write(lines(nl), f3mt)'  %random_sigma_cutoff    = ', beam_init%random_sigma_cutoff
      fmt = '(a, i1, a, es16.8)'
      do i = 1, 3
        if (beam_init%distribution_type(i) == 'ELLIPSE') then
          nl=nl+1; write(lines(nl), iimt) '  %ellipse(', i, ')%part_per_ellipse  = ', beam_init%ellipse(i)%part_per_ellipse
          nl=nl+1; write(lines(nl), iimt) '  %ellipse(', i, ')%n_ellipse         = ', beam_init%ellipse(i)%n_ellipse
          nl=nl+1; write(lines(nl), irmt) '  %ellipse(', i, ')%sigma_cutoff      = ', beam_init%ellipse(i)%sigma_cutoff
        elseif (beam_init%distribution_type(i) == 'GRID') then
          nl=nl+1; write(lines(nl), iimt) '  %grid(', i, ')%n_x            = ', beam_init%grid(i)%n_x
          nl=nl+1; write(lines(nl), iimt) '  %grid(', i, ')%n_px           = ', beam_init%grid(i)%n_px
          nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%x_min          = ', beam_init%grid(i)%x_min
          nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%x_max          = ', beam_init%grid(i)%x_max
          nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%px_min         = ', beam_init%grid(i)%px_min
          nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%px_max         = ', beam_init%grid(i)%px_max
        endif
      enddo
      if (any(beam_init%distribution_type == 'KV')) then
        nl=nl+1; write(lines(nl), imt) '  %kv%part_per_phi(1:2) = ', beam_init%kv%part_per_phi
        nl=nl+1; write(lines(nl), imt) '  %kv%n_i2              = ', beam_init%kv%n_i2
        nl=nl+1; write(lines(nl), rmt) '  %kv%a                 = ', beam_init%kv%a
      endif
    endif
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'bmad_com components:'
    nl=nl+1; write(lines(nl), lmt) '  %sr_wakes_on               = ', bmad_com%sr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %lr_wakes_on               = ', bmad_com%lr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %space_charge_on           = ', bmad_com%space_charge_on
    nl=nl+1; write(lines(nl), lmt) '  %coherent_synch_rad_on     = ', bmad_com%coherent_synch_rad_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_tracking_on          = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_damping_on      = ', bmad_com%radiation_damping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_fluctuations_on = ', bmad_com%radiation_fluctuations_on
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'csr_param components:'
    nl=nl+1; write(lines(nl), rmt) '  %ds_track_step        = ', csr_param%ds_track_step
    nl=nl+1; write(lines(nl), rmt) '  %beam_chamber_height  = ', csr_param%beam_chamber_height
    nl=nl+1; write(lines(nl), rmt) '  %sigma_cutoff         = ', csr_param%sigma_cutoff
    nl=nl+1; write(lines(nl), imt) '  %n_bin                = ', csr_param%n_bin
    nl=nl+1; write(lines(nl), imt) '  %particle_bin_span    = ', csr_param%particle_bin_span
    nl=nl+1; write(lines(nl), imt) '  %n_shield_images      = ', csr_param%n_shield_images
    nl=nl+1; write(lines(nl), imt) '  %ix1_ele_csr          = ', csr_param%ix1_ele_csr
    nl=nl+1; write(lines(nl), imt) '  %ix2_ele_csr          = ', csr_param%ix2_ele_csr
    nl=nl+1; write(lines(nl), lmt) '  %lcsr_component_on    = ', csr_param%lcsr_component_on
    nl=nl+1; write(lines(nl), lmt) '  %lsc_component_on     = ', csr_param%lsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %tsc_component_on     = ', csr_param%tsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %small_angle_approx   = ', csr_param%small_angle_approx
    nl=nl+1; lines(nl) = ''
    if (lat%param%particle == photon$) then
      nl=nl+1; write(lines(nl), rmt) 'model%lat%a%emit               = ', lat%a%emit
      nl=nl+1; write(lines(nl), rmt) 'model%lat%b%emit               = ', lat%b%emit
    else
      call convert_total_energy_to (lat%ele(0)%value(e_tot$), lat%param%particle, gamma = gam)
      nl=nl+1; write(lines(nl), rmt) 'model%lat%a%emit               = ', lat%a%emit
      nl=nl+1; write(lines(nl), rmt) '          a%emit (normalized)  = ', lat%a%emit * gam
      nl=nl+1; write(lines(nl), rmt) 'model%lat%b%emit               = ', lat%b%emit
      nl=nl+1; write(lines(nl), rmt) '          b%emit (normalized)  = ', lat%b%emit * gam
    endif
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), amt) 'global%track_type          = ', s%global%track_type
    nl=nl+1; write(lines(nl), lmt) 'global%beam_timer_on       = ', s%global%beam_timer_on
    nl=nl+1; write(lines(nl), amt) 'track_start                = ', trim(uni_branch%track_start)
    nl=nl+1; write(lines(nl), amt) 'track_end                  = ', trim(uni_branch%track_end)
    nl=nl+1; write(lines(nl), imt) 'ix_track_start             = ', uni_branch%ix_track_start 
    nl=nl+1; write(lines(nl), imt) 'ix_track_end               = ', uni_branch%ix_track_end
    nl=nl+1; write(lines(nl), amt) 'u%beam%saved_at:           = ', trim(u%beam%saved_at)

  ! have element index

  else
    call tao_pick_universe (word1, word1, picked_uni, err, ix_u)
    if (err) return
    u => s%u(ix_u)
    call tao_locate_elements (word1, ix_u, eles, err, ix_dflt_branch = s%com%default_branch)
    if (err .or. size(eles) == 0) return
    ix_ele = eles(1)%ele%ix_ele
    ix_branch = eles(1)%ele%ix_branch
    n = s%global%bunch_to_plot

    bunch_p => u%model%lat_branch(ix_branch)%bunch_params(ix_ele)
    n_live = bunch_p%n_particle_live
    n_tot = bunch_p%n_particle_tot

    if (n_tot == 0) then
      nl=nl+1; lines(nl) = 'Beam has no particles!'
      result_id = 'beam:no_particles'
      return
    endif

    nl=nl+1; lines(nl) = 'Cashed bunch parameters:'
    nl=nl+1; write(lines(nl), imt)  '  Parameters for bunch:       ', n
    nl=nl+1; write(lines(nl), imt)  '  Particles surviving:        ', n_live
    nl=nl+1; write(lines(nl), imt)  '  Particles lost:             ', n_tot - n_live
    nl=nl+1; write(lines(nl), f3mt) '  Particles lost (%):         ', 100 * real(n_tot - n_live) / n_tot
    if (u%model%lat%branch(ix_branch)%param%particle == photon$) then
      nl=nl+1; write(lines(nl), rmt)  '  Intensity:                  ', &
                        bunch_p%centroid%field(1)**2 + bunch_p%centroid%field(2)**2
    endif
    nl=nl+1; write(lines(nl), rmt) '  Centroid:', bunch_p%centroid%vec
    nl=nl+1; write(lines(nl), rmt) '  RMS:     ', &
                      sqrt(bunch_p%sigma(1,1)), sqrt(bunch_p%sigma(2,2)), sqrt(bunch_p%sigma(3,3)), &
                      sqrt(bunch_p%sigma(4,4)), sqrt(bunch_p%sigma(5,5)), sqrt(bunch_p%sigma(6,6))
    if (u%model%lat%branch(eles(1)%ele%ix_branch)%param%particle /= photon$) then
      nl=nl+1; write(lines(nl), rmt) '             norm_emitt           beta             alpha'
      nl=nl+1; write(lines(nl), rmt) '  a:       ', bunch_p%a%norm_emit, bunch_p%a%beta, bunch_p%a%alpha
      nl=nl+1; write(lines(nl), rmt) '  b:       ', bunch_p%b%norm_emit, bunch_p%b%beta, bunch_p%b%alpha
      nl=nl+1; write(lines(nl), rmt) '  x:       ', bunch_p%x%norm_emit, bunch_p%x%beta, bunch_p%x%alpha
      nl=nl+1; write(lines(nl), rmt) '  y:       ', bunch_p%y%norm_emit, bunch_p%y%beta, bunch_p%y%alpha
      nl=nl+1; write(lines(nl), rmt) '  z:       ', bunch_p%z%norm_emit, bunch_p%z%beta
    endif

    beam => u%uni_branch(eles(1)%ele%ix_branch)%ele(ix_ele)%beam
    if (allocated(beam%bunch)) then
      bunch => beam%bunch(n)
      nl=nl+1; lines(nl) = 'Note: The beam distribution is saved at this element.'
    else
      nl=nl+1; lines(nl) = 'Note: The beam distribution is not saved at this element.'
    endif
  
  endif

  result_id = show_what

!----------------------------------------------------------------------
! constraints

case ('branch')

  sub_name = ''

  do 

    call tao_next_switch (stuff2, ['-branch     ', '-universe   '], switch, err, ix_s2)

    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-branch')
      sub_name = stuff2(1:ix_s2)
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-universe')
      read (stuff2(1:ix_s2), *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (ix_s2 == 0 .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" argument'
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    end select
  enddo

  if (sub_name == '') then
    sub_name = stuff2
  elseif (stuff2 /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // stuff2)
    return
  endif

  lat => u%model%lat

  if (size(s%u) > 1) then
    nl=nl+1; write(lines(nl), '()') 'For the lattice of universe: ', ix_u
  endif

  if (sub_name /= '') then
    branch => pointer_to_branch(sub_name, u%model%lat)
    if (.not. associated(branch)) then
      nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
      return
    endif

    nl=nl+1; lines(nl) =             '%name                      = ' // trim(branch%name)
    nl=nl+1; write(lines(nl), imt)   '%ix_branch                 =', branch%ix_branch 
    nl=nl+1; write(lines(nl), imt)   '%ix_from_branch            =', branch%ix_from_branch 
    nl=nl+1; write(lines(nl), imt)   '%ix_from ele               =', branch%ix_from_ele
    nl=nl+1; write(lines(nl), imt)   '%n_ele_track               =', branch%n_ele_track
    nl=nl+1; write(lines(nl), imt)   '%n_ele_max                 =', branch%n_ele_max
    nl=nl+1; write(lines(nl), imt)   '%param%particle            = ' // trim(species_name(branch%param%particle))
    nl=nl+1; write(lines(nl), '(a, f6.1)')  '%param%rel_tracking_charge =', branch%param%rel_tracking_charge
    nl=nl+1; write(lines(nl), amt)   '%param%geometry            = ', geometry_name(branch%param%geometry)
    if (branch%param%particle == photon$) then
      nl=nl+1; write(lines(nl), amt) 'lat%photon_type            = ', photon_type_name(lat%photon_type)
    endif

  else
    nl=nl+1; lines(nl) = 'Branch  Branch                N_ele  N_ele  Ix_From  From                Ix_From  From '
    nl=nl+1; lines(nl) = ' Index  Name                  Track    Max   Branch  Branch                  Ele  Ele'


    branch => lat%branch(0)
    nl=nl+1; write(lines(nl), '(i6, 2x, a21, i6, i7)') &
                  0, branch%name, branch%n_ele_track, branch%n_ele_max

    fmt = '(i6, 2x, a21, i6, i7, i9, 3x, a20, i6, 3x, a)'
    do i = 1, ubound(lat%branch, 1)
      branch => lat%branch(i)
      if (branch%ix_from_ele < 0) then
        nl=nl+1; write(lines(nl), fmt) i, branch%name, branch%n_ele_track, branch%n_ele_max, &
                  branch%ix_from_branch
      else
        ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
        nl=nl+1; write(lines(nl), fmt) i, branch%name, branch%n_ele_track, branch%n_ele_max, &
                  branch%ix_from_branch, lat%branch(ele%ix_branch)%name, ele%ix_ele, ele%name
      endif
    enddo
  endif

  result_id = show_what

!----------------------------------------------------------------------
! building_wall

case ('building_wall')

  if (.not. allocated(s%building_wall%section)) then
    nl=nl+1; lines(nl) = 'No building wall defined.'
    result_id = 'building_wall:none'
    return
  endif

  do i = 1, size(s%building_wall%section)
    section => s%building_wall%section(i)
    n = nl + size(section%point)
    if (n + 10 > size(lines)) call re_allocate (lines, n, .false.)
    nl=nl+1; write(lines(nl), '(a, i0, 2a)') 'Section(', i, ')   constraint: ', section%constraint
    do j = 1, size(section%point)
      pt => section%point(j)
      if (pt%radius == 0) then
        nl=nl+1; write(lines(nl), '(a, i0, a, 3f10.3)') &
              '  point(', j, ') z, x, r: ', pt%z, pt%x, pt%radius
      else
        nl=nl+1; write(lines(nl), '(a, i0, a, 5f10.3)') &
              '  point(', j, ') z, x, r: ', pt%z, pt%x, pt%radius, pt%z_center, pt%x_center
      endif
    enddo
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, '*')
  call tao_show_constraints (0, 'TOP10')
  result_id = show_what

!----------------------------------------------------------------------
! curve

case ('curve')

  ! Look for switches

  show_sym = .false.
  show_line = .false.

  do
    call tao_next_switch (stuff2, ['-symbol', '-line  ' ], switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-symbol') show_sym = .true.
    if (switch == '-line')   show_line = .true.
  enddo

  ! Find particular plot

  call tao_find_plots (err, stuff2, 'BOTH', curve = curve, always_allocate = .true.)
  if (err) return

  ! print info on particular plot, graph, or curve

  if (allocated(curve)) then
    c1 => curve(1)%c
    nl=nl+1; lines(nl) = 'Region.Graph.Curve: ' // trim(tao_curve_name(c1, .true.))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c, .true.))
    enddo
    nl=nl+1; lines(nl) = 'Plot.Graph.Curve:   ' // trim(tao_curve_name(c1))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c))
    enddo
    nl=nl+1; write(lines(nl), amt)  'data_source          = ', c1%data_source
    nl=nl+1; write(lines(nl), amt)  'data_index           = ', c1%data_index
    nl=nl+1; write(lines(nl), amt)  'data_type_x          = ', c1%data_type_x
    nl=nl+1; write(lines(nl), amt)  'data_type            = ', c1%data_type
    nl=nl+1; write(lines(nl), amt)  'legend_text          = ', c1%legend_text
    nl=nl+1; write(lines(nl), amt)  'ele_ref_name         = ', c1%ele_ref_name
    nl=nl+1; write(lines(nl), imt)  'ix_branch            = ', c1%ix_branch
    nl=nl+1; write(lines(nl), imt)  'ix_ele_ref           = ', c1%ix_ele_ref
    nl=nl+1; write(lines(nl), imt)  'ix_ele_ref_track     = ', c1%ix_ele_ref_track
    nl=nl+1; write(lines(nl), imt)  'ix_bunch             = ', c1%ix_bunch
    nl=nl+1; write(lines(nl), imt)  'ix_universe          = ', c1%ix_universe
    nl=nl+1; write(lines(nl), imt)  'symbol_every         = ', c1%symbol_every
    nl=nl+1; write(lines(nl), rmt)  'y_axis_scale_factor  = ', c1%y_axis_scale_factor
    nl=nl+1; write(lines(nl), lmt)  'use_y2               = ', c1%use_y2
    nl=nl+1; write(lines(nl), lmt)  'draw_line            = ', c1%draw_line
    nl=nl+1; write(lines(nl), lmt)  'draw_symbols         = ', c1%draw_symbols
    nl=nl+1; write(lines(nl), lmt)  'draw_symbol_index    = ', c1%draw_symbol_index
    nl=nl+1; write(lines(nl), lmt)  'smooth_line_calc     = ', c1%smooth_line_calc
    nl=nl+1; write(lines(nl), iamt) 'line%width           = ', c1%line%width
    nl=nl+1; write(lines(nl), iamt) 'line%color           = ', c1%line%color, qp_color_name(c1%line%color)
    nl=nl+1; write(lines(nl), iamt) 'line%pattern         = ', c1%line%pattern, qp_line_pattern_name(c1%line%pattern)
    nl=nl+1; write(lines(nl), iamt) 'symbol%type          = ', c1%symbol%type, qp_symbol_type_name(c1%symbol%type)
    nl=nl+1; write(lines(nl), f3mt) 'symbol%height        = ', c1%symbol%height
    nl=nl+1; write(lines(nl), iamt) 'symbol%fill_pattern  = ', c1%symbol%fill_pattern, qp_fill_name(c1%symbol%fill_pattern)
    nl=nl+1; write(lines(nl), iamt) 'symbol%line_width    = ', c1%symbol%line_width
    
    ! Histogram specific components
    if (c1%g%type == 'histogram') then
      nl=nl+1; write(lines(nl), lmt)  'hist%density_normalized = ', c1%hist%density_normalized 
      nl=nl+1; write(lines(nl), lmt)  'hist%weight_by_charge   = ', c1%hist%weight_by_charge
      nl=nl+1; write(lines(nl), rmt)  'hist%minimum            = ', c1%hist%minimum
      nl=nl+1; write(lines(nl), rmt)  'hist%maximum            = ', c1%hist%maximum
      nl=nl+1; write(lines(nl), rmt)  'hist%width              = ', c1%hist%width
      nl=nl+1; write(lines(nl), rmt)  'hist%center             = ', c1%hist%center
      nl=nl+1; write(lines(nl), imt)  'hist%number             = ', c1%hist%number
    endif
    
    if (show_sym) then
      n = nl + size(c1%x_symb) + 10
      if (n > size(lines)) call re_allocate(lines, n, .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Symbol points:'
      err = .false.

      nl=nl+1; lines(nl) = '      i  index        x-axis'
      do j = 1, size(curve)
        str = curve(j)%c%name
        lines(nl) = lines(nl)(1:28+(j-2)*14) // adjustr(str(1:14))
      enddo

      do j = 2, size(curve)
        if (size(curve(j)%c%y_symb) /= size(c1%y_symb)) then
          nl=nl+1; lines(nl) = 'NUMBER OF SYMBOL POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo

      if (.not. err) then
        do i = 1, size(c1%x_symb)
          nl=nl+1; write(lines(nl), '(2i7, 10es14.6)') i, c1%ix_symb(i), &
                      c1%x_symb(i), [ (curve(j)%c%y_symb(i), j = 1, size(curve)) ]
        enddo
      endif
    endif

    if (show_line) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Smooth line points:'
      nl=nl+1; lines(nl) = ' Index        x-axis'

      do j = 1, size(curve)
        str = curve(j)%c%name
        lines(nl) = lines(nl)(1:20+(j-1)*14) // adjustr(str(1:14))
      enddo

      do j = 2, size(curve)
        if (size(curve(j)%c%y_line) /= size(c1%y_line)) then
          nl=nl+1; lines(nl) = 'NUMBER OF LINE POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo

      if (.not. err) then
        call re_allocate (lines, nl+size(c1%x_line)+100, .false.)
        do i = 1, size(c1%x_line)
          nl=nl+1; write(lines(nl), '(i6, 10es14.6)') i, c1%x_line(i), &
                                        [ (curve(j)%c%y_line(i), j = 1, size(curve)) ]
        enddo
      endif
    endif

  else
    nl=1; lines(1) = 'THIS IS NOT A CURVE NAME'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! data

case ('data')


  ! If just "show data" then show all names

  call tao_pick_universe (word1, line1, picked_uni, err)
  if (err) return

  if (line1 == ' ') then  ! just specified a universe

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, t40, a)') 'Name', 'Using for Optimization'

    do iu = lbound(s%u, 1), ubound(s%u, 1)

      if (.not. picked_uni(iu)) cycle

      u => s%u(iu)

      do i = 1, u%n_d2_data_used
        d2_ptr => u%d2_data(i)
        if (d2_ptr%name == ' ') cycle
        call tao_data_show_use (d2_ptr, lines, nl)
      enddo
    enddo

    result_id = 'data:'
    return
  endif

  ! get pointers to the data

  call tao_find_data (err, word1, d2_ptr, d1_array, d_array)
  if (err) return

  n_size = 0
  if (allocated(d_array)) n_size = size(d_array)

  ! If d_ptr points to something then show the datum info.

  if (n_size == 1) then
    d_ptr => d_array(1)%d
    nl=nl+1; lines(nl) = ''
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d_ptr%d1%d2%ix_uni
    endif
    nl=nl+1; write(lines(nl), amt)  '%ele_name          = ', d_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%ele_start_name    = ', d_ptr%ele_start_name
    nl=nl+1; write(lines(nl), amt)  '%ele_ref_name      = ', d_ptr%ele_ref_name
    nl=nl+1; write(lines(nl), amt)  '%data_type         = ', d_ptr%data_type
    nl=nl+1; write(lines(nl), amt)  '%data_source       = ', d_ptr%data_source
    nl=nl+1; write(lines(nl), imt)  '%ix_branch         = ', d_ptr%ix_branch
    nl=nl+1; write(lines(nl), imt)  '%ix_ele            = ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_start      = ', d_ptr%ix_ele_start
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_ref        = ', d_ptr%ix_ele_ref
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_merit      = ', d_ptr%ix_ele_merit
    nl=nl+1; write(lines(nl), imt)  '%ix_dmodel         = ', d_ptr%ix_dModel
    nl=nl+1; write(lines(nl), imt)  '%ix_d1             = ', d_ptr%ix_d1
    nl=nl+1; write(lines(nl), imt)  '%ix_data           = ', d_ptr%ix_data
    nl=nl+1; write(lines(nl), imt)  '%ix_bunch          = ', d_ptr%ix_bunch
    nl=nl+1; write(lines(nl), rmt)  '%model             = ', d_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%design            = ', d_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%meas              = ', d_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%ref               = ', d_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%base              = ', d_ptr%base_value
    nl=nl+1; write(lines(nl), rmt)  '%old               = ', d_ptr%old_value   
    nl=nl+1; write(lines(nl), rmt)  '%invalid           = ', d_ptr%invalid_value
    nl=nl+1; write(lines(nl), rmt)  '%s                 = ', d_ptr%s
    nl=nl+1; write(lines(nl), amt)  '%merit_type        = ', d_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit             = ', d_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit       = ', d_ptr%delta_merit
    nl=nl+1; write(lines(nl), rmt)  '%weight            = ', d_ptr%weight
    nl=nl+1; write(lines(nl), lmt)  '%exists            = ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_model        = ', d_ptr%good_model
    nl=nl+1; write(lines(nl), lmt)  '%good_design       = ', d_ptr%good_design
    nl=nl+1; write(lines(nl), lmt)  '%good_base         = ', d_ptr%good_base 
    nl=nl+1; write(lines(nl), lmt)  '%good_meas         = ', d_ptr%good_meas
    nl=nl+1; write(lines(nl), lmt)  '%good_ref          = ', d_ptr%good_ref
    nl=nl+1; write(lines(nl), lmt)  '%good_user         = ', d_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%good_opt          = ', d_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%good_plot         = ', d_ptr%good_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_plot        = ', d_ptr%useit_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_opt         = ', d_ptr%useit_opt

    if (d_ptr%exists .and. .not. d_ptr%good_model) then
      u => s%u(d_ptr%d1%d2%ix_uni)
      call tao_evaluate_a_datum (d_ptr, u, u%model, val, valid_value, why_invalid)
      nl=nl+1; lines(nl) = 'Model value is invalid since: ' // why_invalid
    endif

  ! Else show the d1_data info.

  elseif (size(d1_array) == 1) then

    d1_ptr => d1_array(1)%d1
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d1_ptr%d2%ix_uni
    endif
    
    nl=nl+1; write(lines(nl), '(2a)') 'Data name: ', trim(d2_ptr%name) // '.' // d1_ptr%name

    ! find string widths
    ! Expressions generally have very long strings so we let this spill over to
    ! the where0 and where fields

    n_name  = 9      ! Set mimimum field widths
    n_start = 10
    n_ref   = 8
    n_ele   = 8

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_constraint_type_name(d_ptr)
      if (d_ptr%data_type(1:11) /= 'expression:') then
        n_name  = max(n_name,  len_trim(name))
        n_start = max(n_start, len_trim(d_ptr%ele_start_name))
        n_ref   = max(n_ref,   len_trim(d_ptr%ele_ref_name))
        n_ele   = max(n_ele,   len_trim(d_ptr%ele_name))
      endif
    enddo

    ! Write header
    ! Element names are left justified and real quantities are right justified

    line1 = ''; line2 = ''
    n=9+n_name;    line2(n:) = 'Ref_Ele'  ! n = i4 + 2x + n_name + 2x + 1
    n=n+n_ref+2;   line2(n:) = 'Start_Ele'
    n=n+n_start+2; line2(n:) = 'Ele'
    n=n+n_ele+10;  line2(n:) = 'Meas         Model        Design | Opt  Plot'
                   line1(n:) = '                                 |   Useit'

    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2

    ! if a range is specified, show the data range   

    call re_allocate (lines, nl+100+size(d_array), .false.)

    fmt  = '(i4, 4(2x, a), 3es14.4, 2l6)'
    fmt2 = '(4x, 4(2x, a), 3es14.4, 2l6)'

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_constraint_type_name(d_ptr)
      if (d_ptr%data_type(1:11) == 'expression:') then
        nl=nl+1; write(lines(nl), fmt) d_ptr%ix_d1, trim(name)
        nl=nl+1; write(lines(nl), fmt2) blank_str(1:n_name), &
                     d_ptr%ele_ref_name(1:n_ref), d_ptr%ele_start_name(1:n_start), &
                     d_ptr%ele_name(1:n_ele), d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
      else
        nl=nl+1; write(lines(nl), fmt) d_ptr%ix_d1, name(1:n_name), &
                     d_ptr%ele_ref_name(1:n_ref), d_ptr%ele_start_name(1:n_start), &
                     d_ptr%ele_name(1:n_ele), d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
      endif
    enddo

    nl=nl+1; lines(nl) = line2
    nl=nl+1; lines(nl) = line1

  ! else we must have a valid d2_ptr.

  elseif (associated(d2_ptr)) then

    call re_allocate (lines, nl+100+size(d2_ptr%d1), .false.)

    nl=nl+1; write(lines(nl), '(t40, a)')     'Using' 

    do i = 1, size(d2_ptr%d1)
      if (size(lines) < nl + 50) call re_allocate (lines, nl+100, .false.)
      call location_encode(line, d2_ptr%d1(i)%d%useit_opt, &
                      d2_ptr%d1(i)%d%exists, lbound(d2_ptr%d1(i)%d, 1))
      nl=nl+1; write(lines(nl), '(2x, 2a, i0, a, i0, a, t40, a)') &
                  trim(tao_d2_d1_name(d2_ptr%d1(i))), '[', lbound(d2_ptr%d1(i)%d, 1), &
                  ':', ubound(d2_ptr%d1(i)%d, 1), ']', trim(line)
    enddo

    if (any(d2_ptr%descrip /= ' ')) then
      call re_allocate (lines, nl+100+size(d2_ptr%descrip), .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Descrip:'
      do i = 1, size(d2_ptr%descrip)
        if (d2_ptr%descrip(i) /= ' ') then
          nl=nl+1; write(lines(nl), '(i4, 2a)') i, ': ', d2_ptr%descrip(i)
        endif
      enddo
    endif

  ! error

  else
    nl=1; lines(1) = 'TRY BEING MORE SPECIFIC.'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! derivative

case ('derivative')

  call tao_find_data (err, word1, d_array = d_array)
  if (err) return

  call string_trim(stuff2(ix_word+1:), stuff2, ix_word)
  call tao_find_var(err, stuff2(:ix_word), v_array = v_array) 
  if (err) return

  i1 = 0
  do id = 1, size(d_array)
    if (d_array(id)%d%ix_dmodel == 0) cycle
    i1 = i1 + 1
  enddo

  i2 = 0
  do iv = 1, size(v_array)
    if (v_array(iv)%v%ix_dvar == 0) cycle
    i2 = i2 + 1
  enddo

  call re_allocate (lines, nl+i1*i2+10, .false.)

  found = .false.
  do id = 1, size(d_array)
    do iv = 1, size(v_array)
      d_ptr => d_array(id)%d
      v_ptr => v_array(iv)%v
      u => s%u(d_ptr%d1%d2%ix_uni)
      jd = d_ptr%ix_dmodel
      jv = v_ptr%ix_dvar
      if (jd /= 0 .and. jv /= 0) then
        nl=nl+1; write(lines(nl), '(2a20, es14.5, 2i5)') tao_datum_name(d_ptr), &
                                  tao_var1_name(v_ptr), u%dModel_dVar(jd, jv), jd, jv
        found = .true.
      endif
    enddo
  enddo

  if (.not. found) then
    nl=nl+1; write(lines(nl), '(a)') 'No Derivatives'
  endif
  
  result_id = show_what

!----------------------------------------------------------------------
! dynamic_aperture

case ('dynamic_aperture')
  if (.not. allocated(u%dynamic_aperture%scan) ) then
    nl=nl+1; lines(nl) ='No dynamic aperture specified for this universe'
    return
  endif

  !call tao_next_switch (stuff2, ['-order'], switch, err, ix)
  !if (err) return
  do i = 1, size(u%dynamic_aperture%scan)
    aperture_scan => u%dynamic_aperture%scan(i) 
    nl=nl+1; write(lines(nl), '(a12, es15.7)')  'pz        : ', u%dynamic_aperture%pz(i)
    nl=nl+1; write(lines(nl), '(a12, i10)')     'n_angle   : ', aperture_scan%n_angle
    nl=nl+1; write(lines(nl), '(a12, f10.6)')   'min_angle : ', aperture_scan%min_angle
    nl=nl+1; write(lines(nl), '(a12, f10.6)')   'max_angle : ', aperture_scan%max_angle
    nl=nl+1; write(lines(nl), '(a20, i10)')      '%praram%n_turn :  ', aperture_scan%param%n_turn
    nl=nl+1; write(lines(nl), '(a20, f10.6)')   '%param%accuracy : ', aperture_scan%param%accuracy
    if (.not. allocated(aperture_scan%aperture)) then
      nl=nl+1; write(lines(nl), '(a20)') 'aperture not calculated for this universe'
    else
      nl=nl+1; write(lines(nl), '(2a15)') 'aperture.x', 'aperture.y' 
      do j = 1, size(aperture_scan%aperture)
        nl=nl+1; write(lines(nl), '(2es15.7)')   aperture_scan%aperture(j)%x, aperture_scan%aperture(j)%y
      enddo
    endif
  enddo



  result_id = show_what

!----------------------------------------------------------------------
! ele

case ('element')

  print_floor = .false.
  print_taylor = .false.
  print_em_field = .false.
  print_all = .false.
  print_data = .false.
  print_wall = .false.
  xfer_mat_print = 0
  print_slaves = .true.
  lat_type = model$
  print_ptc = .false.

  do
    call tao_next_switch (stuff2, ['-taylor        ', '-em_field      ', &
                '-all_attributes', '-data          ', '-design        ', &
                '-no_slaves     ', '-wall          ', '-base          ', &
                '-field         ', '-floor_coords  ', '-xfer_mat      ', &
                '-ptc           ', '-everything    '], switch, err, ix)
    if (err) return
    select case (switch)
    case ('');                exit
    case ('-xfer_mat');       xfer_mat_print = 6
    case ('-floor_coords');   print_floor = .true.
    case ('-taylor');         print_taylor = .true.
    case ('-design');         lat_type = design$
    case ('-base');           lat_type = base$
    case ('-em_field');       print_em_field = .true.  ! Old style. Use "-field".
    case ('-field');          print_em_field = .true.
    case ('-all_attributes'); print_all = .true.
    case ('-data');           print_data = .true.
    case ('-no_slaves');      print_slaves = .false.
    case ('-wall');           print_wall = .true.
    case ('-ptc');            print_ptc = .true.
    case ('-everything')
      print_all = .true.
      xfer_mat_print = 6
      print_taylor = .true.
      print_floor = .true.
      print_em_field = .true.
      print_wall = .true.
    case default
      call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND: ' // switch)
      return
    end select
  enddo

  call str_upcase (ele_name, stuff2)
  call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
  if (err) return
  u => s%u(ix_u)

  ! Wildcard: show all elements.

  if (index(ele_name, '*') /= 0 .or. index(ele_name, '%') /= 0 .or. &
                (ele_name(1:2) /= 'S:' .and. index(ele_name, ':') /= 0) .or. &
                count(picked_uni) > 1) then

    n_tot = 0
    do i_uni = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. picked_uni(i_uni)) cycle
      call tao_locate_elements (ele_name, i_uni, eles, err, &
                                      ignore_blank = .true., ix_dflt_branch = s%com%default_branch)
      if (err) return
      lat => s%u(i_uni)%model%lat
      do i = 1, size(eles)
        ele => eles(i)%ele
        if (.not. print_slaves) then
          if (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) cycle
        endif
        if (size(lines) < nl+100) call re_allocate (lines, nl+200, .false.)
        n_tot = n_tot + 1
        if (size(ele%branch%lat%branch) == 1) then
          write (str, '(i10)') ele%ix_ele
        else
          write (str, '(i0, a, i0)') ele%ix_branch, '>>', ele%ix_ele
          str(1:10) = adjustr(str(1:10))
        endif
        if (count(picked_uni) > 1) then
          nl=nl+1; write(lines(nl), '(a10, 2x, i0, 2a)') str, i_uni, '@', ele%name
        else
          nl=nl+1; write(lines(nl), '(a10, 2x, a)') str, ele%name
        endif
      enddo
    enddo

    deallocate(eles)
    nl=nl+1; write(lines(nl), '(a, i0)') 'Number of Matches: ', n_tot

    if (nl == 0) then
      lines(1) = '*** No Matches to Name Found ***'
      return
    endif

    result_id = 'element:*'
    return

  endif

  ! No wildcard case...
  ! Normal: Show the element info

  call tao_locate_elements (ele_name, ix_u, eles, err, lat_type, ix_dflt_branch = s%com%default_branch)
  if (err) return
  ele => eles(1)%ele

  tao_lat => tao_pointer_to_tao_lat (u, lat_type)
  branch => tao_lat%lat%branch(ele%ix_branch)

  ! Show data associated with this element

  if (print_data) then
    call show_ele_data (u, ele, lines, nl)
    result_id = 'element:data'
    return
  endif

  if (print_ptc) then
    if (.not. associated (ele%ptc_fibre)) then
      nl=nl+1; lines(nl) = 'Creating associated Fibre...'
      call ele_to_fibre (ele, ele%ptc_fibre, ele%branch%param, .true.)
    endif
    call type_ptc_fibre (ele%ptc_fibre, .true., alloc_lines, n)
    if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
    lines(nl+1:nl+n) = alloc_lines(1:n)
    nl = nl + n
    result_id = 'element:ptc'
    return
  endif

  if (s%com%common_lattice) then
    u%calc%lattice = .true.
    call tao_lattice_calc (ok)
  endif

  twiss_out = s%global%phase_units
  if (lat%branch(ele%ix_branch)%param%particle == photon$) twiss_out = 0
  call type2_ele (ele, alloc_lines, n, print_all, xfer_mat_print, print_taylor, &
            twiss_out, .true., .true., print_floor, print_em_field, print_wall)

  if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
  lines(nl+1:nl+n) = alloc_lines(1:n)
  nl = nl + n

  stat = ele%lord_status
  if (stat /= multipass_lord$ .and. stat /= group_lord$ .and. stat /= overlay_lord$ .and. stat /= girder_lord$) then
    orb = tao_lat%lat_branch(ele%ix_branch)%orbit(ele%ix_ele)
    nl=nl+1; lines(nl) = ' '
    nl=nl+1; write(lines(nl), '(4a)') 'Orbit:  ', trim(species_name(orb%species)), '   State: ', trim(coord_state_name(orb%state))
    if (lat%branch(ele%ix_branch)%param%particle == photon$) then
      fmt = '(2x, a, 2f15.8, f15.6, f11.6, 7x, a, f11.3)'
      fmt2 = '(2x, a, 2f15.8, a, es16.8)'
      nl=nl+1; lines(nl) = '         Position[mm]            V/C      Intensity      Phase  '
      nl=nl+1; write(lines(nl), fmt)  'X:  ', orb%vec(1:2), orb%field(1)**2, orb%phase(1), 'E: ', orb%p0c
      nl=nl+1; write(lines(nl), fmt)  'Y:  ', orb%vec(3:4), orb%field(2)**2, orb%phase(2), 'dE:', orb%p0c - ele%value(p0c$)
      nl=nl+1; write(lines(nl), fmt2) 'Z:  ', orb%vec(5:6)
    else
      z = (ele%ref_time - orb%t) * orb%beta * c_light
      fmt = '(2x, a, 2f15.8, a, es16.8, 2x, a, f9.6)'
      nl=nl+1; lines(nl) = '         Position[mm] Momentum[mrad]  |                            Time'
      nl=nl+1; write(lines(nl), fmt) 'X:  ', 1000*orb%vec(1:2),   '  | Absolute [sec]:   ', orb%t
      nl=nl+1; write(lines(nl), fmt) 'Y:  ', 1000*orb%vec(3:4),   '  | Abs-Ref [sec]:    ', orb%t - ele%ref_time
      nl=nl+1; write(lines(nl), fmt) 'Z:  ', 1000*orb%vec(5:6),   '  | (Ref-Abs)*Vel [m]:', z, 'beta:', orb%beta
    endif
  endif

  found = .false.
  do i = 2, size(eles)
    if (size(lines) < nl+2) call re_allocate (lines, nl+10, .false.)
    if (found) then
      nl=nl+1; lines(nl) = ''
      found = .true.
    endif
    nl=nl+1
    if (eles(i)%ele%ix_branch == 0) then
      write(lines(nl), '(a, i0)') &
              'Note: Found another element with same name at: ', eles(i)%ele%ix_ele
    else
      write(lines(nl), '(a, i0, a,i0)') &
              'Note: Found another element with same name at: ', &
              eles(i)%ele%ix_branch, '.', eles(i)%ele%ix_ele
    endif
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! global

case ('global')

  print_global = .true.
  print_optimization = .false.
  print_bmad_com = .false.
  print_csr_param = .false.
  print_ran_state = .false.

  do
    call tao_next_switch (stuff2, ['-optimization', '-bmad_com    ', &
                                   '-csr_param   ', '-ran_state   '], switch, err, ix)
    if (err) return
    if (switch == '') exit
    print_global = .false.
    select case (switch)
    case ('-optimization') 
      print_optimization = .true.
    case ('-bmad_com') 
      print_bmad_com = .true.
    case ('-csr_param') 
      print_csr_param = .true.
    case ('-ran_state')
      print_ran_state = .true.
      print_global = .false.
    end select
  enddo

  if (print_global) then
    nl=nl+1; lines(nl) = 'Global parameters:'
    nl=nl+1; write(lines(nl), imt) '  %bunch_to_plot                 = ', s%global%bunch_to_plot
    nl=nl+1; write(lines(nl), lmt) '  %label_lattice_elements        = ', s%global%label_lattice_elements
    nl=nl+1; write(lines(nl), lmt) '  %label_keys                    = ', s%global%label_keys
    nl=nl+1; write(lines(nl), amt) '  %phase_units                   = ', &
                                                    frequency_units_name(s%global%phase_units)
    nl=nl+1; write(lines(nl), lmt) '  %plot_on                       = ', s%global%plot_on
    nl=nl+1; write(lines(nl), lmt) '  %lattice_calc_on               = ', s%global%lattice_calc_on
    nl=nl+1; write(lines(nl), lmt) '  %command_file_print_on         = ', s%global%command_file_print_on
    nl=nl+1; write(lines(nl), lmt) '  %beam_timer_on                 = ', s%global%beam_timer_on
    nl=nl+1; write(lines(nl), lmt) '  %rf_on                         = ', s%global%rf_on
    nl=nl+1; write(lines(nl), lmt) '  %wait_for_CR_in_single_mode    = ', s%global%wait_for_CR_in_single_mode
    nl=nl+1; write(lines(nl), amt) '  %prompt_string                 = ', s%global%prompt_string
    nl=nl+1; write(lines(nl), amt) '  %print_command                 = ', s%global%print_command
    nl=nl+1; write(lines(nl), amt) '  %random_engine                 = ', s%global%random_engine
    nl=nl+1; write(lines(nl), amt) '  %random_gauss_converter        = ', s%global%random_gauss_converter
    nl=nl+1; write(lines(nl), rmt) '  %delta_e_chrom                 = ', s%global%delta_e_chrom
    nl=nl+1; write(lines(nl), rmt) '  %random_sigma_cutoff           = ', s%global%random_sigma_cutoff
    nl=nl+1; write(lines(nl), imt) '  %random_seed                   = ', s%global%random_seed
    if (s%global%random_seed == 0) then
      call ran_seed_get(ix)
      nl=nl+1; write(lines(nl), imt) '   random_seed (generated)      = ', ix
    endif
    nl=nl+1; write(lines(nl), amt) '  %track_type                    = ', s%global%track_type
    nl=nl+1; write(lines(nl), lmt) '  %var_limits_on                 = ', s%global%var_limits_on
    nl=nl+1; write(lines(nl), amt) '  %var_out_file                  = ', s%global%var_out_file
    nl=nl+1; write(lines(nl), rmt) '  %y_axis_plot_dmin              = ', s%global%y_axis_plot_dmin
    nl=nl+1; write(lines(nl), lmt) '  %draw_curve_off_scale_warn     = ', s%global%draw_curve_off_scale_warn
    nl=nl+1; write(lines(nl), lmt) '  %optimizer_var_limit_warn      = ', s%global%optimizer_var_limit_warn

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Internal Tao Parameters:'
    nl=nl+1; write(lines(nl), imt) 'Universe index range:        = ', lbound(s%u, 1), ubound(s%u, 1)
    nl=nl+1; write(lines(nl), imt) 'default_universe:            = ', s%com%default_universe
    nl=nl+1; write(lines(nl), imt) 'default_branch:              = ', s%com%default_branch
    nl=nl+1; write(lines(nl), lmt) 'common_lattice               = ', s%com%common_lattice
    nl=nl+1; write(lines(nl), amt) 's%com%beam_file              = ', s%com%beam_file
    nl=nl+1; write(lines(nl), amt) 's%com%beam_all_file          = ', s%com%beam_all_file
    nl=nl+1; write(lines(nl), amt) 's%com%beam0_file             = ', s%com%beam0_file
    nl=nl+1; write(lines(nl), amt) 's%com%data_file              = ', s%com%data_file
    nl=nl+1; write(lines(nl), amt) 's%com%init_tao_file          = ', s%com%init_tao_file
    nl=nl+1; write(lines(nl), amt) 's%com%lat_file               = ', s%com%lat_file
    nl=nl+1; write(lines(nl), amt) 's%com%plot_file              = ', s%com%plot_file
    nl=nl+1; write(lines(nl), amt) 's%com%var_file               = ', s%com%var_file
    nl=nl+1; write(lines(nl), amt) 's%com%startup_file           = ', s%com%startup_file
    nl=nl+1; write(lines(nl), lmt) 's%com%combine_consecutive_elements_of_like_name = ', &
                                                s%com%combine_consecutive_elements_of_like_name
    nl=nl+1; write(lines(nl), imt) 'Number paused command files    = ', count(s%com%cmd_file%paused)
  endif

  if (print_ran_state) then
    call ran_default_state (get_state = ran_state)
    nl=nl+1; write(lines(nl), '(a, i0, 2x, i0, l3, es26.16)') 'ran_state = ', ran_state
  endif

  if (print_optimization) then
    call show_opt()
  endif

  if (print_bmad_com) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Bmad_com Parameters:'
    nl=nl+1; write(lines(nl), imt) '  ptc_com%taylor_order_ptc   = ', ptc_com%taylor_order_ptc
    nl=nl+1; write(lines(nl), imt) '  %ptc_max_fringe_order      = ', bmad_com%ptc_max_fringe_order
    nl=nl+1; write(lines(nl), lmt) '  %auto_bookkeeper           = ', bmad_com%auto_bookkeeper
    nl=nl+1; write(lines(nl), lmt) '  %space_charge_on           = ', bmad_com%space_charge_on
    nl=nl+1; write(lines(nl), lmt) '  %coherent_synch_rad_on     = ', bmad_com%coherent_synch_rad_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_tracking_on          = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_damping_on      = ', bmad_com%radiation_damping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_fluctuations_on = ', bmad_com%radiation_fluctuations_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_tracking_on          = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) '  %use_hard_edge_drifts      = ', bmad_com%use_hard_edge_drifts
  endif

  if (print_csr_param) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'CSR_param Parameters:'
    nl=nl+1; write(lines(nl), rmt) '  %ds_track_step        = ', csr_param%ds_track_step
    nl=nl+1; write(lines(nl), rmt) '  %beam_chamber_height  = ', csr_param%beam_chamber_height
    nl=nl+1; write(lines(nl), rmt) '  %sigma_cutoff         = ', csr_param%sigma_cutoff
    nl=nl+1; write(lines(nl), imt) '  %n_bin                = ', csr_param%n_bin
    nl=nl+1; write(lines(nl), imt) '  %particle_bin_span    = ', csr_param%particle_bin_span
    nl=nl+1; write(lines(nl), imt) '  %n_shield_images      = ', csr_param%n_shield_images
    nl=nl+1; write(lines(nl), imt) '  %ix1_ele_csr          = ', csr_param%ix1_ele_csr
    nl=nl+1; write(lines(nl), imt) '  %ix2_ele_csr          = ', csr_param%ix2_ele_csr
    nl=nl+1; write(lines(nl), lmt) '  %lcsr_component_on    = ', csr_param%lcsr_component_on
    nl=nl+1; write(lines(nl), lmt) '  %lsc_component_on     = ', csr_param%lsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %tsc_component_on     = ', csr_param%tsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %small_angle_approx   = ', csr_param%small_angle_approx
  endif

  result_id = show_what

!----------------------------------------------------------------------
! graph

case ('graph')

  ! Look for switches

  call string_trim(stuff, stuff2, ix)

  ! Find particular graph

  call tao_find_plots (err, stuff2, 'BOTH', graph = graph, always_allocate = .true.)
  if (err) return

  if (allocated(graph)) then
    g => graph(1)%g
    if (associated(g%p%r)) then
      nl=nl+1; lines(nl) = 'Region.Graph: ' // trim(g%p%name) // '.' // trim(g%name)
    endif
    nl=nl+1; lines(nl) = 'Plot.Graph:   ' // trim(g%p%name) // '.' // trim(g%name)
    nl=nl+1; write(lines(nl), amt) 'type                  = ', g%type
    nl=nl+1; write(lines(nl), amt) 'title                 = ', g%title
    nl=nl+1; write(lines(nl), amt) 'title_suffix          = ', g%title_suffix
    nl=nl+1; write(lines(nl), amt) 'component             = ', g%component
    nl=nl+1; write(lines(nl), '(a, 4f10.2, 2x, a)') &
                                    'margin                = ', g%margin
    nl=nl+1; write(lines(nl), '(a, 4f10.2, 2x, a)') &
                                    'scale_margin          = ', g%scale_margin
    nl=nl+1; write(lines(nl), imt) 'box                   = ', g%box
    nl=nl+1; write(lines(nl), imt) 'ix_universe           = ', g%ix_universe
    nl=nl+1; write(lines(nl), lmt) 'valid                 = ', g%valid

    nl=nl+1; write(lines(nl), rmt) 'x_axis_scale_factor   = ', g%x_axis_scale_factor
    nl=nl+1; write(lines(nl), rmt) 'symbol_size_scale     = ', g%symbol_size_scale
    nl=nl+1; write(lines(nl), amt) 'x%label               = ', g%x%label
    nl=nl+1; write(lines(nl), rmt) 'x%max                 = ', g%x%max
    nl=nl+1; write(lines(nl), rmt) 'x%min                 = ', g%x%min
    nl=nl+1; write(lines(nl), imt) 'x%major_div           = ', g%x%major_div
    nl=nl+1; write(lines(nl), imt) 'x%major_div_nominal   = ', g%x%major_div_nominal
    nl=nl+1; write(lines(nl), imt) 'x%places              = ', g%x%places
    nl=nl+1; write(lines(nl), lmt) 'x%draw_label          = ', g%x%draw_label
    nl=nl+1; write(lines(nl), lmt) 'x%draw_numbers        = ', g%x%draw_numbers

    nl=nl+1; write(lines(nl), lmt) 'y2_mirrors_y          = ', g%y2_mirrors_y
    nl=nl+1; write(lines(nl), amt) 'y%label               = ', g%y%label
    nl=nl+1; write(lines(nl), rmt) 'y%label_offset        = ', g%y%label_offset
    nl=nl+1; write(lines(nl), rmt) 'y%max                 = ', g%y%max
    nl=nl+1; write(lines(nl), rmt) 'y%min                 = ', g%y%min
    nl=nl+1; write(lines(nl), imt) 'y%major_div           = ', g%y%major_div
    nl=nl+1; write(lines(nl), imt) 'y%major_div_nominal   = ', g%y%major_div_nominal
    nl=nl+1; write(lines(nl), imt) 'y%places              = ', g%y%places
    nl=nl+1; write(lines(nl), lmt) 'y%draw_label          = ', g%y%draw_label
    nl=nl+1; write(lines(nl), lmt) 'y%draw_numbers        = ', g%y%draw_numbers

    nl=nl+1; write(lines(nl), amt) 'y2%label              = ', g%y2%label
    nl=nl+1; write(lines(nl), rmt) 'y2%label_offset       = ', g%y2%label_offset
    nl=nl+1; write(lines(nl), rmt) 'y2%max                = ', g%y2%max
    nl=nl+1; write(lines(nl), rmt) 'y2%min                = ', g%y2%min
    nl=nl+1; write(lines(nl), imt) 'y2%major_div          = ', g%y2%major_div
    nl=nl+1; write(lines(nl), imt) 'y2%major_div_nominal  = ', g%y2%major_div_nominal
    nl=nl+1; write(lines(nl), imt) 'y2%places             = ', g%y2%places
    nl=nl+1; write(lines(nl), lmt) 'y2%draw_label         = ', g%y2%draw_label
    nl=nl+1; write(lines(nl), lmt) 'y2%draw_numbers       = ', g%y2%draw_numbers
    nl=nl+1; write(lines(nl), lmt) 'limited               = ', g%limited
    nl=nl+1; write(lines(nl), lmt) 'clip                  = ', g%clip
    nl=nl+1; write(lines(nl), lmt) 'draw_axes             = ', g%draw_axes
    nl=nl+1; write(lines(nl), lmt) 'draw_grid             = ', g%draw_grid
    nl=nl+1; write(lines(nl), lmt) 'correct_xy_distortion = ', g%correct_xy_distortion
    nl=nl+1; write(lines(nl), lmt) 'draw_only_good_user_data_or_vars = ', g%draw_only_good_user_data_or_vars
    nl=nl+1; lines(nl) = 'Curves:'
    do i = 1, size(g%curve)
      nl=nl+1; write(lines(nl), amt) '   ', g%curve(i)%name
    enddo

  else
    nl=1; lines(1) = 'This is not a graph'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! hom

case ('hom')

  nl=nl+1; lines(nl) = &
        '       #        Freq         R/Q           Q   m  Polarization_Angle'
  do i = 1, ubound(lat%ele, 1)
    ele => lat%ele(i)
    if (ele%key /= lcavity$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    nl=nl+1; write(lines(nl), '(a, i6)') ele%name, i
    do j = 1, size(ele%wake%lr)
      lr => ele%wake%lr(j)
      angle_str = '-'
      if (lr%polarized) write (angle_str, '(f9.4)') lr%angle
      nl=nl+1; write(lines(nl), '(i8, 3es12.4, i4, a)') j, &
                  lr%freq, lr%R_over_Q, lr%Q, lr%m, angle_str
    enddo
    nl=nl+1; lines(nl) = ' '
  enddo
  nl=nl+1; lines(nl) = '       #        Freq         R/Q           Q   m  Polarization_Angle'

  result_id = show_what

!----------------------------------------------------------------------
! keys

case ('key_bindings')

  call tao_key_info_to_str (1, 1, size(s%key), str, header)
  nl=nl+1; lines(nl) = ' Ix  ' // header

  do i = 1, size(s%key)
    call tao_key_info_to_str (i, 1, size(s%key), str, header)
    nl=nl+1; write(lines(nl), '(i3, 2x, a)') i, str
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! lattice

case ('lattice')
  
  limited = .false.
  all_lat = .false.
  at_ends = .true.
  by_s = .false.
  print_header_lines = .true.
  print_tail_lines = .true.
  replacement_for_blank = ''
  ix_branch = s%com%default_branch
  undef_str = '---'
  show_lords = .false.
  what_to_print = 'standard'
  ix_remove = -1
  lat_type = model$

  column(:)%name = ""
  column(:)%label = ""
  column(:)%remove_line_if_zero = .false.

  ! get command line switches

  do
    call tao_next_switch (stuff2, [ &
        '-branch             ', '-blank_replacement  ', '-lords              ', '-middle             ', &
        '-all_tracking       ', '-0undef             ', '-no_label_lines     ', '-no_tail_lines      ', &
        '-custom             ', '-s                  ', '-radiation_integrals', '-remove_line_if_zero', &
        '-base               ', '-design             ', '-floor_coords       '], &
              switch, err, ix_s2)
    if (err) return
    if (switch == '') exit
    select case (switch)

    case ('-0undef')
      undef_str = '  0'

    case ('-all_tracking')
      all_lat = .true. 

    case ('-base')
      lat_type = base$

    case ('-blank_replacement')
      replacement_for_blank = stuff2(1:ix_s2)
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-branch')
      branch => pointer_to_branch(stuff2(1:ix_s2), u%model%lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif
      ix_branch = branch%ix_branch
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-custom')
      what_to_print = 'custom'
      file_name = stuff2(1:ix_s2)
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)
      iu = lunget()
      open (iu, file = file_name, status = 'old', iostat = ios)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT OPEN FILE: ' // file_name
        return
      endif
      column(:)%name = ""
      column(:)%label = ""
      read (iu, nml = custom_show_list, iostat = ios)
      close (iu)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ "CUSTOM_SHOW_LIST" NAMELIST IN FILE: ' // file_name
        return
      endif

    case ('-design')
      lat_type = design$

    case ('-floor_coords')
      what_to_print = 'floor_coords'

    case ('-lords')
      show_lords = .true.

    case ('-middle')
      at_ends = .false.

    case ('-no_label_lines')
      print_header_lines = .false.
      print_tail_lines = .false.

    case ('-no_tail_lines')
      print_tail_lines = .false.

    case ('-radiation_integrals')
      what_to_print = 'rad_int'

    case ('-remove_line_if_zero')
      read (stuff2(1:ix_s2), *, iostat = ios) ix_remove
      if (ios /= 0 .or. ix_remove < 1 .or. ix_remove > size(column)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-remove_line_if_zero" argument'
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-s')
      by_s = .true.
    end select

  enddo

  !

  tao_lat => tao_pointer_to_tao_lat(u, lat_type)
  lat => tao_lat%lat
  branch => lat%branch(ix_branch)

  ! Construct columns if needed.

  select case (what_to_print)
  case ('floor_coords')
    column( 1)  = show_lat_column_struct('#',                      'i6',        6, '', .false.)
    column( 2)  = show_lat_column_struct('x',                      'x',         2, '', .false.)
    column( 3)  = show_lat_column_struct('ele::#[name]',           'a',         0, '', .false.)
    column( 4)  = show_lat_column_struct('ele::#[key]',            'a16',      16, '', .false.)
    column( 5)  = show_lat_column_struct('ele::#[s]',              'f10.3',    10, '', .false.)
    column( 6)  = show_lat_column_struct('ele::#[x_position]',     'f12.5',    12, 'X',     .false.)
    column( 7)  = show_lat_column_struct('ele::#[y_position]',     'f12.5',    12, 'Y',     .false.)
    column( 8)  = show_lat_column_struct('ele::#[z_position]',     'f12.5',    12, 'Z',     .false.)
    column( 9)  = show_lat_column_struct('ele::#[theta_position]', 'f12.5',    12, 'Theta', .false.)
    column(10)  = show_lat_column_struct('ele::#[phi_position]',   'f12.5',    12, 'Phi',   .false.)
    column(11)  = show_lat_column_struct('ele::#[psi_position]',   'f12.5',    12, 'Psi',   .false.)

  case ('rad_int')
    column(1)  = show_lat_column_struct('#',                     'i6',        6, '', .false.)
    column(2)  = show_lat_column_struct('x',                     'x',         2, '', .false.)
    column(3)  = show_lat_column_struct('ele::#[name]',          'a',         0, '', .false.)
    column(4)  = show_lat_column_struct('ele::#[key]',           'a16',      16, '', .false.)
    column(5)  = show_lat_column_struct('ele::#[s]',             'f10.3',    10, '', .false.)
    if (branch%param%geometry == open$) then
      column(6)  = show_lat_column_struct('lat::rad_int1.i1[#]',     'es10.2',  10, '', .true.)
      column(7)  = show_lat_column_struct('lat::rad_int1.i2_e4[#]',  'es10.2',  10, '', .false.)
      column(8)  = show_lat_column_struct('lat::rad_int1.i3_e7[#]',  'es10.2',  10, '', .false.)
      column(9)  = show_lat_column_struct('lat::rad_int1.i5a_e6[#]', 'es10.2',  10, '', .false.)
      column(10) = show_lat_column_struct('lat::rad_int1.i5b_e6[#]', 'es10.2',  10, '', .false.)
    else
      column(6)  = show_lat_column_struct('lat::rad_int1.i1[#]',     'es10.2',  10, '', .true.)
      column(7)  = show_lat_column_struct('lat::rad_int1.i2[#]',     'es10.2',  10, '', .false.)
      column(8)  = show_lat_column_struct('lat::rad_int1.i3[#]',     'es10.2',  10, '', .false.)
      column(9)  = show_lat_column_struct('lat::rad_int1.i4a[#]',    'es10.2',  10, '', .false.)
      column(10) = show_lat_column_struct('lat::rad_int1.i5a[#]',    'es10.2',  10, '', .false.)
      column(11) = show_lat_column_struct('lat::rad_int1.i4b[#]',    'es10.2',  10, '', .false.)
      column(12) = show_lat_column_struct('lat::rad_int1.i5b[#]',    'es10.2',  10, '', .false.)
      column(13) = show_lat_column_struct('lat::rad_int1.i6b[#]',    'es10.2',  10, '', .false.)
    endif

  case ('standard')
    column(1)  = show_lat_column_struct('#',                 'i6',        6, '', .false.)
    column(2)  = show_lat_column_struct('x',                 'x',         2, '', .false.)
    column(3)  = show_lat_column_struct('ele::#[name]',      'a',         0, '', .false.)
    column(4)  = show_lat_column_struct('ele::#[key]',       'a16',      16, '', .false.)
    column(5)  = show_lat_column_struct('ele::#[s]',         'f10.3',    10, '', .false.)
    column(6)  = show_lat_column_struct('ele::#[l]',         'f8.3',      8, '', .false.)
    column(7)  = show_lat_column_struct('ele::#[beta_a]',    'f8.2',      8, 'beta|  a', .false.)
    column(8)  = show_lat_column_struct('ele::#[phi_a]',     'f8.3',      8, '', .false.)
    column(9)  = show_lat_column_struct('ele::#[eta_a]',     'f5.1',      5, '', .false.)
    column(10) = show_lat_column_struct('ele::#[orbit_x]',   '3p, f8.3',  8, '', .false.)
    column(11) = show_lat_column_struct('ele::#[beta_b]',    'f8.2',      8, '', .false.)
    column(12) = show_lat_column_struct('ele::#[phi_b]',     'f8.3',      8, '', .false.)
    column(13) = show_lat_column_struct('ele::#[eta_b]',     'f5.1',      5, '', .false.)
    column(14) = show_lat_column_struct('ele::#[orbit_y]',   '3p, f8.3',  8, '', .false.)
  end select

  if (what_to_print /= 'custom' .and. show_lords) then
    n = size(column)
    column(8:n) = column(6:n-2)
    column(6)  = show_lat_column_struct('x',                   'x',         2, '', .false.)
    column(7)  = show_lat_column_struct("ele::#[lord_status]", 'a16',      16, '', .false.) 
  endif


  ! remove_line_if_zero bookkeeping. Ignore space lines (name = 'x')

  if (ix_remove > 0) then
    j = 0
    do i = 1, size(column)
      if (column(i)%name == 'x') cycle
      j = j + 1
      if (j == ix_remove) then
        column(i)%remove_line_if_zero = .true.
        exit
      endif
      if (i == size(column)) then
        nl=1; lines(1) = 'ARGUMENT FOR "-remove_line_if_zero" OUT OF RANGE!'
        return
      endif
    enddo
  endif

  ! Need to compute radiation integrals?

  do i = 1, size(column)
    if (index(column(i)%name, 'rad_int') /= 0) then
      ix = s%com%default_branch
      call radiation_integrals (u%model%lat, u%model%lat_branch(ix)%orbit, u%model%modes, u%model%ix_rad_int_cache, 0, u%model%rad_int)
      call radiation_integrals (u%design%lat, u%design%lat_branch(ix)%orbit, u%design%modes, u%design%ix_rad_int_cache, 0, u%design%rad_int)
      call radiation_integrals (u%base%lat, u%base%lat_branch(ix)%orbit, u%base%modes, u%base%ix_rad_int_cache, 0, u%base%rad_int)
      exit
    endif
  enddo

  ! Find elements to use

  if (allocated (picked_ele)) deallocate (picked_ele)
  allocate (picked_ele(0:branch%n_ele_max))
  if (show_lords) then
    picked_ele(0:branch%n_ele_track) = .false.
    picked_ele(branch%n_ele_track+1:branch%n_ele_max) = .true.
  else
    picked_ele(0:branch%n_ele_track) = .true.
    picked_ele(branch%n_ele_track+1:branch%n_ele_max) = .false.
  endif

  if (by_s) then
    ix_s2 = index(stuff2, ':')
    if (ix_s2 == 0) then
      nl=1; lines(1) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (stuff2(1:ix_s2-1), *, iostat = ios1) s1
    read (stuff2(ix_s2+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // stuff2
      return
    endif

    picked_ele = .false.
    do ie = 1, branch%n_ele_track
      if (at_ends) then;
        s_ele = branch%ele(ie)%s
      else
        s_ele = (branch%ele(ie-1)%s + branch%ele(ie)%s) / 2
      endif
      if (s_ele >= s1 .and. s_ele <= s2) picked_ele(ie) = .true.
    enddo

  elseif (stuff2(1:ix_s2) == '*' .or. all_lat) then
    ! picked_ele already set

  elseif (ix_s2 /= 0) then
    call tao_locate_elements (stuff2, u%ix_uni, eles, err, lat_type, &
                  ignore_blank = .true., above_ubound_is_err = .false., ix_dflt_branch = s%com%default_branch)
    if (err) return
    picked_ele = .false.
    do i = 1, size(eles)
      picked_ele(eles(i)%ele%ix_ele) = .true.
    enddo

  elseif (.not. show_lords) then
    if (count(picked_ele) > 300) then
      picked_ele(201:) = .false.
      limited = .true.
    endif

  endif

  !

  if (at_ends) then
    write (line1, '(6x, a)') 'Values at End of Element:'
  else
    write (line1, '(6x, a)') 'Values at Center of Element:'
  endif

  ! Setup columns

  ix1 = 1
  line2 = ""
  line3 = ""
  do i = 1, size(column)
    if (column(i)%name == "") cycle

    ! Convert from old 'dat::' format to 'lat::' format.
    ix = index(column(i)%name, 'dat::')
    if (ix /= 0) column(i)%name = column(i)%name(1:ix-1) // 'lat' // column(i)%name(ix+3:)

    column(i)%format = '(' // trim(column(i)%format) // ')'

    if (column(i)%field_width == 0) then
      if (column(i)%name /= 'ele::#[name]') then
        call out_io (s_error$, r_name, &
            'FIELD_WIDTH = 0 CAN ONLY BE USED WITH "ele::#[name]" TYPE COLUMNS')
        return
      endif
      column(i)%field_width = 5
      do ie = 0, branch%n_ele_max
        if (.not. picked_ele(ie)) cycle
        column(i)%field_width = max(column(i)%field_width, len_trim(branch%ele(ie)%name)+1)
      enddo
    endif

    ix2 = ix1 + column(i)%field_width

    if (column(i)%label == '') then
      name = column(i)%name
      if (index(name, 'ele::') /= 0) then
        i1 = index(name, '[')
        i2 = index(name, ']')
        name = name(i1+1:i2-1)
      elseif (index(name, 'beam::') /= 0) then
        ix = index(name, 'beam::')
        name = name(ix+6:)
        i2 = index(name, '[')
        name = name(1:i2-1)
      elseif (index(name, 'lat::') /= 0) then
        ix = index(name, 'lat::')
        name = name(ix+5:)
        i2 = index(name, '[')
        name = name(1:i2-1)
      elseif  (name == '#') then
        line2(ix2-5:) = 'Index'
        ix1 = ix2
        cycle    
      elseif (name == 'x') then
        ix1 = ix2
        cycle    
      else
        name = ''
      endif

      ix = index(name, '.')
      if (ix == 0) ix = index(name, '_')
      n = len_trim(name)
      if (column(i)%format(2:2) == 'a') then
        line2(ix1:) = name
      elseif (ix == 0) then
        line2(ix2-n:) = name
      else
        line2(ix2-ix+1:) = name(1:ix-1)
        line3(ix2-n+ix:) = name(ix+1:)
      endif

    else
      name = column(i)%label
      ix = index(name, '|')
      if (ix == 0) then
        j = len_trim(name)
        line2(ix2-j:) = name(1:j)
      else
        j = max(ix-1, len_trim(name(ix+1:)))
        line2(ix2-j:) = name(1:ix-1)
        line3(ix2-j:) = trim(name(ix+1:))
      endif
    endif

    ix1 = ix2
  enddo

  ! Collect lines

  if (print_header_lines) then
    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
  endif

  line_loop: do ie = 0, branch%n_ele_max
    if (.not. picked_ele(ie)) cycle
    if (size(lines) < nl+100) call re_allocate (lines, nl+200, .false.)
    line = ''
    nc = 1
    ele => branch%ele(ie)
    do i = 1, size(column)
        
      name = column(i)%name
      if (name == '') cycle

      if (name == '#') then
        write (line(nc:), column(i)%format, iostat = ios) ie

      elseif (name == 'ele::#[name]') then
        write (line(nc:), column(i)%format, iostat = ios) ele%name

      elseif (name == 'ele::#[key]') then
        write (line(nc:), column(i)%format, iostat = ios) key_name(ele%key)

      elseif (name == 'ele::#[slave_status]') then
        write (line(nc:), column(i)%format, iostat = ios) control_name(ele%slave_status)

      elseif (name == 'ele::#[lord_status]') then
        write (line(nc:), column(i)%format, iostat = ios) control_name(ele%lord_status)

      elseif (name == 'ele::#[type]') then
        if (ele%type == '') then
          write (line(nc:), column(i)%format, iostat = ios) replacement_for_blank
        else
          write (line(nc:), column(i)%format, iostat = ios) ele%type
        endif

      elseif (name == 'x') then
        ios = 0

      else
        write (nam, '(i0, a, i0)') ix_branch, '>>', ie
        call str_substitute (name, '#', trim(nam))
        ix = index(name, 'ele::')
        if (.not. at_ends .and. ix /= 0) then
          name = name(:ix+2) // '_mid' // trim(name(ix+3:))
        endif
        call tao_evaluate_expression (name, 1, .false., value, good, err, .false., &
                                                  dflt_component = lat_type_name(lat_type))
        if (err .or. .not. allocated(value) .or. size(value) /= 1) then
          if (column(i)%remove_line_if_zero) cycle line_loop
          n = len(undef_str)
          k = min(n, column(i)%field_width - 1)
          j = nc + column(i)%field_width - k
          line(j:) = undef_str(n-k+1:n)
        else
          if (index(column(i)%format, 'l') /= 0 .or. index(column(i)%format, 'L') /= 0) then
            if (value(1) == 0) then
              write (line(nc:), column(i)%format, iostat = ios) .false.
            else
              write (line(nc:), column(i)%format, iostat = ios) .true.
            endif
          elseif (index(column(i)%format, 'i') /= 0 .or. index(column(i)%format, 'I') /= 0) then
            write (line(nc:), column(i)%format, iostat = ios) nint(value(1))
            if (column(i)%remove_line_if_zero .and. nint(value(1)) == 0) cycle line_loop
          else
            call write_real (line(nc:), column(i)%format, value(1))
            if (column(i)%remove_line_if_zero .and. value(1) == 0) cycle line_loop
          endif
        endif
      endif

      if (ios /= 0) then
        lines(1) = 'WIDTH TOO SMALL FOR NUMBER OR BAD FORMAT: ' // column(i)%format
        lines(2) = 'FOR DISPLAYING: ' // column(i)%name
        nl = 2
        return
      endif

      nc  = nc + column(i)%field_width

    enddo

    nl=nl+1; lines(nl) = line

  enddo line_loop

  if (print_tail_lines) then
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
    nl=nl+1; lines(nl) = line1
  endif

  if (limited) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, i0)') &
          'NOTE: Since no range given, the number of elements shown is first 200 of ', branch%n_ele_track
  endif

  deallocate(picked_ele)

  result_id = show_what

!----------------------------------------------------------------------
! normal_form

case ('normal_form')
  
  if (s%global%rf_on) then
    normal_form => branch%normal_form_with_rf
  else
    normal_form => branch%normal_form_no_rf
  endif
  
  if (.not. associated(normal_form%ele_origin) ) then
    nl=nl+1; lines(nl) ='One-turn-map has not been computed'
    return
  endif

  call tao_next_switch (stuff2, ['-order'], switch, err, ix)
  if (err) return
  
  n_order = ptc_com%taylor_order_ptc
  select case (switch)
    case ('-order')
      read (stuff2(:ix), *, iostat = ios) n_order
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ORDER NUMBER!'
        return
      endif
    call string_trim (stuff2(ix+1:), stuff2, ix)
  end select
  
  nl=nl+1; lines(nl) = 'normal_form: '//stuff2(1:5)
  
  select case(stuff2(1:5))
    case ('dhdj ')
      call type2_taylors (normal_form%dhdj, lines, nl, max_order = n_order)
    case ('A    ')
      call type2_taylors (normal_form%A, lines, nl, max_order = n_order)
    case ('A_inv')
      call type2_taylors (normal_form%A_inv, lines, nl, max_order = n_order)
    case ('M    ')
      call type2_taylors (normal_form%M, lines, nl, max_order = n_order)
    case ('F  ')
      call type2_complex_taylors (normal_form%F, lines, nl, max_order = n_order)
    case ('L  ')
      call type2_complex_taylors (normal_form%L, lines, nl, max_order = n_order)      
    case default 
      nl=nl+1; lines(nl) = 'bad normal_form map: '//trim(stuff2)
      nl=nl+1; lines(nl) = 'Must be one of: M A A_inv dhdj F L'
  end select

  result_id = show_what

!----------------------------------------------------------------------
! optimizer

case ('optimizer')

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    nl=nl+1; lines(nl) = 'Data Used:'
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe: ', i
    endif
    do j = 1, u%n_d2_data_used
      if (u%d2_data(j)%name == ' ') cycle
      call tao_data_show_use (u%d2_data(j), lines, nl)
    enddo
  enddo

  nl=nl+1; lines(nl) = 'Variables Used:'
  do j = 1, s%n_v1_var_used
    if (s%v1_var(j)%name == ' ') cycle
    call tao_var_show_use (s%v1_var(j), lines, nl)
  enddo

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; write(lines(nl), amt) 'optimizer:        ', s%global%optimizer
  call show_opt
  call out_io (s_blank$, r_name, lines(1:nl))
  nl = 0

  result_id = show_what

!----------------------------------------------------------------------
! particle

case ('orbit')

  call tao_locate_elements (word1, u%ix_uni, eles, err, ix_dflt_branch = s%com%default_branch)
  if (err) return
  do i = 1, 6
    nl=nl+1; write(lines(nl), rmt) '     ', &
                u%model%lat_branch(eles(1)%ele%ix_branch)%orbit%vec(eles(1)%ele%ix_ele)
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! particle

case ('particle')

  nb = s%global%bunch_to_plot
  ix_branch = s%com%default_branch
  show_all = .false.
  show_lost = .false.
  ele_name = ''
  ix_ele = 0

  do

    call tao_next_switch (stuff2, ['-element ', '-particle', '-bunch   ', '-lost    ', '-all     '], switch, err, ix_word)
    if (err) return

    select case (switch)
    case ('') 
      exit

    case ('-lost') 
      show_lost = .true.

    case ('-all')
      show_all = .true.

    case ('-element')
      ele_name = stuff2(:ix_word)
      call string_trim (stuff2(ix_word+1:), stuff2, ix_word)

      call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
      if (err) return
      call tao_locate_elements (ele_name, ix_u, eles, err, ix_dflt_branch = s%com%default_branch)
      if (err) return
      ix_ele = eles(1)%ele%ix_ele
      ix_branch = eles(1)%ele%ix_branch

    case ('-particle')
      read (stuff2(:ix_word), *, iostat = ios) ix_p
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ PARTICLE INDEX!'
        return
      endif
      call string_trim (stuff2(ix_word+1:), stuff2, ix_word)

    case ('-bunch')
      read (stuff2(:ix_word), *, iostat = ios) nb
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ BUNCH INDEX!'
        return
      endif
      call string_trim (stuff2(ix_word+1:), stuff2, ix_word)

    end select

  enddo

  uni_branch => u%uni_branch(ix_branch)
  branch => u%model%lat%branch(ix_branch)

  ! show lost

  if (show_lost) then
    bunch => u%uni_branch(ix_branch)%ele(lat%n_ele_track)%beam%bunch(nb)
    nl=nl+1; write(lines(nl), *) 'Bunch:', nb
    nl=nl+1; lines(nl) = 'Particles lost at:'
    nl=nl+1; lines(nl) = '    Ix Ix_Ele  Ele_Name '
    do i = 1, size(bunch%particle)
      if (bunch%particle(i)%state == alive$) cycle
      if (nl == size(lines)) call re_allocate (lines, nl+100, .false.)
      ie = bunch%particle(i)%ix_ele
      nl=nl+1; write(lines(nl), '(i6, i7, 2x, a)') i, ie, lat%ele(ie)%name
    enddo
    result_id = 'particle:lost'
    return
  endif

  ! check

  if (.not. allocated(uni_branch%ele(ix_ele)%beam%bunch)) then
    call out_io (s_error$, r_name, 'BUNCH NOT ASSOCIATED WITH THIS ELEMENT.')
    return
  endif

  if (nb < 1 .or. nb > size(uni_branch%ele(ix_ele)%beam%bunch)) then
    call out_io (s_error$, r_name, 'BUNCH INDEX OUT OF RANGE: \i0\ ', i_array = [ nb ])
    return
  endif

  bunch => uni_branch%ele(ix_ele)%beam%bunch(nb)

  ! show all

  if (show_all) then
    nl=nl+1; write(lines(nl), *) 'Element:', ix_ele, '  ', branch%ele(ix_ele)%name
    nl=nl+1; write(lines(nl), '(a, 6(12x, a))') '  Ix', '  x', 'px', '  y', 'py', '  z', 'pz'
    do i = 1, size(bunch%particle)
      if (nl == size(lines)) call re_allocate (lines, nl+100, .false.)
      nl=nl+1; write(lines(nl), '(i6, 6es15.7)') i, (bunch%particle(i)%vec(j), j = 1, 6)
    enddo
    result_id = 'particle:lost'
    return
  endif

  !

  if (ix_p < 1 .or. ix_p > size(bunch%particle)) then
    call out_io (s_error$, r_name, 'PARTICLE INDEX OUT OF RANGE: \i0\ ', i_array = [ ix_p ])
    return
  endif

  nl=nl+1; write(lines(nl), imt) 'At lattice element:', ix_ele
  nl=nl+1; write(lines(nl), imt) 'Bunch:       ', nb
  nl=nl+1; write(lines(nl), imt) 'Particle:    ', ix_p
  nl=nl+1; write(lines(nl), lmt) 'Is Alive?    ', bunch%particle(ix_p)%state == alive$
  if (u%model%lat%branch(ix_branch)%param%particle == photon$) then
    nl=nl+1; write(lines(nl), rmt) 'Intensity_x: ', bunch%particle(ix_p)%field(1)**2
    nl=nl+1; write(lines(nl), rmt) 'Intensity_y: ', bunch%particle(ix_p)%field(2)**2
  else
    nl=nl+1; write(lines(nl), rmt) 'Charge:      ', bunch%particle(ix_p)%charge
  endif
  nl=nl+1; write(lines(nl), lmt) 'Coords: '
  nl=nl+1; write(lines(nl), '(a, 6es13.5)') '  ', bunch%particle(ix_p)%vec

  result_id = show_what

!----------------------------------------------------------------------
! plot

case ('plot')

  ! Look for switches

  what = ''

  do
    call tao_next_switch (stuff2, ['-floor_plan', '-lat_layout'], switch, err, ix)
    if (err) return
    if (switch == '') exit
    what = switch
  enddo

  ! Floor plan info

  if (what == '-floor_plan') then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Plot_page parameters:'
    nl=nl+1; write(lines(nl), f3mt) '%floor_plan_rotation        = ', s%plotting%floor_plan_rotation
    nl=nl+1; write(lines(nl), amt)  '%floor_plan_view            = ', s%plotting%floor_plan_view
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Shape_Name  Ele_Name                        Shape         Color           Size  Label  Draw'
    nl=nl+1; lines(nl) = &
          '----------  ----------------------------    --------      -----           ----  -----  ----'

    do i = 1, size(s%plotting%floor_plan%ele_shape)
      shape => s%plotting%floor_plan%ele_shape(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write(lines(nl), '(a, i0, t13, 3a, f10.1, 2x, a6, 1x, l2, 4x, a)') &
                'shape', i, shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%size, shape%label, shape%draw
    enddo

    result_id = 'plot:floor_plan'
    return
  endif

  ! lat_layout info

  if (what == '-lat_layout') then
    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Shape_Name  Ele_Name                        Shape         Color           Size  Label  Draw'
    nl=nl+1; lines(nl) = &
          '----------  ----------------------------    --------      -----           ----  -----  ----'

    do i = 1, size(s%plotting%lat_layout%ele_shape)
      shape => s%plotting%lat_layout%ele_shape(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write(lines(nl), '(a, i0, t13, 3a, f10.1, 2x, a6, 1x, l2, 4x, a)') &
                'shape', i, shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%size, shape%label, shape%draw
    enddo

    result_id = 'plot:lat_layout'
    return 
  endif

  ! stuff2 is blank => print overall info

  if (stuff2 == ' ') then

    nl=nl+1; lines(nl) = 'plot_page parameters:'
    nl=nl+1; write(lines(nl), imt)  '%size                       = ', nint(s%plotting%size)
    nl=nl+1; write(lines(nl), imt)  '%n_curve_pts                = ', s%plotting%n_curve_pts
    nl=nl+1; write(lines(nl), f3mt) '%text_height                = ', s%plotting%text_height 
    nl=nl+1; write(lines(nl), f3mt) '%main_title_text_scale      = ', s%plotting%main_title_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%graph_title_text_scale     = ', s%plotting%graph_title_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%axis_number_text_scale     = ', s%plotting%axis_number_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%axis_label_text_scale      = ', s%plotting%axis_label_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%key_table_text_scale       = ', s%plotting%key_table_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%legend_text_scale          = ', s%plotting%legend_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '%floor_plan_rotation        = ', s%plotting%floor_plan_rotation
    nl=nl+1; write(lines(nl), amt)  '%floor_plan_view            = ', s%plotting%floor_plan_view

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Templates:'
    nl=nl+1; lines(nl) = '   Plot                          .Graph                   Description        '
    nl=nl+1; lines(nl) = '   ----------------------------  -----------------------  -------------------'
    do i = 1, size(s%plotting%template)
      p => s%plotting%template(i)
      if (p%name == '') cycle
      if (p%name == 'scratch') cycle
      if (allocated(p%graph)) then
          write (fmt, '(a, i0, a)') '(3x, a30, ', size(p%graph), '(2a, 2x), t57, a)'
          nl=nl+1; write(lines(nl), fmt) p%name, &
                      ('.', trim(p%graph(j)%name), j = 1, size(p%graph)), p%description
      else
        nl=nl+1; write(lines(nl), '(3x, a)') p%name 
      endif
    enddo

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = '                    ' // &
                         '                                    Location on Page'
    nl=nl+1; lines(nl) = 'Visible  Plot Region' // &
                         '         <-->  Template             x1    x2    y1    y2'  
    nl=nl+1; lines(nl) = '-------  -----------' // &
                         '               -----------------------------------------'
    do i = 1, size(s%plotting%region)
      region => s%plotting%region(i)
      if (region%name == '') cycle
      nl=nl+1; write(lines(nl), '(3x l1, 5x, a20, a, a18, 4f6.2)') region%visible, &
                                    region%name, '<-->  ', region%plot%name, region%location
    enddo

    result_id = 'plot:'
    return
  endif

! Find particular plot

  call tao_find_plots (err, stuff2, 'BOTH', plot, print_flag = .false.)
  if (err) return

  if (allocated(plot)) then
    p => plot(1)%p
    if (associated(p%r)) then
      nl=nl+1; lines(nl) = 'Region:  ' // trim(p%name)
    endif
    nl=nl+1; lines(nl) = 'Plot:  ' // p%name
    nl=nl+1; write(lines(nl), amt) 'x_axis_type          = ', p%x_axis_type
    nl=nl+1; write(lines(nl), amt) 'x%label              = ', p%x%label
    nl=nl+1; write(lines(nl), rmt) 'x%max                = ', p%x%max
    nl=nl+1; write(lines(nl), rmt) 'x%min                = ', p%x%min
    nl=nl+1; write(lines(nl), imt) 'x%major_div          = ', p%x%major_div
    nl=nl+1; write(lines(nl), imt) 'x%major_div_nominal  = ', p%x%major_div_nominal
    nl=nl+1; write(lines(nl), imt) 'x%places             = ', p%x%places
    nl=nl+1; write(lines(nl), lmt) 'x%draw_label         = ', p%x%draw_label
    nl=nl+1; write(lines(nl), lmt) 'x%draw_numbers       = ', p%x%draw_numbers
    nl=nl+1; write(lines(nl), lmt) 'autoscale_x          = ', p%autoscale_x
    nl=nl+1; write(lines(nl), lmt) 'autoscale_y          = ', p%autoscale_y
    nl=nl+1; write(lines(nl), lmt) 'autoscale_gang_x     = ', p%autoscale_gang_x
    nl=nl+1; write(lines(nl), lmt) 'autoscale_gang_y     = ', p%autoscale_gang_y
    
    nl=nl+1; lines(nl) = 'Graphs:'
    do i = 1, size(p%graph)
      nl=nl+1; write(lines(nl), amt) '   ', p%graph(i)%name
    enddo

  else
    nl=1; lines(1) = 'This is not a name of a plot'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! taylor_map

case ('taylor_map', 'matrix')

  by_s = .false.
  if (show_what == 'matrix') then
    n_order = 1
  else
    n_order = -1
  endif

  do
    call tao_next_switch (stuff2, ['-order', '-s    '], switch, err, ix)
    if (err) return
    if (switch == '') exit
    select case (switch)
    case ('-s')
      by_s = .true.
    case ('-order')
      read (stuff2(:ix), *, iostat = ios) n_order
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ORDER NUMBER!'
        return
      endif
      call string_trim (stuff2(ix+1:), stuff2, ix)
      if (n_order > ptc_com%taylor_order_ptc) then
        nl=1; write(lines(nl), '(a, i0)') &
                  'TAYLOR ORDER CANNOT BE ABOVE ORDER USED IN CALCULATIONS WHICH IS \i0\ ', &
                  ptc_com%taylor_order_ptc
        return
      endif        
    end select
  enddo

  ele1_name = stuff2(:ix)
  call string_trim(stuff2(ix+1:), stuff2, ix)
  ele2_name = stuff2(:ix)
  if (stuff2(ix+1:) /= '') then
    nl=1; lines(1) = 'EXTRA STUFF ON LINE!'
    return
  endif

  ! By s

  if (by_s) then
    if (ele1_name == '') then
      s2 = lat%ele(lat%n_ele_track)%s
      s1 = 0
    else
      read (ele1_name, *, iostat = ios) s1
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD S1 VALUE:' // ele1_name
        return
      endif
    endif

    if (ele2_name == '') then
      s2 = s1
      if (lat%param%geometry == open$) s1 = 0
    else 
      read (ele2_name, *, iostat = ios) s2
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD S2 VALUE:' // ele2_name
        return
      endif
    endif

    if (n_order > 1) then
      call transfer_map_from_s_to_s (lat, taylor, s1, s2, ix_branch, one_turn = .true.)
    else
      call twiss_and_track_at_s (lat, s1, ele0, u%model%lat_branch(ix_branch)%orbit, orb, ix_branch)
      call mat6_from_s_to_s (lat, mat6, vec0, s1, s2, orb, ix_branch, one_turn = .true.)
    endif

  ! By element

  else

    if (ele1_name == '') then
      ix2 = lat%n_ele_track
      ix1 = 0
    else
      call tao_locate_elements (ele1_name, u%ix_uni, eles, err, ix_dflt_branch = s%com%default_branch)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele1_name
        return
      endif
      if (err .or. size(eles) == 0) return
      ix1 = eles(1)%ele%ix_ele
    endif

    if (ele2_name == '') then
      ix2 = ix1
      if (lat%param%geometry == open$) then
        ix1 = 0
      else
        ix1 = ix1 - 1
      endif
    else
      call tao_locate_elements (ele2_name, u%ix_uni, eles, err, ix_dflt_branch = s%com%default_branch)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele2_name
        return
      endif
      if (err .or. size(eles) == 0) return
      ix2 = eles(1)%ele%ix_ele
    endif

    if (n_order > 1) then
      call transfer_map_calc (lat, taylor, ix1, ix2, one_turn = .true.)
    else
      call transfer_matrix_calc (lat, .true., mat6, vec0, ix1, ix2, ix_branch, one_turn = .true.)
    endif

  endif

  ! Print results

  if (n_order > 1) call truncate_taylor_to_order (taylor, n_order, taylor)

  if (n_order > 1) then
    call type2_taylors (taylor, lines, nl)
  else
    vec_in = 0
    if (n_order == 0) then 
      nl = nl+1; write(lines(nl), '(6f11.6)') vec0
    else
      if (any(abs(mat6(1:n_order,1:n_order)) >= 1000)) then
        fmt = '(6es12.4, a, es12.4)'
      else
        fmt = '(6f12.5, a, es12.4)'
      endif

      do i = 1, 6
        nl=nl+1; write(lines(nl), fmt) mat6(i,:), '   : ', vec0(i)
      enddo
    endif
  endif

  result_id = show_what

!----------------------------------------------------------------------
! top10

case ('top10')

  call string_trim(stuff, stuff2, ix)
  if (ix == 0) then
    call tao_show_constraints (0, 'TOP10')
    call tao_top10_merit_categories_print (0)
  elseif (index('-derivative', trim(stuff2)) == 1) then 
    call tao_top10_derivative_print ()
  else
    nl=1; lines(1) = 'UNKNOWN SWITCH: ' // stuff2
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! tune

case ('tune')

  nl=nl+1; lines(nl) = 'Use "show universe" instead.'

  result_id = show_what

!----------------------------------------------------------------------
! twiss
    
case ('twiss_and_orbit')

  tao_lat => u%model
  branch => tao_lat%lat%branch(s%com%default_branch)
  lat_type = model$

  do 

    call tao_next_switch (stuff2, [ &
        '-branch     ', '-universe   ', '-design     ', '-base       '], &
              switch, err, ix_s2)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-base')
      lat_type = base$

    case ('-branch')
      branch => pointer_to_branch(stuff2(1:ix_s2), lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-design')
      lat_type = design$

    case ('-universe')
      read (stuff2(1:ix_s2), *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (ix_s2 == 0 .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" argument'
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    end select

  enddo

  !

  tao_lat => tao_pointer_to_tao_lat (u, lat_type)
  lat => tao_lat%lat
  branch => lat%branch(branch%ix_branch)
  ix_branch = branch%ix_branch

  call string_trim(stuff2, stuff2, ix)
  if (ix == 0) then
    s_pos = 0
  else
    if (.not. is_real(stuff2)) then
      nl=1; lines(1) = 'NOT A REAL NUMBER: ' // stuff2
      return
    endif
    read (stuff2, *) s_pos
  endif

  call twiss_and_track_at_s (lat, s_pos, ele0, tao_lat%lat_branch(ix_branch)%orbit, orb, ix_branch, err)
  if (err) return 

  nl=nl+1; write(lines(nl), '(a, f10.5)') 'At S =', s_pos
  nl=nl+1; write(lines(nl), '(2a)')       'In Element: ', ele0%name

  call type2_twiss (ele0, lines(nl+1:), n, s%global%phase_units)
  nl = nl + n

  fmt = '(2x, a, 3p2f11.4)'
  nl=nl+1; write(lines(nl), *) ' '
  nl=nl+1; write(lines(nl), *)   'Orbit: [mm, mrad]'
  nl=nl+1; write(lines(nl), fmt) "X  X':", orb%vec(1), orb%vec(2)
  nl=nl+1; write(lines(nl), fmt) "Y  Y':", orb%vec(3), orb%vec(4)
  nl=nl+1; write(lines(nl), fmt) "Z  Z':", orb%vec(5), orb%vec(6)

  result_id = show_what

!----------------------------------------------------------------------
! universe
    
case ('universe')

  if (word1 == ' ') then
    ix_u = s%com%default_universe
  else
    read (word1, *, iostat = ios) ix_u
    if (ios /= 0) then
      nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
      return
    endif
    if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
      nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
      return
    endif
  endif

  u => s%u(ix_u)
  ix_branch = s%com%default_branch
  uni_branch => u%uni_branch(ix_branch)
  branch => lat%branch(ix_branch)
  lat_branch => u%model%lat_branch(ix_branch)

  nl = 0
  nl=nl+1; write(lines(nl), imt) 'Universe: ', ix_u
  nl=nl+1; write(lines(nl), imt) 'Branch:   ', ix_branch
  nl=nl+1; write(lines(nl), imt) '%n_d2_data_used        = ', u%n_d2_data_used
  nl=nl+1; write(lines(nl), imt) '%n_data_used           = ', u%n_data_used
  nl=nl+1; write(lines(nl), lmt) '%do_rad_int_calc       = ', u%calc%rad_int_for_data .or. u%calc%rad_int_for_plotting
  nl=nl+1; write(lines(nl), lmt) '%do_chrom_calc         = ', u%calc%chrom_for_data .or. u%calc%chrom_for_plotting
  nl=nl+1; write(lines(nl), lmt) '%calc%mat6             = ', u%calc%mat6
  nl=nl+1; write(lines(nl), lmt) '%calc%dynamic_aperture = ', u%calc%dynamic_aperture
  nl=nl+1; write(lines(nl), lmt) '%calc%one_turn_map     = ', u%calc%one_turn_map
  nl=nl+1; write(lines(nl), lmt) '%calc%track            = ', u%calc%track
  nl=nl+1; write(lines(nl), lmt) '%is_on                 = ', u%is_on
  nl=nl+1; write(lines(nl), amt) '%beam0_file            = ', trim(u%beam%beam0_file)
  nl=nl+1; write(lines(nl), amt) '%beam_all_file         = ', trim(u%beam%beam_all_file)
  nl=nl+1; write(lines(nl), amt) '%beam%saved_at:        = ', trim(u%beam%saved_at)
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), amt) 'Lattice name:           ', lat%lattice
  nl=nl+1; write(lines(nl), amt) 'Input_file_name:        ', lat%input_file_name
  nl=nl+1; write(lines(nl), amt) 'photon_type:            ', photon_type_name(lat%photon_type)
  nl=nl+1; write(lines(nl), lmt) 'Auto_Scale_Field_Phase: ', lat%auto_scale_field_phase
  nl=nl+1; write(lines(nl), lmt) 'Auto_scale_Field_Amp:   ', lat%auto_scale_field_amp
  nl=nl+1; write(lines(nl), lmt) 'Absolute_Time_Tracking: ', lat%absolute_time_tracking
  nl=nl+1; lines(nl) =           'Geometry:               ' // geometry_name(branch%param%geometry)
  nl=nl+1; write(lines(nl), lmt) 'global%rf_on:           ', s%global%rf_on
  nl=nl+1; write(lines(nl), imt) &
                'Elements used in tracking: From 1 through ', branch%n_ele_track
  if (branch%n_ele_max .gt. branch%n_ele_track) then
    nl=nl+1; write(lines(nl), '(a, i0, a, i0)') 'Lord elements:   ', &
                      branch%n_ele_track+1, '  through ', branch%n_ele_max
  else
    nl=nl+1; write(lines(nl), '(a)') 'There are NO Lord elements'
  endif

  nl=nl+1; write(lines(nl), '(a, f0.3)')   'Lattice length:             ', branch%param%total_length
  nl=nl+1; write(lines(nl), lmt)           'Aperture limits on?:        ', branch%param%aperture_limit_on

  if (branch%param%geometry == open$ .and. lat_branch%track_state /= moving_forward$) then
    if (s%global%track_type == 'beam') then
      nl=nl+1; write(lines(nl), '(a, i0)') 'Tracking: Lost beam at:     ', lat_branch%track_state
    else
      nl=nl+1; write(lines(nl), '(a, i0)') 'Tracking: Lost particle at: ', lat_branch%track_state
    endif
  endif

  if (.not. branch%param%stable) then
    nl=nl+1; write(lines(nl), '(a, l)') 'Model lattice stability: ', branch%param%stable
    nl=nl+1; write(lines(nl), '(a, l)') 'Design lattice stability:', u%design%lat%param%stable
    result_id = 'universe:unstable'
    return
  endif
 
  call radiation_integrals (u%model%lat, u%model%lat_branch(ix_branch)%orbit, u%model%modes, u%model%ix_rad_int_cache)
  call radiation_integrals (u%design%lat, u%design%lat_branch(ix_branch)%orbit, u%design%modes, u%design%ix_rad_int_cache)
  if (lat%param%geometry == closed$) then
    call chrom_calc (lat, s%global%delta_e_chrom, u%model%a%chrom, u%model%b%chrom)
    call chrom_calc (u%design%lat, s%global%delta_e_chrom, u%design%a%chrom, u%design%b%chrom)
  endif

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(17x, a)') '       X          |            Y'
  nl=nl+1; write(lines(nl), '(17x, a)') 'Model     Design  |     Model     Design'
  fmt  = '(1x, a10, 2es11.3, 2x, 2es11.3, 2x, a)'
  fmt2 = '(1x, a10, 2f11.3, 2x, 2f11.3, 2x, a)'
  fmt3 = '(1x, a10,        24x, 2es11.3, 2x, a)'
  f_phi = 1 / twopi
  l_lat = lat%param%total_length
  gamma2 = (lat%ele(0)%value(e_tot$) / mass_of(lat%param%particle))**2
  n = lat%n_ele_track
  if (lat%param%geometry == closed$) then
    nl=nl+1; write(lines(nl), fmt2) 'Q', f_phi*lat%ele(n)%a%phi, &
          f_phi*u%design%lat%ele(n)%a%phi, f_phi*lat%ele(n)%b%phi, &
          f_phi*u%design%lat%ele(n)%b%phi,  '! Tune'
    nl=nl+1; write(lines(nl), fmt2) 'Chrom', u%model%a%chrom, & 
          u%design%a%chrom, u%model%b%chrom, u%design%b%chrom, '! dQ/(dE/E)'
    nl=nl+1; write(lines(nl), fmt2) 'J_damp', u%model%modes%a%j_damp, &
        u%design%modes%a%j_damp, u%model%modes%b%j_damp, &
        u%design%modes%b%j_damp, '! Damping Partition #'
    nl=nl+1; write(lines(nl), fmt) 'Emittance', u%model%modes%a%emittance, &
        u%design%modes%a%emittance, u%model%modes%b%emittance, &
        u%design%modes%b%emittance, '! Meters'
  endif
  nl=nl+1; write(lines(nl), fmt) 'Alpha_damp', u%model%modes%a%alpha_damp, &
        u%design%modes%a%alpha_damp, u%model%modes%b%alpha_damp, &
        u%design%modes%b%alpha_damp, '! Damping per turn'
  nl=nl+1; write(lines(nl), fmt) 'I4', u%model%modes%a%synch_int(4), &
        u%design%modes%a%synch_int(4), u%model%modes%b%synch_int(4), &
        u%design%modes%b%synch_int(4), '! Radiation Integral'
  nl=nl+1; write(lines(nl), fmt) 'I5', u%model%modes%a%synch_int(5), &
        u%design%modes%a%synch_int(5), u%model%modes%b%synch_int(5), &
        u%design%modes%b%synch_int(5), '! Radiation Integral'
  nl=nl+1; write(lines(nl), fmt3) 'I6/gamma^2', u%model%modes%b%synch_int(6) / gamma2, &
        u%design%modes%b%synch_int(6) / gamma2, '! Radiation Integral'

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(19x, a)') 'Model     Design'
  fmt = '(1x, a12, 1p2e11.3, 3x, a)'
  if (lat%param%geometry == closed$) then
    call calc_z_tune(u%model%lat)
    nl=nl+1; write(lines(nl), '(1x, a12, 2f11.4, 3x, a)') 'Z_tune:', &
         -u%model%lat%z%tune/twopi, -u%design%lat%z%tune/twopi, '! The design value is calculated with RF on'
  endif
  nl=nl+1; write(lines(nl), fmt) 'Sig_E/E:', u%model%modes%sigE_E, &
            u%design%modes%sigE_E
  nl=nl+1; write(lines(nl), fmt) 'Energy Loss:', u%model%modes%e_loss, &
            u%design%modes%e_loss, '! Energy_Loss (eV / Turn)'
  nl=nl+1; write(lines(nl), fmt) 'J_damp:', u%model%modes%z%j_damp, &
        u%design%modes%z%j_damp, '! Longitudinal Damping Partition #'
  nl=nl+1; write(lines(nl), fmt) 'Alpha_damp:', u%model%modes%z%alpha_damp, &
        u%design%modes%z%alpha_damp, '! Longitudinal Damping per turn'
  nl=nl+1; write(lines(nl), fmt) 'Alpha_p:', u%model%modes%synch_int(1)/l_lat, &
               u%design%modes%synch_int(1)/l_lat, '! Momentum Compaction'
  nl=nl+1; write(lines(nl), fmt) 'I0:', u%model%modes%synch_int(0), &
               u%design%modes%synch_int(0), '! Radiation Integral'
  nl=nl+1; write(lines(nl), fmt) 'I1:', u%model%modes%synch_int(1), &
               u%design%modes%synch_int(1), '! Radiation Integral'
  nl=nl+1; write(lines(nl), fmt) 'I2:', u%model%modes%synch_int(2), &
               u%design%modes%synch_int(2), '! Radiation Integral'
  nl=nl+1; write(lines(nl), fmt) 'I3:', u%model%modes%synch_int(3), &
               u%design%modes%synch_int(3), '! Radiation Integral'

  result_id = show_what

!----------------------------------------------------------------------
! variable
    
case ('use')  

  nl=nl+1; lines(nl) = 'veto data *@*'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    do j = 1, s%u(i)%n_d2_data_used
      d2_ptr => s%u(i)%d2_data(j)
      call re_allocate (lines, nl+size(d2_ptr%d1)+10, .false.)
      do k = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
        d1_ptr => d2_ptr%d1(k)
        call location_encode(line, d1_ptr%d%useit_opt, d1_ptr%d%exists, lbound(d1_ptr%d, 1), err_flag = err)
        if (err .or. len_trim(line) + 50 > len(line)) then
          nl=nl+1; lines(nl) = 'veto data ' // trim(d1_ptr%name)
          do n = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
            if (.not. d1_ptr%d(n)%useit_opt) cycle
            if (nl + 100 > size(lines)) call re_allocate(lines, nl+100, .false.)
            nl=nl+1; write(lines(nl), '(3a, i0, a)') 'restore data ', trim(tao_d2_d1_name(d1_ptr)), '[', n, ']'
          enddo
        else
          if (line == '') cycle
          nl=nl+1; write(lines(nl), '(5a)') 'use data ', trim(tao_d2_d1_name(d1_ptr)), '[', trim(line), ']'
        endif
      enddo
    enddo
  enddo
  nl=nl+1; lines(nl) = ''

  call re_allocate (lines, nl+s%n_v1_var_used+10, .false.)
  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == ' ') cycle
    call re_allocate (lines, nl+200, .false.)
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1), err_flag = err)
    if (err .or. len_trim(line) + 50 > len(line)) then
      nl=nl+1; lines(nl) = 'veto var ' // trim(v1_ptr%name)
      do j = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
        if (.not. v1_ptr%v(j)%useit_opt) cycle
        if (nl + 100 > size(lines)) call re_allocate(lines, nl+100, .false.)
        nl=nl+1; write(lines(nl), '(3a, i0, a)') 'restore var ', trim(v1_ptr%name), '[', j, ']'
      enddo
    else
      nl=nl+1; write(lines(nl), '(5a)') 'use var ', trim(v1_ptr%name), '[', trim(line), ']'
    endif
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! variable

case ('value')

  call tao_evaluate_expression (stuff2, 0, .false., value, good, err)
  if (err) return

  if (size(value) == 1) then
    nl=nl+1; write(lines(nl), '(3x, es17.8)') value(1)
  else
    call re_allocate (lines, size(value)+100, .false.)
    do i = 1, size(value)
      nl=nl+1; write(lines(nl), '(i4, a, es17.8)') i, ':  ', value(i)
    enddo
  endif

  result_id = show_what

!----------------------------------------------------------------------
! variable

case ('variable')

  good_opt_only = .false.
  bmad_format = .false.
  print_header_lines = .true.

  do
    call tao_next_switch (stuff2, [ '-bmad_format     ', '-good_opt_only   ', & 
                                    '-no_label_lines  '], switch, err, ix_word)
    if (err) return

    select case (switch)  
    case ('') 
      exit
    case ('-bmad_format') 
      bmad_format = .true.
    case ('-good_opt_only') 
      good_opt_only = .true.
    case ('-no_label_lines')
      print_header_lines = .false.
    end select

  enddo

  word1 = stuff2(1:ix_word)
  
  if (.not. allocated (s%v1_var)) then
    nl=1; lines(1) = 'NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!'
    return 
  endif

  ! Bmad format

  if (bmad_format) then
    call tao_var_write (' ', good_opt_only)
    result_id = 'variable:bmad'
    return
  endif

  ! If 'n@' is present then write out stuff for universe n

  ix = index(word1, '@')
  if (ix /= 0) then
    if (ix == 1) then
      ix_u = s%com%default_universe
    else
      read (word1(:ix-1), *, iostat = ios) ix_u
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
        return
      endif
      if (ix_u == 0) ix_u = s%com%default_universe
      if (ix_u < 1 .or. ix_u > ubound(s%u, 1)) then
        nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
        return
      endif
    endif

    if (print_header_lines) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Variables controlling universe:', ix_u
      nl=nl+1; write(lines(nl), '(5x, a)') '                    '
      nl=nl+1; write(lines(nl), '(5x, a)') 'Name                '
    endif

    do i = 1, s%n_var_used
      if (.not. s%var(i)%exists) cycle
      found = .false.
      do j = 1, size(s%var(i)%this)
        if (s%var(i)%this(j)%ix_uni == ix_u) found = .true.
      enddo
      if (.not. found) cycle
      nl=nl+1; write(lines(nl), '(5x, a25, a40)') tao_var1_name(s%var(i)), &
                                                      tao_var_attrib_name(s%var(i))
    enddo
    result_id = 'variable:@'
    return
  endif

  ! If just "show var" then show all names

  if (word1 == ' ') then
    if (print_header_lines) then
      nl=nl+1; write(lines(nl), '(7x, a, t50, a)') 'Name', 'Using for Optimization'
    endif
    do i = 1, s%n_v1_var_used
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      call re_allocate (lines, nl+200, .false.)
      call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
      nl=nl+1; write(lines(nl), '(4x, 2a, i0, a, i0, a, t50, a)') &
                      trim(v1_ptr%name), '[', lbound(v1_ptr%v, 1), ':', &
                      ubound(v1_ptr%v, 1), ']', trim(line)
      
    enddo

    result_id = 'variable:'
    return
  endif

  ! are we looking at a range of locations?

  call tao_find_var(err, word1, v1_array, v_array) 
  if (err) return
  n_size = 0
  if (allocated(v_array)) n_size = size(v_array)

  ! v_ptr is valid then show the variable info.

  if (n_size == 1) then

    v_ptr => v_array(1)%v

    nl=nl+1; write(lines(nl), amt)  '%ele_name         = ', v_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%attrib_name      = ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  '%ix_attrib        = ', v_ptr%ix_attrib 
    nl=nl+1; write(lines(nl), imt)  '%ix_var           = ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  '%ix_dvar          = ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  '%ix_v1            = ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), rmt)  '%model            = ', v_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%base             = ', v_ptr%base_value

    if (.not. allocated (v_ptr%this)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      n = nl + 3*size(v_ptr%this) + 100
      if (size(lines) < n) call re_allocate(lines, n)
      do i = 1, size(v_ptr%this)
        nl=nl+1; write(lines(nl), '(4(a, i0))')  '%this(', i, ')%uni@branch>>ele:        ', &
                        v_ptr%this(i)%ix_uni, '@', v_ptr%this(i)%ix_branch, '>>', v_ptr%this(i)%ix_ele
        if (associated (v_ptr%this(i)%model_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_value: ', &
                                                            v_ptr%this(i)%model_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_value: <not associated>'
        endif
        if (associated (v_ptr%this(i)%base_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_value:  ', &
                                                            v_ptr%this(i)%base_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_value:  <not associated>'
        endif
      enddo
    endif

    if (associated (v_ptr%common%model_value)) then
      nl=nl+1; write(lines(nl), imt)  '%common%ix_uni:      ', v_ptr%common%ix_uni
      nl=nl+1; write(lines(nl), imt)  '%common%ix_ele:      ', v_ptr%common%ix_ele
      nl=nl+1; write(lines(nl), rmt)  '%common%Model_value: ', v_ptr%common%model_value
      nl=nl+1; write(lines(nl), rmt)  '%common%Base_value:  ', v_ptr%common%base_value
    endif

    nl=nl+1; write(lines(nl), rmt)  '%design           = ', v_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%old              = ', v_ptr%old_value
    nl=nl+1; write(lines(nl), rmt)  '%meas             = ', v_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%ref              = ', v_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%correction       = ', v_ptr%correction_value
    nl=nl+1; write(lines(nl), rmt)  '%high_lim         = ', v_ptr%high_lim
    nl=nl+1; write(lines(nl), rmt)  '%low_lim          = ', v_ptr%low_lim
    nl=nl+1; write(lines(nl), rmt)  '%step             = ', v_ptr%step
    nl=nl+1; write(lines(nl), rmt)  '%weight           = ', v_ptr%weight
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit      = ', v_ptr%delta_merit
    nl=nl+1; write(lines(nl), amt)  '%merit_type       = ', v_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit            = ', v_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%dmerit_dvar      = ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), rmt)  '%s                = ', v_ptr%s
    nl=nl+1; write(lines(nl), imt)  '%ix_key_table     = ', v_ptr%ix_key_table
    if (v_ptr%ix_key_table > 0 ) then
      nl=nl+1; write(lines(nl), rmt)  '%key_val0         = ', v_ptr%key_val0
      nl=nl+1; write(lines(nl), rmt)  '%key_delta        = ', v_ptr%key_delta
    endif
    nl=nl+1; write(lines(nl), lmt)  '%exists           = ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_var         = ', v_ptr%good_var
    nl=nl+1; write(lines(nl), lmt)  '%good_user        = ', v_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%good_opt         = ', v_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%useit_opt        = ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%useit_plot       = ', v_ptr%useit_plot

    result_id = 'variable:1:' // word1

  ! check if there is a variable number
  ! if no variable number requested, show a range

  elseif (size(v1_array) == 1) then

    nc = 0
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      nc = max(nc, len_trim(tao_var_attrib_name(v_ptr)))
    enddo

    if (print_header_lines) then
      line1 = '       Name'
      line1(nc+17:) = 'Meas         Model        Design  Useit_opt'
      nl=nl+1; write(lines(nl), '(2a)') 'Variable name:  ', v1_array(1)%v1%name
      nl=nl+1; lines(nl) = ' '
      nl=nl+1; lines(nl) = line1
    endif

    ! if a range is specified, show the variable range   
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      call re_allocate (lines, nl+200, .false.)
      nl=nl+1
      write(lines(nl), '(i6, 2x, a)') v_ptr%ix_v1, tao_var_attrib_name(v_ptr)
      write(lines(nl)(nc+9:), '(3es14.4, 7x, l)') v_ptr%meas_value, &
                 v_ptr%model_value, v_ptr%design_value, v_ptr%useit_opt
    enddo

    if (print_header_lines) then
      nl=nl+1; lines(nl) = line1
    endif

  else
    nl=1; lines(1) = '???'
    result_id = 'variable:?'
  endif

  result_id = show_what

!----------------------------------------------------------------------
! wall

case ('wall')

  by_s = .false.
  ix_sec = -1
  angle = 0
  sub_name = ''

  do
    call tao_next_switch (stuff2, ['-section  ', '-element  ', &
                                   '-angle    ', '-s        '], switch, err, ix_s2)
    if (err) return
    if (switch == '') exit
    select case (switch)

    case ('-angle')
      read (stuff2(1:ix_s2), *, iostat = ios) angle
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ANGLE.'
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-branch')
      branch => pointer_to_branch(stuff2(1:ix_s2), u%model%lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif
      ix_branch = branch%ix_branch
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-section')
      read (stuff2(1:ix_s2), *, iostat = ios) ix_sec
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ SECTION INDEX.'
        return
      endif
      call string_trim(stuff2(ix_s2+1:), stuff2, ix_s2)

    case ('-s')
      by_s = .true.
    end select

  enddo

  !-------

  if (.not. associated(branch%wall3d)) then
    nl=1; lines(nl) = 'No associated vacuum chamber wall.'
    result_id = 'wall:none'
    return
  endif

  wall => branch%wall3d

  if (ix_sec > 0) then 
    if (ix_sec > size(wall%section)) then
      nl=1; write(lines(nl), '(a, i0)') 'Section index larger than number of sections: ', size(wall%section)
      result_id = 'wall:sec:large'
      return
    endif

    wall_sec => wall%section(ix_sec)
    ele => pointer_to_ele(lat, wall_sec%ix_ele, wall_sec%ix_branch)
    nl=nl+1; write(lines(nl), '(5a)')            'ele:    ', trim(ele%name), '   (', trim(ele_loc_to_string(ele)), ')'
    nl=nl+1; write(lines(nl), '(2a)')            'type:   ', trim(wall3d_section_type_name(wall_sec%type))
    nl=nl+1; write(lines(nl), '(a, f14.6)')      'S:      ', wall_sec%s
    nl=nl+1; write(lines(nl), '(3(a, f10.6))')  '(x0, y0):  (', wall_sec%x0, ',', wall_sec%y0, ')'
    if (wall_sec%dr_ds == real_garbage$) then
      nl=nl+1; write(lines(nl), '(3(a, f10.6))')  'dr_ds:       Not-Set'
    else
      nl=nl+1; write(lines(nl), '(3(a, f10.6))')  'dr_ds:      ', wall_sec%dr_ds
    endif

    do j = 1, size(wall_sec%v)
      v => wall_sec%v(j)
      nl=nl+1; write(lines(nl), '(a, i0, a, 5f11.6)') &
                            'v(', j, ') =', v%x, v%y, v%radius_x, v%radius_y, v%tilt
    enddo

    return

  endif

  !-------

  if (by_s) then
    ix_s2 = index(stuff2, ':')
    if (ix_s2 == 0) then
      nl=1; lines(nl) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (stuff2(1:ix_s2-1), *, iostat = ios1) s1
    read (stuff2(ix_s2+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // stuff2
      return
    endif

    call bracket_index (wall%section%s, 1, size(wall%section), s1 - 1d-10, ix1)
    call bracket_index (wall%section%s, 1, size(wall%section), s2 + 1d-10, ix2)
    ix1 = ix1 + 1

  elseif (ix_s2 /= 0) then
    ix_s2 = index(stuff2, ':')
    if (ix_s2 == 0) then
      nl=1; lines(nl) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (stuff2(1:ix_s2-1), *, iostat = ios1) ix1
    read (stuff2(ix_s2+1:), *, iostat = ios2) ix2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // stuff2
      return
    endif

  else
    ix1 = 1; ix2 = min(200, size(wall%section))
  endif

  nl=nl+1; lines(nl) = '    Ix             S    ix_ele                 Ele      Type   Radius (mm)'

  do i = ix1, ix2
    wall_sec => wall%section(i)
    ele => pointer_to_ele (lat, wall_sec%ix_ele, wall_sec%ix_branch)
    
    call calc_wall_radius (wall%section(i)%v, cos(angle), sin(angle), r, z)
    nl=nl+1; write(lines(nl), '(i6, f14.6, a10, a20, a10, f14.3)') i, wall_sec%s, &
                trim(ele_loc_to_string(ele)), trim(ele%name), trim(wall3d_section_type_name(wall_sec%type)), 1000*r
  enddo

  nl=nl+1; lines(nl) = '    Ix             S    ix_ele                 Ele      Type   Radius (mm)'

  result_id = show_what

!----------------------------------------------------------------------
! wave

case ('wave')

  nl=nl+1; write(lines(nl), '(a, 2i4)') 'ix_a:', s%wave%ix_a1, s%wave%ix_a2
  nl=nl+1; write(lines(nl), '(a, 2i4)') 'ix_b:', s%wave%ix_b1, s%wave%ix_b2

  select case (s%wave%data_type)
  case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'beta.a', 'beta.b')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_a
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_b
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_Kick/Kick: ', s%wave%rms_rel_k
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi:       ', s%wave%rms_phi
    nl=nl+1; lines(nl) = ' '
    if (s%wave%data_type(1:4) == 'beta') then
      nl=nl+1; lines(nl) = 'Normalized Kick = kick * beta  [urad * meter]'
    else
      nl=nl+1; lines(nl) = 'Normalized Kick = kick * sqrt(beta)  [urad * sqrt(meter)]'
    endif
    nl=nl+1; lines(nl) = 'After Dat#    Norm_K       phi'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write(lines(nl), '(i9, f12.2, 1f10.3)') s%wave%kick(i)%ix_dat, &
                  1e6*s%wave%kick(i)%amp, s%wave%kick(i)%phi
    enddo

  case ('phase.a', 'phase.b')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_a
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_b
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_Kick/Kick: ', s%wave%rms_rel_k
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi:       ', s%wave%rms_phi
    nl=nl+1; write(lines(nl), '(a, f8.3, a)') &
                                    'Chi_C:           ', s%wave%chi_c, ' [Figure of Merit]'
    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Normalized Kick = k * l * beta [dimensionless]'
    nl=nl+1; lines(nl) = '   where k = quadrupole gradient [rad/m^2].'
    nl=nl+1; lines(nl) = 'After Dat#     Norm_K       phi'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write(lines(nl), '(i9, f12.4, f10.3)') s%wave%kick(i)%ix_dat, &
                  s%wave%kick(i)%amp, s%wave%kick(i)%phi
    enddo

  case ('cbar.11', 'cbar.12', 'cbar.22')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_+/Amp_+:  ', s%wave%rms_rel_as
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_-/Amp_-:  ', s%wave%rms_rel_ar
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_+/Amp_+:  ', s%wave%rms_rel_bs
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_-/Amp_-:  ', s%wave%rms_rel_br
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Kick |K+|   = ', 2*s%wave%amp_ba_s
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_K+/K+ = ', 2*s%wave%rms_rel_ks
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Kick |K-|   = ', 2*s%wave%amp_ba_r
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_K-/K- = ', 2*s%wave%rms_rel_kr
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi+:      ', s%wave%rms_phi_s
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi-:      ', s%wave%rms_phi_r
    nl=nl+1; write(lines(nl), '(a, f8.3, a)') &
                                    'Chi_a:           ', s%wave%chi_a, ' [Figure of Merit]'

    nl=nl+1; lines(nl) = 'After Dat#     Norm_K    phi+    phi-   phi_a   phi_b'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write(lines(nl), '(i11, f10.4, 4f8.3, 2f10.3)') &
            s%wave%kick(i)%ix_dat, &
            s%wave%kick(i)%amp, s%wave%kick(i)%phi_s, s%wave%kick(i)%phi_r, &
            (s%wave%kick(i)%phi_s+s%wave%kick(i)%phi_r)/2, &
            (s%wave%kick(i)%phi_s-s%wave%kick(i)%phi_r)/2
    enddo

  end select

  if (s%wave%n_kick > 10) then
    nl=nl+1; lines(nl) = ' etc...'
  endif


  result_id = show_what

!----------------------------------------------------------------------

case default

  nl=1; lines(1) = "INTERNAL ERROR, SHOULDN'T BE HERE!"

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine show_ele_data (u, ele, lines, nl)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum
type (ele_struct) ele

character(*) :: lines(:)
character(100) l1
integer nl, i

logical found_one

!

nl=nl+1; write(lines(nl), '(a)') "  "
write (l1, '(a, 20x, a)') "Data Name", &
          "Data Type             |  Model Value  |  Design Value |  Base Value"
nl=nl+1; lines(nl) = l1

found_one = .false.
do i = 1, size(u%data)
  if (u%data(i)%ix_ele == ele%ix_ele .and. u%data(i)%ix_branch == ele%ix_branch) then
    found_one = .true.
    datum => u%data(i)
    nl=nl+1; write(lines(nl), "(a, t30, a20, 3(1x, es15.5))") &
                trim(tao_datum_name(datum)),  datum%data_type, datum%model_value, &
                datum%design_value, datum%base_value 
  endif
enddo

if (found_one) then
  nl=nl+1; lines(nl) = l1
else
  write(lines(nl), '(a)') "No data associated with this element."
endif

end subroutine show_ele_data

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! contains

subroutine show_opt ()

implicit none

nl=nl+1; lines(nl) = 'Global optimization parameters:'
nl=nl+1; write(lines(nl), rmt) '  %de_lm_step_ratio              = ', s%global%de_lm_step_ratio
nl=nl+1; write(lines(nl), rmt) '  %de_var_to_population_factor   = ', s%global%de_var_to_population_factor
nl=nl+1; write(lines(nl), rmt) '  %lm_opt_deriv_reinit           = ', s%global%lm_opt_deriv_reinit
nl=nl+1; write(lines(nl), rmt) '  %lmdif_eps                     = ', s%global%lmdif_eps
nl=nl+1; write(lines(nl), rmt) '  %merit_stop_value              = ', s%global%merit_stop_value
nl=nl+1; write(lines(nl), rmt) '  %svd_cutoff                    = ', s%global%svd_cutoff
nl=nl+1; write(lines(nl), imt) '  %n_top10                       = ', s%global%n_top10
nl=nl+1; write(lines(nl), imt) '  %n_opti_loops                  = ', s%global%n_opti_loops
nl=nl+1; write(lines(nl), imt) '  %n_opti_cycles                 = ', s%global%n_opti_cycles
nl=nl+1; write(lines(nl), lmt) '  %derivative_recalc             = ', s%global%derivative_recalc
nl=nl+1; write(lines(nl), lmt) '  %svd_retreat_on_merit_increase = ', s%global%svd_retreat_on_merit_increase 
nl=nl+1; write(lines(nl), lmt) '  %derivative_uses_design        = ', s%global%derivative_uses_design
nl=nl+1; write(lines(nl), lmt) '  %opt_with_ref                  = ', s%global%opt_with_ref 
nl=nl+1; write(lines(nl), lmt) '  %opt_with_base                 = ', s%global%opt_with_base
nl=nl+1; write(lines(nl), amt) '  %optimizer                     = ', s%global%optimizer
nl=nl+1; lines(nl) = ''
nl=nl+1; lines(nl) = 'opti_de_param Parameters:'
nl=nl+1; write(lines(nl), rmt) '  %CR                   = ', opti_de_param%CR
nl=nl+1; write(lines(nl), rmt) '  %F                    = ', opti_de_param%F
nl=nl+1; write(lines(nl), rmt) '  %l_best               = ', opti_de_param%l_best
nl=nl+1; write(lines(nl), lmt) '  %binomial_cross       = ', opti_de_param%binomial_cross
nl=nl+1; write(lines(nl), lmt) '  %use_2nd_diff         = ', opti_de_param%use_2nd_diff
nl=nl+1; write(lines(nl), lmt) '  %randomize_F          = ', opti_de_param%randomize_F
nl=nl+1; write(lines(nl), lmt) '  %minimize_merit       = ', opti_de_param%minimize_merit

end subroutine show_opt

end subroutine tao_show_this

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine write_real (line, fmt, value)

implicit none

real(rp) value

integer id, ip, ix, wid, pl, wid_want, pl_want

character(16) fmt2, num_str
character(*) line, fmt

!

write (line, fmt) value

if (value == 0) return
if (line(1:1) == ' ') return
if (fmt(2:2) /= 'f' .and. fmt(2:2) /= 'F') return

! Always have a blank as first character.

id = index(fmt, '.')
ip = index(fmt, ')')

read (fmt(3:id-1), *) wid
read (fmt(id+1:ip-1), *) pl

wid_want = log10((1 + 1/10.0**(wid-2)) * abs(value)) + 1
if (value < 0) wid_want = wid_want + 1

pl_want = wid - wid_want - 2

if (pl_want == 0) then
  write (fmt2, '(a, i0, a, i0, a)') '(1x, f', wid, '.', pl_want, ')'
  write (line, fmt2) value
  line(wid+1:wid+1) = ' ' ! Get rid of '.'  
  return
endif

if (pl_want > 0) then  ! Can use F format.
  write (fmt2, '(a, i0, a, i0, a)') '(1x, f', wid-1, '.', pl_want, ')'
  write (line, fmt2) value
  return
endif

! Number is too large so switch to ES format.

wid_want = 4
if (value < 0) wid_want = wid_want + 1
if (abs(value) > 0.99d10) wid_want = wid_want + 1
if (abs(value) > 0.99d100) wid_want = wid_want + 1

if (wid < wid_want) return  ! Cannot do anything

pl = wid - wid_want - 1
if (pl < 0) pl = 0
write (fmt2, '(a, i0, a, i0, a)') '(es16.', pl, ')'
write (num_str, fmt2) value

id = index(num_str, '+')  ! Get rid of "+" sign in "E+nnn"
if (id /= 0) num_str = num_str(1:id-1) // num_str(id+1:)

id = index(num_str, 'E0')  ! Get rid of "0" in "E005"
if (id /= 0) num_str = num_str(1:id) // num_str(id+2:)

id = index(num_str, 'E0')  ! Get rid of "0" in "E005"
if (id /= 0) num_str = num_str(1:id) // num_str(id+2:)

call string_trim (num_str, num_str, ix)
line = ''
line(wid-ix+1:) = num_str

end subroutine

end module
