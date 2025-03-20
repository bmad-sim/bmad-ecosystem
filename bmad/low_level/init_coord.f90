!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

use bmad_routine_interface
implicit none

type (coord_struct) orb, orb_temp
type (ele_struct), optional :: ele
real(rp) :: vec(6)
real(rp), optional :: t_offset, E_photon, spin(3), s_pos
integer, optional :: element_end, particle, direction
logical, optional :: shift_vec6, random_on

!

orb_temp = coord_struct()
orb_temp%vec = vec

call init_coord2 (orb, orb_temp, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

end subroutine init_coord1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord2 (orb_out, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord2 (orb_out, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

use bmad_routine_interface
implicit none

type (coord_struct) orb_out, orb_in, orb
type (ele_struct), optional, target :: ele
type (branch_struct), pointer :: branch

real(rp), optional :: E_photon, t_offset, spin(3), s_pos
real(rp) p0c, e_tot, ref_time

integer, optional :: element_end, particle, direction
integer species, dir

logical, optional :: shift_vec6, random_on

character(*), parameter :: r_name = 'init_coord2'

!

orb = orb_in
if (present(ele)) then
  branch => pointer_to_branch(ele)
else
  branch => null()
endif

species = orb_in%species
if (present(particle)) then
  if (particle /= not_set$) species = particle
endif

if (present(particle)) orb%species = particle

if (orb%species == not_set$ .and. present(ele) .and. associated (branch)) then
  orb%species = default_tracking_species(branch%param)
endif

if (orb%species == not_set$) then
  call out_io (s_warn$, r_name, 'NO PARTICLE SPECIES GIVEN. USING POSITRONS AS DEFAULT!')
  orb%species = positron$
endif

orb%state = alive$

dir = integer_option(0, direction)
if (dir == 0) then
  if (orb%species == photon$ .and. orb%vec(6) > 0) orb%direction = 1
  if (orb%species == photon$ .and. orb%vec(6) < 0) orb%direction = -1
else
  orb%direction = dir
  if (orb%species == photon$) orb%vec(6) = orb%direction * abs(orb%vec(6))
endif

! Set location and species

if (present(element_end)) orb%location = element_end

if (orb%location == start_end$) then
  if (orb%direction == 1) then
    orb%location = upstream_end$
  else
    orb%location = downstream_end$
  endif
endif

! spin

if (present(spin)) orb%spin = spin

! Energy values

if (present(ele)) then
  if (.not. present(element_end)) then
    call out_io (s_fatal$, r_name, 'RULE: "ELEMENT_END" ARGUMENT MUST BE PRESENT IF "ELE" ARGUMENT IS.')
    call err_exit
  endif

  if (orb%location == downstream_end$ .or. ele%key == beginning_ele$) then
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    ref_time = ele%ref_time
    orb%s = ele%s
  elseif (orb%location == upstream_end$) then
    p0c = ele%value(p0c_start$)
    e_tot = ele%value(e_tot_start$)
    ref_time = ele%value(ref_time_start$)
    orb%s = ele%s_start
  else
    p0c = ele%value(p0c$)
    e_tot = ele%value(e_tot$)
    orb%s = real_option(ele%s, s_pos)
    ref_time = (ele%value(ref_time_start$) * (ele%s - orb%s) + ele%ref_time * (orb%s - ele%s_start)) / ele%value(l$)
  endif
endif

! Photon

if (orb%species == photon$) then
  if (present(ele)) then
    orb%p0c = p0c
    if (ele%key == photon_init$) then
      call init_photon_from_a_photon_init_ele (ele, branch%param, orb, random_on)
    endif
  endif

  orb%dt_ref = 0
  orb%beta = 1

  if (present(E_photon)) then
    if (E_photon /= 0) orb%p0c = E_photon
  endif

  ! If original particle is not a photon, photon is launched in same direciton as particle 
  if (orb_in%species /= photon$ .and. orb_in%species /= not_set$) then
    orb%vec(2:4:2) = orb%vec(2:4:2) / (1 + orb%vec(6))
    if (dir == 0) orb%direction = orb_in%direction
  endif
  orb%vec(6) = orb%direction * sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  if (orb%location == downstream_end$) then
    orb%vec(5) = ele%value(l$)
  else
    orb%vec(5) = 0
  endif
endif

! If ele is present...

orb%ix_ele = -1
orb%ix_branch = -1

if (present(ele)) then

  if (ele%slave_status == slice_slave$) then
    orb%ix_ele = ele%lord%ix_ele
    orb%ix_branch = ele%lord%ix_branch
  else
    orb%ix_ele = ele%ix_ele
    orb%ix_branch = ele%ix_branch
  endif

  if (ele%key == beginning_ele$) orb%location = downstream_end$

  if (orb%species /= photon$) then

    orb%p0c = p0c

    ! E_gun shift. Only time p0c_start /= p0c for an init_ele is when there is an e_gun present in the branch.
    if (logic_option(.true., shift_vec6)) then
      if (ele%key == beginning_ele$) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_start$) - ele%value(p0c$)) / ele%value(p0c$)
      elseif ((ele%key == e_gun$ .and. orb%location == upstream_end$) .or. (ele%key == marker$ .and. ele%value(e_tot_ref_init$) /= 0)) then
        orb%vec(6) = orb%vec(6) + (ele%value(p0c_ref_init$) - ele%value(p0c$)) / ele%value(p0c$)
      endif
    endif

    if (orb%vec(6) == 0) then
      orb%beta = p0c / e_tot
    else
      call convert_pc_to (p0c * (1 + orb%vec(6)), orb%species, beta = orb%beta)
    endif

    ! Do not set %t if %beta = 0 since %t may be a good value.

    if (orb%vec(5) == real_garbage$) then
      orb%vec(5) = orb%beta * c_light * (ref_time - orb%t)

    elseif (orb%beta == 0) then
      if (orb%vec(5) /= 0) then
        call out_io (s_error$, r_name, 'Z-POSITION IS NONZERO WITH BETA = 0.', &
                                       'THIS IS NONSENSE SO SETTING Z TO ZERO.')
        orb%vec(5) = 0
      endif
    else
      orb%t = ref_time - orb%vec(5) / (orb%beta * c_light)
      if (present(t_offset)) orb%t = orb%t + t_offset
      if (orb%ix_turn /= 0 .and. associated(ele%branch)) then
        orb%t = orb%t + orb%ix_turn * (ele%branch%ele(ele%branch%n_ele_track)%ref_time - ele%branch%ele(0)%ref_time)
      endif
    endif
  endif

endif

orb_out = orb

end subroutine init_coord2

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)
! 
! Subroutine to initialize a coord_struct. 
! This subroutine is overloaded by init_coord. See init_coord for more details.
!-

subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)

use bmad_routine_interface
implicit none

type (coord_struct) orb
type (ele_struct), optional :: ele
real(rp), optional :: t_offset, E_photon, spin(3)
integer, optional :: element_end, particle, direction
logical, optional :: shift_vec6

!

call init_coord2 (orb, coord_struct(), ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin)

end subroutine init_coord3
