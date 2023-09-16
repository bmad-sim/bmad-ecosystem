!+
! Subroutine init_fringe_info (fringe_info, ele, orbit, leng_sign)
!
! Routine to initalize a fringe_info_struct for a particular lattice element.
!
! Note: An element that has a longitudinal field component or an electric field is considered to 
! have a fringe at both ends even though the slice ends may not actually abut the ends of the full element.
! This is done to make computed transfer matrices symplectic. Non-symplectic matrices
! will mess up Twiss calculations. EG: Think of the case of a solenoid with a superimposed marker.
!
! Input:
!   ele           -- ele_struct: Lattice element associated with fringe_info.
!   orbit         -- coord_struct, optional: Particle position. Must be present for a full init.
!                     If not full init only fringe_info%has_fringe will be set.
!   leng_sign     -- integer, optional: Is element length positive (+1) or negative (-1)? 
!                     Must be present if orbit is present.
!
! Output:
!   fringe_info   -- fringe_info_struct: Fringe information.
!-

subroutine init_fringe_info (fringe_info, ele, orbit, leng_sign)

use bmad_interface, dummy => init_fringe_info

implicit none

type (fringe_field_info_struct) fringe_info
type (ele_struct) ele
type (coord_struct), optional :: orbit
type (ele_struct), pointer :: lord

real(rp) s_orb

integer, optional :: leng_sign
integer i

character(*), parameter :: r_name = 'init_finge_info'

! Calc fringe_info%has_fringe

select case (ele%key)
case (em_field$, solenoid$, sol_quad$, wiggler$, rfcavity$, lcavity$, crab_cavity$, &
                                                                 elseparator$, e_gun$, sad_mult$, custom$)
  fringe_info%has_fringe = .true.
case default
  if (attribute_name(ele, fringe_at$) /= 'FRINGE_AT') then
    fringe_info%has_fringe = .false.
  elseif ((nint(ele%value(fringe_type$)) == none$ .or. nint(ele%value(fringe_at$)) == no_end$) .and. &
                                                                        .not. associated(ele%a_pole_elec)) then
    fringe_info%has_fringe = .false.
  else
    fringe_info%has_fringe = .true.
  endif
end select

nullify (fringe_info%hard_ele)
fringe_info%particle_at = first_track_edge$

! Full init

if (.not. fringe_info%has_fringe .or. .not. present(orbit)) return

s_orb = orbit%s - ele%s_start

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  call re_allocate(fringe_info%location, ele%n_lord)
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%key == overlay$ .or. lord%key == group$) cycle
    call init_this_ele (lord, i, leng_sign)
  enddo
  if (all(fringe_info%location == nowhere$)) fringe_info%has_fringe = .false.

else
  call re_allocate(fringe_info%location, 1)
  call init_this_ele (ele, 1, leng_sign)
endif

!-------------------------------------------------------------------------
contains 

subroutine init_this_ele (this_ele, ix_loc, leng_sign)

type (ele_struct) this_ele
real(rp) s_off, s1, s2, s_hard_entrance, s_hard_exit, ds_small, leng, sa, sb
integer ix_loc, leng_sign

!

s_off = this_ele%s_start - ele%s_start

leng = this_ele%value(l$)
select case (this_ele%key)
case (rfcavity$, lcavity$)
  s1 = s_off + (leng - this_ele%value(l_active$)) / 2  ! Distance from entrance end to active edge
  s2 = s_off + (leng + this_ele%value(l_active$)) / 2  ! Distance from entrance end to the other active edge
case default
  s1 = s_off         ! Distance from entrance end to hard edge
  s2 = s_off + leng  ! Distance from entrance end to the other hard edge
end select

if (leng_sign > 0) then
  s_hard_entrance = min(s1, s2)
  s_hard_exit     = max(s1, s2)
else
  s_hard_entrance = max(s1, s2)
  s_hard_exit     = min(s1, s2)
endif

ds_small = bmad_com%significant_length

! If there are no longitudinal fields nor electric fields then mark the element as having
! no fringe if the there is a lord with edges outside of the slice. 

if (ele%key /= rfcavity$ .and. ele%key /= lcavity$ .and. ele%key /= crab_cavity$ .and. &
      (ele%field_calc == bmad_standard$ .or. ele%tracking_method == bmad_standard$) .and. &
      .not. associated(ele%a_pole_elec)) then
  sa = min(0.0_rp, ele%value(l$)) - bmad_com%significant_length
  sb = max(0.0_rp, ele%value(l$)) + bmad_com%significant_length

  if ((s1 < sa .or. s1 > sb) .and. (s2 < sa .or. s2 > sb)) then
    fringe_info%location(ix_loc) = nowhere$
    return
  endif
endif

!

if (orbit%direction*orbit%time_dir * ele%orientation == 1) then
  if ((orbit%location == upstream_end$ .and. ele%orientation == 1) .or. &
                                    (orbit%location == downstream_end$ .and. ele%orientation == -1)) then
    fringe_info%location(ix_loc) = entrance_end$

  elseif (orbit%location == inside$) then
    if (leng_sign * s_hard_entrance > leng_sign * (s_orb + ds_small)) then
      fringe_info%location(ix_loc) = entrance_end$
    elseif (leng_sign * s_hard_exit > leng_sign * s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = exit_end$
    endif

  else
    call out_io (s_fatal$, r_name, 'CONFUSED FORWARD DIRECTION INITIALIZATION!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

else  ! orbit%direction*orbit%time_dir * ele%orientation = -1
  if ((orbit%location == downstream_end$ .and. ele%orientation == 1) .or. &
                                      (orbit%location == upstream_end$ .and. ele%orientation == -1) ) then
      fringe_info%location(ix_loc) = exit_end$

  elseif (orbit%location == inside$) then
    if (leng_sign * s_hard_exit < leng_sign * (s_orb - ds_small)) then
      fringe_info%location(ix_loc) = exit_end$
    elseif (leng_sign * s_hard_entrance < leng_sign * s_orb) then
      fringe_info%location(ix_loc) = inside$
    else
      fringe_info%location(ix_loc) = entrance_end$
    endif

  else
    call out_io (s_fatal$, r_name, 'CONFUSED REVERSE DIRECTION INITIALIZATION!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

end subroutine init_this_ele

end subroutine init_fringe_info
