!+
! Subroutine offset_particle (ele, param, coord, set, 
!                        set_canonical, set_tilt, set_multipoles, set_hvkicks)
! Subroutine to effectively offset an element by instead offsetting
! the particle position to correspond to the local element coordinates.
! Options:
!   Conversion between angle (x', y') and canonical (P_x, P_y) momenta.
!   Using the element tilt in the offset.
!   Using the HV kicks.
!   Using the multipoles.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!     %value(x_offset$) -- Horizontal offset of element.
!     %value(x_pitch$)  -- Horizontal roll of element.
!     %value(tilt$)     -- titlt of element.
!   coord
!     %z%vel            -- Energy deviation dE/E. 
!                          Used to modify %x%vel and %y%vel
!   param     -- Param_struct:
!     %particle   -- What kind of particle (for elseparator elements).
!   set       -- Logical: 
!                   T -> Translate from lab coords to the local element coords.
!                   F -> Translate back to lab coords.
!   set_canonical  -- Logical, optional: Default is True.
!                   T -> Convert between (P_x, P_y) and (x', y') also.
!                   F -> No conversion between (P_x, P_y) and (x', y').
!   set_tilt       -- Logical, optional: Default is True.
!                   T -> Rotate using ele%value(tilt$) and 
!                            ele%value(roll$) for sbends.
!                   F -> Do not rotate
!   set_multipoles -- Logical, optional: Default is True.
!                   T -> 1/2 of the multipole is applied.
!   set_hvkicks    -- Logical, optional: Default is True.
!                                               
! Output:
!     coord -- Coord_struct: Coordinates of particle.
!-

!$Id$
!$Log$
!Revision 1.2  2002/02/23 20:32:21  dcs
!Double/Single Real toggle added
!
!Revision 1.1  2002/01/08 21:44:42  dcs
!Aligned with VMS version  -- DCS
!

#include "CESR_platform.inc"
                                                              
subroutine offset_particle (ele, param, coord, set, &
              set_canonical, set_tilt, set_multipoles, set_hvkicks)

  use bmad

  implicit none

  type (ele_struct) ele
  type (coord_struct) coord
  type (param_struct) param

  real(rdef) E_rel, knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rdef) del_x_vel, del_y_vel

  integer n

  logical set
  logical, optional :: set_canonical, set_tilt, set_multipoles, set_hvkicks
  logical set_canon, set_multi, set_hv, set_t

!---------------------------------------------------------------         
! E_rel               

  E_rel = (1 + coord%z%vel)

! init

  if (present(set_canonical)) then
    set_canon = set_canonical
  else
    set_canon = .true.
  endif

  if (present(set_multipoles)) then
    set_multi = set_multipoles
  else
    set_multi = .true.
  endif

  if (present(set_hvkicks)) then
    set_hv = set_hvkicks
  else
    set_hv = .true.
  endif

  if (present(set_tilt)) then
    set_t = set_tilt
  else
    set_t = .true.
  endif

!----------------------------------------------------------------
! Set...

  if (set) then

! Set s_offset

    if (ele%value(s_offset$) /= 0) then
      coord%x%pos = coord%x%pos + ele%value(s_offset$) * coord%x%vel / E_rel
      coord%y%pos = coord%y%pos + ele%value(s_offset$) * coord%y%vel / E_rel
    endif

! Set: Offset and pitch

    if (ele%value(x_offset$) /= 0 .or. ele%value(y_offset$) /= 0 .or. &
              ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) then
      if (ele%value(l$) == 0) then
        coord%x%pos = coord%x%pos - ele%value(x_offset$)
        coord%y%pos = coord%y%pos - ele%value(y_offset$)
      else
        coord%x%pos = coord%x%pos - ele%value(x_offset$) +  &
                           ele%value(x_pitch$) * ele%value(l$) / 2
        coord%x%vel = coord%x%vel - ele%value(x_pitch$) * E_rel
        coord%y%pos = coord%y%pos - ele%value(y_offset$) +  &
                           ele%value(y_pitch$) * ele%value(l$) / 2
        coord%y%vel = coord%y%vel - ele%value(y_pitch$) * E_rel
      endif
    endif

! Set: HV kicks.
! HV kicks must come after s_offset but before any tilts are applied.
! Note: Change in %vel is NOT dependent upon energy since we are using
! canonical momentum.

    if (set_hv) then
      if (ele%is_on) then
        if (ele%key == elseparator$ .and. param%particle < 0) then
          coord%x%vel = coord%x%vel - ele%value(hkick$) / 2
          coord%y%vel = coord%y%vel - ele%value(vkick$) / 2
        else
          coord%x%vel = coord%x%vel + ele%value(hkick$) / 2
          coord%y%vel = coord%y%vel + ele%value(vkick$) / 2
        endif
      endif
    endif

! Set: Multipoles

    if (set_multi .and. ele%nonzero_multipoles) then
      call multipole_ele_to_kt(ele, param%particle, knl, tilt, .true.)
      knl = knl / 2
      do n = 0, n_pole_maxx
        call multipole_kick (knl(n), tilt(n), n, coord)
      enddo
    endif

! Set: Tilt
! A non-zero roll has a zeroth order effect that must be included  

    if (set_t) then

      if (ele%key == sbend$) then
        if (ele%value(roll$) /= 0) then
          if (abs(ele%value(roll$)) < 0.01) then
            del_x_vel = ele%value(angle$) * ele%value(roll$)**2 / 4
          else
            del_x_vel = ele%value(angle$) * (1 - cos(ele%value(roll$))) / 2
          endif
          del_y_vel = -ele%value(angle$) * sin(ele%value(roll$)) / 2
          coord%x%vel = coord%x%vel + del_x_vel
          coord%y%vel = coord%y%vel + del_y_vel
        endif
        call tilt_coords (ele%value(tilt$)+ele%value(roll$), coord%vec, set$)
      else
        call tilt_coords (ele%value(tilt$), coord%vec, set$)
      endif

    endif

! Set: P_x, P_y -> x', y' 

    if (set_canon .and. coord%z%vel /= 0) then
      coord%x%vel = coord%x%vel / E_rel
      coord%y%vel = coord%y%vel / E_rel
    endif

!----------------------------------------------------------------
! Unset... 

  else

! Unset: x', y' -> P_x, P_y 

    if (set_canon .and. coord%z%vel /= 0) then
      coord%x%vel = coord%x%vel * E_rel
      coord%y%vel = coord%y%vel * E_rel
    endif

! Unset: Tilt

    if (set_t) then

      if (ele%key == sbend$) then
        call tilt_coords (ele%value(tilt$)+ele%value(roll$), coord%vec, unset$) 
        if (ele%value(roll$) /= 0) then  
          coord%x%vel = coord%x%vel + del_x_vel
          coord%y%vel = coord%y%vel + del_y_vel
        endif
      else
        call tilt_coords (ele%value(tilt$), coord%vec, unset$)   
      endif

    endif

! Unset: Multipoles

    if (set_multi .and. ele%nonzero_multipoles) then
      do n = 0, n_pole_maxx
        call multipole_kick (knl(n), tilt(n), n, coord)
      enddo
    endif

! unset: HV kicks.
! HV kicks must come after s_offset but before any tilts are applied.
! Note: Change in %vel is NOT dependent upon energy since we are using
! canonical momentum.

    if (set_hv) then
      if (ele%is_on) then
        if (ele%key == elseparator$ .and. param%particle < 0) then
          coord%x%vel = coord%x%vel - ele%value(hkick$) / 2
          coord%y%vel = coord%y%vel - ele%value(vkick$) / 2
        else
          coord%x%vel = coord%x%vel + ele%value(hkick$) / 2
          coord%y%vel = coord%y%vel + ele%value(vkick$) / 2
        endif
      endif
    endif

! Unset: Offset and pitch

    if (ele%value(x_offset$) /= 0 .or. ele%value(y_offset$) /= 0 .or. &
              ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) then
      if (ele%key == multipole$ .or. ele%key == ab_multipole$) then
        coord%x%pos = coord%x%pos + ele%value(x_offset$)
        coord%y%pos = coord%y%pos + ele%value(y_offset$)
      else
        coord%x%pos = coord%x%pos + ele%value(x_offset$) + &
                           ele%value(x_pitch$) * ele%value(l$) / 2
        coord%x%vel = coord%x%vel + ele%value(x_pitch$) * E_rel
        coord%y%pos = coord%y%pos + ele%value(y_offset$) +  &
                           ele%value(y_pitch$) * ele%value(l$) / 2
        coord%y%vel = coord%y%vel + ele%value(y_pitch$) * E_rel
      endif
    endif

    if (ele%value(s_offset$) /= 0) then
      coord%x%pos = coord%x%pos - ele%value(s_offset$) * coord%x%vel / E_rel
      coord%y%pos = coord%y%pos - ele%value(s_offset$) * coord%y%vel / E_rel
    endif

  endif

end subroutine
                          

