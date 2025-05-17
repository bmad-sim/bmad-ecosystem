!+
! Subroutine multipole_ele_to_kt (ele, use_ele_tilt, ix_pole_max, knl, tilt, pole_type, include_kicks)
!
! Subroutine to put the multipole components (strength and tilt)
! into 2 vectors along with the appropriate scaling for the relative tracking charge, etc.
!
! Note: tilt(:) does includes ele%value(tilt_tot$).
!
! Note: To save time, if the element does not have ele%a/b_pole allocated, the knl array 
!       will NOT be set to zero (ix_pole_max will be set to -1). 
!       That is, the ix_pole_max argument needs to be tested before any calculations!
!
! Input:
!   ele          -- Ele_struct: Lattice element.
!   use_ele_tilt -- Logical: If True then include ele%value(tilt_tot$) in calculations. 
!                     use_ele_tilt is ignored in the case of multipole$ elements.
!   pole_type    -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   include_kicks -- integer, optional: Possibilities are: 
!            no$                      -- Default. Do not include any kick components in a and b multipoles. 
!            include_kicks$           -- Include hkick/vkick/dg in the n = 0 components.
!                                           Also included are quad k1, sextupole k2 and octupole k3 components.
!
! Output:
!   ix_pole_max         -- Integer: Index of largest nonzero pole.
!   knl(0:n_pole_maxx)  -- Real(rp): Vector of strengths, MAD units.
!   tilt(0:n_pole_maxx) -- Real(rp): Vector of tilts.
!-

subroutine multipole_ele_to_kt (ele, use_ele_tilt, ix_pole_max, knl, tilt, pole_type, include_kicks)

use multipole_mod, dummy => multipole_ele_to_kt

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord

real(rp) knl(0:), tilt(0:), a(0:n_pole_maxx), b(0:n_pole_maxx), factor
real(rp), pointer :: knl_pole(:), tilt_pole(:), ksl_pole(:)

integer ix_pole_max
integer, optional :: pole_type, include_kicks
integer n, ix_max

logical use_ele_tilt

! Multipole type element case

if (ele%key == multipole$) then
  if (integer_option(magnetic$, pole_type) == electric$) then
    knl = 0
    tilt = 0
    ix_pole_max = 0
    return
  endif

  if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
    lord => pointer_to_super_lord(ele)
    call pointer_to_ele_multipole (lord, knl_pole, tilt_pole, ksl_pole, pole_type)
  else
    call pointer_to_ele_multipole (ele, knl_pole, tilt_pole, ksl_pole, pole_type)
  endif

  if (all(ksl_pole == 0)) then
    knl  = knl_pole
    tilt = tilt_pole + ele%value(tilt_tot$)
  else
    do n = 0, n_pole_maxx
      knl(n)  = sqrt(knl_pole(n)**2 + ksl_pole(n)**2)
      tilt(n) = tilt_pole(n) - atan2(ksl_pole(n), knl_pole(n)) / (n + 1)
      ! In case the user looks at this, make tilt(n) to be in the range [-pi, pi]/(n+1)
      if (2 * (n + 1) * abs(tilt(n)) > pi) then
        knl(n) = -knl(n)
        tilt(n) = sign_of(tilt(n)) * (abs(tilt(n)) - pi/(n+1))
      endif
    enddo
  endif

  if (ele%field_master .and. ele%value(p0c$) /= 0) then
    factor = charge_of(ele%ref_species) * c_light / ele%value(p0c$)
    knl = factor * knl
  endif

  ix_pole_max = max_nonzero(0, knl)
  return
endif

! Everything else

call multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks)

if (ix_pole_max > -1) then
  call multipole_ab_to_kt (a, b, knl, tilt)
else
  knl = 0
  tilt = 0
endif

end subroutine multipole_ele_to_kt

