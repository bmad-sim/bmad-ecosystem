!+
! Subroutine multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks, b1)
!                             
! Subroutine to extract the ab multipole values of an element.
!
! Note: The ab values will be scalled by the strength of the element if scale_multipole = T is set for the element.
!
! Input:
!   ele           -- ele_struct: Element.
!     %value()      -- ab_multipole values.
!   use_ele_tilt  -- logical: If True then include ele%value(tilt_tot$) in calculations.
!                      use_ele_tilt is ignored in the case of multipole$ elements.
!   pole_type     -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   include_kicks -- integer, optional: Ignored for for pole_type == electric$ for non-elseparator elements.
!                      Possibilities are: 
!            no$                      -- Default. Do not include any kick components in a and b multipoles. 
!            include_kicks$           -- Include hkick/vkick/dg in the n = 0 components.
!                                           Also included are quad k1, sextupole k2 and octupole k3 components.
!
! Output:
!   ix_pole_max       -- Integer: Index of largest nonzero pole. Set to -1 if all multipoles are zero.
!   a(0:n_pole_maxx)  -- real(rp): Array of scalled multipole values.
!   b(0:n_pole_maxx)  -- real(rp): Array of scalled multipole values.
!   b1                -- real(rp), optional: If present and True then b1 is set to the value of the b(1) component
!                         of the b(:) array and b(1) is set to zero. Also ix_pole_max is ajusted as needed.
!                         This is used by routines that want to handle b(1) in a special way in tracking.
!-

recursive subroutine multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks, b1)

use bmad_interface, dummy => multipole_ele_to_ab

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (multipole_cache_struct), pointer :: cache

real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp), optional :: b1
real(rp) const, radius, factor, k1_lord, this_a_kick(0:3), this_b_kick(0:3), a_kick(0:3), b_kick(0:3)
real(rp) this_a(0:n_pole_maxx), this_b(0:n_pole_maxx)
real(rp), pointer :: a_pole(:), b_pole(:), ksl_pole(:)

integer ix_pole_max
integer, optional :: pole_type, include_kicks
integer i, n, p_type, include_kck

logical use_ele_tilt, can_use_cache

character(*), parameter :: r_name = 'multipole_ele_to_ab'

! Init

ix_pole_max = -1

if (.not. ele%is_on) then
  a = 0;  b = 0
  if (present(b1)) b1 = 0
  return
endif

! Use cache if possible. 
! Caching requires intelligent bookkeeping to mark when the cache goes stale.

p_type = integer_option(magnetic$, pole_type)
include_kck = integer_option(no$, include_kicks)
can_use_cache = (.not. bmad_com%auto_bookkeeper .and. allocated(ele%multipole_cache))

! Note: slice_slave and super_slave elements have multipoles stored in their lords.

cache => ele%multipole_cache
if (can_use_cache) then
  if (p_type == magnetic$ .and. cache%ix_pole_mag_max /= invalid$) then
    ix_pole_max = cache%ix_pole_mag_max
    if (ix_pole_max == invalid$ .and. (include_kck == no$ .or. cache%ix_kick_mag_max == invalid$)) then
      a = 0;  b = 0
      if (present(b1)) b1 = 0
      return
    endif

    a = cache%a_pole_mag
    b = cache%b_pole_mag

    if (cache%ix_kick_mag_max > -1 .and. include_kck == include_kicks$) then
      a(0:3) = a(0:3) + cache%a_kick_mag
      b(0:3) = b(0:3) + cache%b_kick_mag
      ix_pole_max = max(ix_pole_max, cache%ix_kick_mag_max)
    endif

    if (use_ele_tilt) call tilt_this_multipole(ele, +1, a, b, ix_pole_max)
    if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)
    return

  elseif (p_type == electric$ .and. cache%ix_pole_elec_max /= invalid$) then
    a = cache%a_pole_elec
    b = cache%b_pole_elec
    ix_pole_max = cache%ix_pole_elec_max
    if (cache%ix_kick_elec_max > -1 .and. include_kck == include_kicks$) then
      a(0) = a(0) + cache%a_kick_elec(0)
      b(0) = b(0) + cache%b_kick_elec(0)
      ix_pole_max = max(ix_pole_max, cache%ix_kick_elec_max)
    endif

    if (use_ele_tilt) call tilt_this_multipole(ele, +1, a, b, ix_pole_max)
    if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)
    return
  endif
endif

! Multipole ele 

if (ele%key == multipole$) then
  if (p_type == electric$) then
    a = 0;  b = 0
    if (present(b1)) b1 = 0
    return
  endif

  if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
    lord => pointer_to_lord(ele, 1)
  else
    lord => ele
  endif

  call pointer_to_ele_multipole (ele, a_pole, b_pole, ksl_pole, pole_type)
  call multipole_kt_to_ab (a_pole, ksl_pole, b_pole, a, b)
  ix_pole_max = max_nonzero(0, a, b)

  if (ele%field_master .and. ele%value(p0c$) /= 0) then
    factor = charge_of(ele%ref_species) * c_light / ele%value(p0c$)
    a = factor * a
    b = factor * b
  endif

  if (can_use_cache) then
    a_kick = 0;  b_kick = 0
    call load_this_cache(cache, p_type, ix_pole_max, a, b, a_kick, b_kick)
  endif
  if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)
  return
endif

! Slice slaves and super slaves have their associated multipoles stored in the lord.

a = 0
b = 0

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  ! Easy case
  if (ele%n_lord == 1) then
    lord => pointer_to_lord(ele, 1)
    call multipole_ele_to_ab (lord, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks)
    if (p_type == magnetic$ .and. ix_pole_max > -1) then
      n = ix_pole_max
      a(0:n) = a(0:n) * (ele%value(l$) / lord%value(l$))
      b(0:n) = b(0:n) * (ele%value(l$) / lord%value(l$))
    endif
    if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)
    return
  endif
  ! Not so easy case
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%key == group$ .or. lord%key == overlay$ .or. lord%key == girder$) cycle
    if (.not. lord%is_on) cycle
    call multipole_ele_to_ab (lord, .true., n, this_a, this_b, pole_type, include_kicks)
    if (n > -1) then
      if (p_type == magnetic$) then
        a(0:n) = a(0:n) + this_a(0:n) * (ele%value(l$) / lord%value(l$))
        b(0:n) = b(0:n) + this_b(0:n) * (ele%value(l$) / lord%value(l$))
      else
        a(0:n) = a(0:n) + this_a(0:n)
        b(0:n) = b(0:n) + this_b(0:n)
      endif
      ix_pole_max = max(ix_pole_max, n)
    endif
  enddo

  if (.not. use_ele_tilt) call tilt_this_multipole(ele, -1, a, b, ix_pole_max)
  if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)
  return
endif

! Not a slave

call pointer_to_ele_multipole (ele, a_pole, b_pole, ksl_pole, pole_type)
if (ele%multipoles_on .and. associated(a_pole)) then
  call convert_this_ab (ele, p_type, a_pole, b_pole, a, b)
endif

call this_ele_non_multipoles (ele, can_use_cache, a_kick, b_kick)

ix_pole_max = max_nonzero(0, a, b)
if (can_use_cache) call load_this_cache(cache, p_type, ix_pole_max, a, b, a_kick, b_kick)

if (include_kck == include_kicks$) then
  a(0:3) = a(0:3) + a_kick
  b(0:3) = b(0:3) + b_kick
endif

ix_pole_max = max_nonzero(0, a, b)

if (use_ele_tilt) call tilt_this_multipole(ele, +1, a, b, ix_pole_max)
if (present(b1)) b1 = pull_this_b1(a, b, ix_pole_max)

!---------------------------------------------
contains

subroutine this_ele_non_multipoles (ele, can_use_cache, a_kick, b_kick)

type (ele_struct) ele

real(rp) a_kick(0:3), b_kick(0:3)
real(rp) an, bn, hk, vk, tilt, sin_t, cos_t

logical can_use_cache

! Express non-multipoles in element body frame.

a_kick = 0
b_kick = 0

if (.not. can_use_cache .and. include_kck == no$) return

! Magnetic elements
select case (p_type)
case (magnetic$)
  select case (ele%key)
  case (hkicker$)
    b_kick(0) = -ele%value(kick$)
  case (vkicker$)
    a_kick(0) = ele%value(kick$)
  case (ac_kicker$, kicker$)
    b_kick(0) = -ele%value(hkick$)
    a_kick(0) =  ele%value(vkick$)
  case (elseparator$)
    ! Kicks are electric
  case default
    ! tilt to element body coords
    hk = ele%value(hkick$)
    vk = ele%value(vkick$)
    if (has_hkick_attributes(ele%key) .and. (hk /= 0 .or. vk /= 0)) then
      if (ele%key == sbend$) then
        tilt = ele%value(ref_tilt_tot$)
      else
        tilt = ele%value(tilt_tot$) 
      endif

      if (tilt /= 0) then
        cos_t = cos(tilt)
        sin_t = sin(tilt)
        b_kick(0) = -hk * cos_t - vk * sin_t
        a_kick(0) = -hk * sin_t + vk * cos_t
      else
        b_kick(0) = -hk
        a_kick(0) =  vk
      endif
    endif
  end select

  select case (ele%key)
  case (quadrupole$, sol_quad$)
    b_kick(1) = ele%value(k1$) * ele%value(l$)

  case (sbend$)
    b_kick(0) = b_kick(0) + ele%value(dg$) * ele%value(l$)
    b_kick(1) = ele%value(k1$) * ele%value(l$)
    b_kick(2) = ele%value(k2$) * ele%value(l$) / 2.0_rp

  case (sextupole$)
    b_kick(2) = ele%value(k2$) * ele%value(l$) / 2.0_rp

  case (octupole$)
    b_kick(3) = ele%value(k3$) * ele%value(l$) / 6.0_rp
  end select

! Electric element

case default
  select case (ele%key)
  case (elseparator$)
    if (ele%value(l$) == 0) then
      b_kick(0) = 0
      a_kick(0) = 0
    else
      b_kick(0) = -ele%value(hkick$) * ele%value(p0c$) / ele%value(l$)
      a_kick(0) =  ele%value(vkick$) * ele%value(p0c$) / ele%value(l$)
    endif
  end select
end select

end subroutine this_ele_non_multipoles

!---------------------------------------------
! contains

! Tilt from element body coords to laboratory coords

subroutine tilt_this_multipole (ele, sgn, a, b, ix_pole_max)

type (ele_struct) ele
real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) tilt, an, bn, sin_t, cos_t
integer sgn, ix_pole_max, n

!

if (ele%key == sbend$) then
  tilt = sgn * ele%value(ref_tilt_tot$)
else
  tilt = sgn * ele%value(tilt_tot$) 
endif

if (tilt == 0) return

do n = 0, ix_pole_max
  if (a(n) == 0 .and. b(n) == 0) cycle
  an = a(n); bn = b(n)
  cos_t = cos((n+1)*tilt)
  sin_t = sin((n+1)*tilt)
  b(n) =  bn * cos_t + an * sin_t
  a(n) = -bn * sin_t + an * cos_t
enddo

end subroutine tilt_this_multipole

!---------------------------------------------
! contains

subroutine convert_this_ab (this_ele, p_type, a_in, b_in, this_a, this_b)

type (ele_struct) this_ele
type (branch_struct), pointer :: branch
real(rp) a_in(0:n_pole_maxx), b_in(0:n_pole_maxx), this_a(0:n_pole_maxx), this_b(0:n_pole_maxx)
real(rp) tilt, an, bn, sin_t, cos_t
integer p_type, ix_max, ref_exp
logical has_nonzero
logical a, b ! protect symbols

!

this_a = a_pole
this_b = b_pole

! All zero then we do not need to scale.
! Also if scaling is turned off

ix_max = max_nonzero(0, this_a, this_b)
if (ix_max < 0) return

! tilt?

tilt = this_ele%value(tilt_tot$) - ele%value(tilt_tot$) 
if (this_ele%key /= sbend$ .and. tilt /= 0) then
  do n = 0, n_pole_maxx
    if (this_a(n) /= 0 .or. this_b(n) /= 0) then
      an = this_a(n); bn = this_b(n)
      cos_t = cos((n+1)*tilt)
      sin_t = sin((n+1)*tilt)
      this_b(n) =  bn * cos_t + an * sin_t
      this_a(n) = -bn * sin_t + an * cos_t
    endif
  enddo
endif

! field_master = T -> Convert to normalized strength.

if (this_ele%field_master .and. this_ele%value(p0c$) /= 0 .and. &
        p_type == magnetic$ .and. this_ele%key /= sad_mult$) then
  factor = charge_of(this_ele%ref_species) * c_light / this_ele%value(p0c$)
  this_a = factor * this_a
  this_b = factor * this_b
endif

! radius = 0 defaults to radius = 1

if (.not. this_ele%scale_multipoles) return
if (p_type == electric$) then
  radius = this_ele%value(r0_elec$)
  if (radius /= 0 .and. radius /= 1) then
    factor = radius
    do n = 0, n_pole_maxx
      factor = factor / radius
      this_a(n) = factor * this_a(n)
      this_b(n) = factor * this_b(n)
    enddo
  endif
  return
endif

! normal case...

select case (this_ele%key)

case (sbend$, rbend$)
  const = this_ele%value(l$) * (this_ele%value(g$) + this_ele%value(dg$))
  ref_exp = 0

case (elseparator$, kicker$, ac_kicker$)
  if (this_ele%value(hkick$) == 0) then
    const = this_ele%value(vkick$)
  elseif (this_ele%value(vkick$) == 0) then
    const = this_ele%value(hkick$)
  else
    const = sqrt(this_ele%value(hkick$)**2 + this_ele%value(vkick$)**2)
  endif
  ref_exp = 0

case (hkicker$, vkicker$)
  const = this_ele%value(kick$)
  ref_exp = 0

case (wiggler$, undulator$)
  const = c_light * this_ele%value(b_max$) * this_ele%value(l_period$) / (pi * this_ele%value(p0c$))
  ref_exp = 0

case (quadrupole$, sol_quad$)
  const = this_ele%value(k1$) * this_ele%value(l$)
  ref_exp = 1

case (solenoid$)
  const = this_ele%value(ks$) * this_ele%value(l$)
  ref_exp = 1

case (sextupole$)
  const = this_ele%value(k2$) * this_ele%value(l$)
  ref_exp = 2
 
case (octupole$)
  const = this_ele%value(k3$) * this_ele%value(l$)
  ref_exp = 3
  
case (ab_multipole$, sad_mult$) ! multipoles do not scale
  return

case default
  call out_io (s_fatal$, r_name, 'CONFUSION! TRYING TO SCALE NON-SCALABLE ELEMENT: ' // this_ele%name)
  if (global_com%exit_on_error) call err_exit

end select

! scale multipole values

radius = this_ele%value(r0_mag$)
if (radius == 0) then
  this_a = const * this_a
  this_b = const * this_b

else
  factor = const * radius**(ref_exp+1)
  do n = 0, n_pole_maxx
    factor = factor / radius
    this_a(n) = factor * this_a(n)
    this_b(n) = factor * this_b(n)
  enddo
endif

end subroutine convert_this_ab

!---------------------------------------------
! contains

subroutine load_this_cache(cache, p_type, ix_pole_max, a, b, a_kick, b_kick)

type (multipole_cache_struct), pointer :: cache
real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx), a_kick(0:3), b_kick(0:3)
integer p_type, ix_pole_max

!

cache => ele%multipole_cache

select case (p_type)
case (magnetic$)
  cache%a_pole_mag = a
  cache%b_pole_mag = b
  cache%ix_pole_mag_max = ix_pole_max
  cache%a_kick_mag = a_kick
  cache%b_kick_mag = b_kick
  cache%ix_kick_mag_max = max_nonzero(0, a_kick, b_kick)

case default
  cache%a_pole_elec = a
  cache%b_pole_elec = b
  cache%ix_pole_elec_max = ix_pole_max
  cache%a_kick_elec = a_kick
  cache%b_kick_elec = b_kick
  cache%ix_kick_elec_max = max_nonzero(0, a_kick, b_kick)
end select

end subroutine load_this_cache

!---------------------------------------------
! contains

function pull_this_b1 (a, b, ix_pole_max) result (b1)

real(rp) a(0:), b(0:), b1
integer ix_pole_max

!

b1 = b(1)
b(1) = 0
if (ix_pole_max > 1) return
ix_pole_max = max_nonzero(0, a(0:1), b(0:1))

end function pull_this_b1

end subroutine multipole_ele_to_ab

