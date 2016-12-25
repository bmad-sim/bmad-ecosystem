module multipole_mod

use bmad_utils_mod

implicit none

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ab_to_kt (an, bn, knl, tn)
!
! Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.
! Also see: multipole1_ab_to_kt.
!
! Input:
!   an(0:n_pole_maxx) -- Real(rp): Skew multipole component.
!   bn(0:n_pole_maxx) -- Real(rp): Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real(rp): Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real(rp): Multipole angle.
!-

subroutine multipole_ab_to_kt (an, bn, knl, tn)

real(rp) an(0:), bn(0:)
real(rp) knl(0:), tn(0:)

integer n

!

do n = 0, n_pole_maxx
  call multipole1_ab_to_kt (an(n), bn(n), n, knl(n), tn(n))
enddo

end subroutine multipole_ab_to_kt

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole1_ab_to_kt (an, bn, n, knl, tn)
!
! Subroutine to convert ab type multipole to kt (MAD standard) multipole.
! Also see: multipole_ab_to_kt.
!
! Input:
!   an -- Real(rp): Skew multipole component.
!   bn -- Real(rp): Normal multipole component.
!   n  -- Integer: Order of multipole. 
!
! Output:
!   knl -- Real(rp): Multitude magnatude.
!   tn  -- Real(rp): Multipole angle.
!-

subroutine multipole1_ab_to_kt (an, bn, n, knl, tn)

real(rp) an, bn, knl, tn

integer n

!

real(rp) a, b

if (an == 0 .and. bn == 0) then
  knl = 0
  tn = 0
else
  ! Use temp a, b to avoid problems when actual (knl, tn) args are the same as (an, bn).
  a = an
  b = bn
  knl  = factorial(n) * sqrt(a**2 + b**2)
  tn = -atan2(a, b) / (n + 1)
  ! In case the user looks at this, make tn to be in the range [-pi, pi]/(n+1)
  if ((n + 1) * abs(tn) > pi) then
    knl = -knl
    tn = sign_of(tn) * (abs(tn) - pi/(n+1))
  endif

endif

end subroutine multipole1_ab_to_kt

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ele_to_kt (ele, use_ele_tilt, ix_pole_max, knl, tilt, pole_type)
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
!
! Output:
!   ix_pole_max         -- Integer: Index of largest nonzero pole.
!   knl(0:n_pole_maxx)  -- Real(rp): Vector of strengths, MAD units.
!   tilt(0:n_pole_maxx) -- Real(rp): Vector of tilts.
!-

subroutine multipole_ele_to_kt (ele, use_ele_tilt, ix_pole_max, knl, tilt, pole_type)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord

real(rp) knl(0:), tilt(0:), a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) this_a(0:n_pole_maxx), this_b(0:n_pole_maxx)
real(rp) tilt1
real(rp), pointer :: a_pole(:), b_pole(:)

integer ix_pole_max
integer, optional :: pole_type
integer i, ix_max

logical use_ele_tilt

! Init

ix_pole_max = -1
call pointer_to_ele_multipole (ele, a_pole, b_pole, pole_type)

! Note: For multipoles, use_ele_tilt arg is ignored.
! Slice slaves and super slaves have their associated multipoles stored in the lord.

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    call pointer_to_ele_multipole (lord, a_pole, b_pole, pole_type)
    if (.not. associated(a_pole)) cycle
    call multipole_ele_to_ab (lord, use_ele_tilt, ix_max, this_a, this_b, pole_type)
    if (ix_max == -1) cycle
    if (ix_pole_max == -1) then
      a = this_a * (ele%value(l$) / lord%value(l$))
      b = this_b * (ele%value(l$) / lord%value(l$))
    else
      a = a + this_a * (ele%value(l$) / lord%value(l$))
      b = b + this_b * (ele%value(l$) / lord%value(l$))
    endif
    ix_pole_max = max(ix_pole_max, ix_max)
  enddo

! Not a slave

else
  if (.not. associated(a_pole)) return
  if (ele%key == multipole$) then
    knl  = a_pole
    tilt = b_pole + ele%value(tilt_tot$)
    ix_pole_max = max_nonzero(0, knl)
    return
  else
    call multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type)
  endif
endif

!

if (ix_pole_max > -1) then
  call multipole_ab_to_kt (a, b, knl, tilt)
endif

end subroutine multipole_ele_to_kt

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kt_to_ab (knl, tn, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
! Also see: multipole1_kt_to_ab.
!
! Input:
!   knl(0:) -- Real(rp): Multitude magnatude.
!   tn(0:)  -- Real(rp): Multipole angle.
!
! Output:
!   an(0:) -- Real(rp): Skew multipole component.
!   bn(0:) -- Real(rp): Normal multipole component.
!-

subroutine multipole_kt_to_ab (knl, tn, an, bn)

real(rp) an(0:), bn(0:)
real(rp) knl(0:), tn(0:)

integer n

!

do n = lbound(an, 1), ubound(an, 1)
  call multipole1_kt_to_ab (knl(n), tn(n), n, an(n), bn(n))
enddo

end subroutine multipole_kt_to_ab

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole1_kt_to_ab (knl, tn, n, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
! Also see: multipole_kt_to_ab.
!
! Input:
!   knl -- Real(rp): Multitude magnatude.
!   tn  -- Real(rp): Multipole angle.
!   n   -- Integer: Multipole order.
!
! Output:
!   an -- Real(rp): Skew multipole component.
!   bn -- Real(rp): Normal multipole component.
!-

subroutine multipole1_kt_to_ab (knl, tn, n, an, bn)

real(rp) an, bn
real(rp) knl, tn
real(rp) angle, kl

integer n

!

if (knl == 0) then
  an = 0
  bn = 0
else
  kl = knl / factorial(n)
  angle = -tn * (n + 1)
  an = kl * sin(angle)
  bn = kl * cos(angle)
endif

end subroutine multipole1_kt_to_ab

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks)
!                             
! Subroutine to extract the ab multipole values of an element.
! Note: The ab values will be scalled by the strength of the element.
!
! Input:
!   ele           -- ele_struct: Element.
!     %value()      -- ab_multipole values.
!   use_ele_tilt  -- logical: If True then include ele%value(tilt_tot$) in calculations.
!                      use_ele_tilt is ignored in the case of multipole$ elements.
!   pole_type     -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   include_kicks -- logical, optional: If True, include hkick/vkick in the n = 0 components.
!                      Default is False.
!
! Output:
!   ix_pole_max      -- Integer: Index of largest nonzero pole.
!   a(0:n_pole_maxx) -- real(rp): Array of scalled multipole values.
!   b(0:n_pole_maxx) -- real(rp): Array of scalled multipole values.
!-

subroutine multipole_ele_to_ab (ele, use_ele_tilt, ix_pole_max, a, b, pole_type, include_kicks)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord

real(rp) const, radius, factor, a(0:), b(0:)
real(rp) an, bn, cos_t, sin_t, tilt, hk, vk
real(rp) this_a(0:n_pole_maxx), this_b(0:n_pole_maxx)
real(rp), pointer :: a_pole(:), b_pole(:)

integer ix_pole_max
integer, optional :: pole_type
integer i, ref_exp, n, tilt_dir

logical, optional :: include_kicks
logical use_ele_tilt

character(*), parameter :: r_name = 'multipole_ele_to_ab'

! Init

a = 0
b = 0
ix_pole_max = -1

call pointer_to_ele_multipole (ele, a_pole, b_pole, pole_type)

! Multipole type element case. Note: use_ele_tilt is ignored in this case.

! Slice slaves and super slaves have their associated multipoles stored in the lord

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%key == group$ .or. lord%key == overlay$ .or. lord%key == girder$) cycle
    call pointer_to_ele_multipole (lord, a_pole, b_pole, pole_type)
    if (.not. associated(a_pole)) cycle
    if (.not. (lord%multipoles_on .and. lord%is_on)) cycle

    if (lord%key == multipole$) then
      ix_pole_max = max_nonzero(0, a_pole)
      if (ix_pole_max == -1) return
      call multipole_kt_to_ab (a_pole, b_pole, a, b)
      return
    endif

    call convert_this_ab (lord, this_a, this_b)
    a = a + this_a * (ele%value(l$) / lord%value(l$))
    b = b + this_b * (ele%value(l$) / lord%value(l$))
  enddo

  ix_pole_max = max_nonzero(0, a, b)

! Not a slave
else
  if (ele%key == multipole$) then
    if (integer_option(magnetic$, pole_type) == electric$) return
    ix_pole_max = max_nonzero(0, a_pole, b_pole)
    if (ix_pole_max == -1) return
    call multipole_kt_to_ab (a_pole, b_pole, a, b)
  else
    call convert_this_ab (ele, a, b)
  endif
endif

! Include h/v kicks?

if (logic_option(.false., include_kicks) .and. ele%is_on) then
  hk = 0; vk = 0

  if (integer_option(magnetic$, pole_type) == magnetic$) then
    select case (ele%key)
    case (hkicker$)
      tilt_dir = 0
      hk = -ele%value(kick$)
    case (vkicker$)
      tilt_dir = 0
      vk = ele%value(kick$)
    case (elseparator$)
      ! Kicks are electric
    case default
      tilt_dir = -1
      hk = -ele%value(hkick$)
      vk =  ele%value(vkick$)
    end select

  else ! electric
    select case (ele%key)
    case (elseparator$)
      tilt_dir = 0
      hk = -ele%value(hkick$) * ele%value(p0c$) / ele%value(l$)
      vk =  ele%value(vkick$) * ele%value(p0c$) / ele%value(l$)
    end select
  endif

  if (hk /= 0 .or. vk /= 0) then
    if (use_ele_tilt) tilt_dir = tilt_dir + 1

    if (ele%key == sbend$) then
      tilt = tilt_dir * ele%value(ref_tilt_tot$)
    else
      tilt = tilt_dir * ele%value(tilt_tot$) 
    endif

    if (tilt /= 0) then
      cos_t = cos(tilt)
      sin_t = sin(tilt)
      b(0) = b(0) + hk * cos_t + vk * sin_t
      a(0) = a(0) - hk * sin_t + vk * cos_t
    else
      b(0) = b(0) + hk 
      a(0) = a(0) + vk
    endif
  endif

endif

ix_pole_max = max_nonzero(0, a, b)

!---------------------------------------------
contains

subroutine convert_this_ab (this_ele, this_a, this_b)

type (ele_struct) this_ele
type (branch_struct), pointer :: branch
real(rp) this_a(0:n_pole_maxx), this_b(0:n_pole_maxx), tilt
integer ix_max
logical has_nonzero
logical a, b ! protect symbols

!

if (.not. associated(a_pole)) return
if (.not. (this_ele%multipoles_on .and. this_ele%is_on)) return

this_a = a_pole
this_b = b_pole

! All zero then we do not need to scale.
! Also if scaling is turned off

ix_max = max_nonzero(0, this_a, this_b)
if (ix_max == -1) return

! use tilt?

if (this_ele%key == sbend$) then
  tilt = ele%value(ref_tilt_tot$)
else
  tilt = this_ele%value(tilt_tot$) 
endif

if (use_ele_tilt .and. tilt /= 0) then
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
        integer_option(magnetic$, pole_type) == magnetic$ .and. this_ele%key /= sad_mult$) then
  branch => pointer_to_branch(this_ele)
  if (.not. associated(branch)) then
    call out_io (s_fatal$, r_name, 'ELEMENT WITH MULTIPOLES AND FIELD_MASTER = T NOT ASSOCIATED WITH ANY LATTICE!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
  factor = charge_of(branch%param%particle) * c_light / this_ele%value(p0c$)
  this_a = factor * this_a
  this_b = factor * this_b
endif

! radius = 0 defaults to radius = 1

if (.not. this_ele%scale_multipoles) return
if (integer_option(magnetic$, pole_type) == electric$) then
  radius = this_ele%value(r0_elec$)
  if (radius /= 0) then
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
  const = this_ele%value(l$) * (this_ele%value(g$) + this_ele%value(g_err$))
  ref_exp = 0

case (elseparator$, kicker$)
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
  const = 2 * c_light * this_ele%value(b_max$) * this_ele%value(l_pole$) / (pi * this_ele%value(p0c$))
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

end subroutine multipole_ele_to_ab

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kicks (knl, tilt, ref_species, coord, pole_type, ref_orb_offset)
!
! Subroutine to put in the kick due to a multipole.
!
! Input:
!   knl(0:)        -- real(rp): Multipole strengths (mad units).
!   tilt(0:)       -- real(rp): Multipole tilts.
!   ref_species    -- integer: Reference species.
!   coord          -- coord_struct: Particle position.
!   pole_type      -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   ref_orb_offset -- logical, optional: If present and n = 0 then the
!                       multipole simulates a zero length bend with bending
!                       angle knl.
!
! Output:
!   coord -- coord_struct: Kicked particle.
!-

subroutine multipole_kicks (knl, tilt, ref_species, coord, pole_type, ref_orb_offset)

type (coord_struct)  coord

real(rp) knl(0:), tilt(0:)

integer ref_species, n

integer, optional :: pole_type
logical, optional :: ref_orb_offset

!

do n = 0, n_pole_maxx
  if (knl(n) == 0) cycle
  call multipole_kick (knl(n), tilt(n), n, ref_species, coord, pole_type, ref_orb_offset)
enddo

end subroutine multipole_kicks

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ab_multipole_kicks (an, bn, ref_species, coord, pole_type, scale, mat6, make_matrix)
!
! Routine to put in the kick due to ab_multipole components.
!
! Input:
!   an(0:)       -- real(rp): Skew multipole strengths.
!   bn(0:)       -- real(rp): Normal multipole tilts.
!   ref_species    -- integer: Reference species.
!   coord        -- coord_struct: Particle position.
!   pole_type    -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   scale        -- real(rp), optional: Factor to scale the kicks. Default is 1.
!                     For pole_type = electric$, set scale to the longitudinal length of the field region
!   mat6(6,6)    -- Real(rp), optional: Transfer matrix before the multipole.
!   make_matrix  -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   coord      -- coord_struct: Kicked particle.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including multipole.
!-

subroutine ab_multipole_kicks (an, bn, ref_species, coord, pole_type, scale, mat6, make_matrix)

type (coord_struct)  coord

real(rp) an(0:), bn(0:)
real(rp) kx, ky, pz2, rel_p2, dk(2,2), alpha, kx_tot, ky_tot
real(rp), optional :: scale, mat6(6,6)

integer, optional :: pole_type
integer ref_species, n

logical, optional :: make_matrix

!

if (integer_option(magnetic$, pole_type) == electric$) then
  pz2 = (1 + coord%vec(6))**2 - coord%vec(2)**2 - coord%vec(4)**2
endif

kx_tot = 0
ky_tot = 0

do n = 0, n_pole_maxx
  if (an(n) == 0 .and. bn(n) == 0) cycle

  if (logic_option(.false., make_matrix)) then
    call ab_multipole_kick (an(n), bn(n), n, ref_species, coord, kx, ky, dk, pole_type = pole_type, scale = scale)
    mat6(2,:) = mat6(2,:) + dk(1,1) * mat6(1,:) + dk(1,2) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + dk(2,1) * mat6(1,:) + dk(2,2) * mat6(3,:)
  else
    call ab_multipole_kick (an(n), bn(n), n, ref_species, coord, kx, ky, pole_type = pole_type, scale = scale)
  endif

  kx_tot = kx_tot + kx
  ky_tot = ky_tot + ky
enddo

coord%vec(2) = coord%vec(2) + kx_tot
coord%vec(4) = coord%vec(4) + ky_tot

if (integer_option(magnetic$, pole_type) == electric$) then
  alpha = (kx_tot * (2*coord%vec(2) - kx_tot) + ky_tot * (2*coord%vec(4) - ky_tot)) / (1 + coord%vec(6))**2
  if (alpha < -1) then
    coord%state = lost_z_aperture$
    return
  endif
  coord%vec(6) = coord%vec(6) + (1 + coord%vec(6)) * sqrt1(alpha)
  rel_p2 = pz2 + coord%vec(2)**2 + coord%vec(4)**2
  coord%beta = (1 + coord%vec(6)) / sqrt(rel_p2 + (mass_of(coord%species)/coord%p0c)**2)
endif

end subroutine ab_multipole_kicks

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kick (knl, tilt, n, ref_species, coord, pole_type, ref_orb_offset)
!
! Subroutine to put in the kick due to a multipole.
! Note: The kick for an electric multipole does not include any energy change.
!
! Input:
!   knl            -- real(rp): Multipole integrated strength.
!   tilt           -- real(rp): Multipole tilt.
!   n              -- real(rp): Multipole order.
!   ref_species    -- integer: Reference species.
!   coord          -- coord_struct:
!     %vec(1)         -- X position.
!     %vec(3)         -- Y position.
!   pole_type      -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   ref_orb_offset -- logical, optional: If True and n = 0 then the
!                       multipole simulates a zero length bend with bending
!                       angle knl. Default is False.
!
! Output:
!   coord -- coord_struct: 
!     %vec(2) -- X kick.
!     %vec(4) -- Y kick.
!-

subroutine multipole_kick (knl, tilt, n, ref_species, coord, pole_type, ref_orb_offset)

type (coord_struct) coord

real(rp) knl, tilt, x, y, sin_ang, cos_ang
real(rp) x_vel, y_vel, charge
real(rp) x_value, y_value
real(rp) cval, rp_dummy, t
real(rp) x_terms(0:n)
real(rp) y_terms(0:n)
real(rp), save :: cc(0:n_pole_maxx, 0:n_pole_maxx)

logical, save :: first_call = .true.

integer, optional :: pole_type
integer ref_species, n, m

logical, optional :: ref_orb_offset

! simple case

if (knl == 0) return

! normal case

t = tilt
if (integer_option(magnetic$, pole_type) == magnetic$) then
  charge = charge_to_mass_of(coord%species) / charge_to_mass_of(ref_species)
else
  charge = charge_of(coord%species) / (coord%beta * coord%p0c)
  t = pi/(n+1) - t
endif

if (t == 0) then
  sin_ang = 0
  cos_ang = 1
  x = coord%vec(1)
  y = coord%vec(3)
else
  sin_ang = sin(t)
  cos_ang = cos(t)
  x =  coord%vec(1) * cos_ang + coord%vec(3) * sin_ang
  y = -coord%vec(1) * sin_ang + coord%vec(3) * cos_ang
endif

! ref_orb_offset with n = 0 means that we are simulating a zero length dipole.

if (n == 0 .and. logic_option(.false., ref_orb_offset)) then
  coord%vec(2) = coord%vec(2) + charge * knl * cos_ang * coord%vec(6)
  coord%vec(4) = coord%vec(4) + charge * knl * sin_ang * coord%vec(6)
  coord%vec(5) = coord%vec(5) - charge * knl * (cos_ang * coord%vec(1) + sin_ang * coord%vec(3))
  return
endif

! normal case

x_terms(n)=1.0
y_terms(0)=1.0
do m = 1, n
  x_terms(n-m) = x_terms(n-m+1)*x
  y_terms(m) = y_terms(m-1)*y
enddo

if (first_call) then
  ! populate cc 
  rp_dummy = c_multi(0,0,c_full=cc)
  first_call = .false.
endif

x_value = SUM(cc(n,0:n:2) * x_terms(0:n:2) * y_terms(0:n:2))
y_value = SUM(cc(n,1:n:2) * x_terms(1:n:2) * y_terms(1:n:2))

x_vel = charge * knl * x_value
y_vel = charge * knl * y_value

if (t == 0) then
  coord%vec(2) = coord%vec(2) + x_vel
  coord%vec(4) = coord%vec(4) + y_vel
else
  coord%vec(2) = coord%vec(2) + x_vel * cos_ang - y_vel * sin_ang
  coord%vec(4) = coord%vec(4) + x_vel * sin_ang + y_vel * cos_ang
endif

end subroutine multipole_kick

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ab_multipole_kick (a, b, n, ref_species, coord, kx, ky, dk, pole_type, scale)
!
! Subroutine to put in the kick due to an ab_multipole.
!                          
! Input:
!   a            -- Real(rp): Multipole skew component.
!   b            -- Real(rp): Multipole normal component.
!   n            -- Real(rp): Multipole order.
!   ref_species  -- integer: Reference species.
!   coord        -- Coord_struct:
!   pole_type    -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!   scale        -- real(rp), optional: Factor to scale the kicks. Default is 1.
!                     For pole_type = electric$, set scale to the longitudinal length of the field region.
!
! Output:
!   kx      -- Real(rp): X kick.
!   ky      -- Real(rp): Y kick.
!   dk(2,2) -- Real(rp), optional: Kick derivative: dkick(x,y)/d(x,y).
!-

subroutine ab_multipole_kick (a, b, n, ref_species, coord, kx, ky, dk, pole_type, scale)

type (coord_struct)  coord

real(rp) a, b, x, y
real(rp), optional :: dk(2,2)
real(rp), optional :: scale
real(rp) kx, ky, f, a2, b2

integer, optional :: pole_type
integer ref_species, n, m, n1

! Init

kx = 0
ky = 0

if (present(dk)) dk = 0

! simple case

if (a == 0 .and. b == 0) return

! normal case
! Note that c_multi can be + or -

if (integer_option(magnetic$, pole_type) == electric$) then
  f = charge_of(coord%species) / (coord%beta * coord%p0c)
  a2 =  a * f
  b2 = -b * f
else
  f = charge_to_mass_of(coord%species) / charge_to_mass_of(ref_species)
  a2 = a * f
  b2 = b * f
endif

if (present(scale)) then
  a2 = scale * a2
  b2 = scale * b2
endif

x = coord%vec(1)
y = coord%vec(3)

do m = 0, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  kx = kx + b2 * f
  ky = ky - a2 * f
enddo

do m = 1, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  kx = kx + a2 * f
  ky = ky + b2 * f
enddo

! dk calc

if (present(dk)) then
  n1 = n - 1
  
  do m = 0, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dk(1,1) = dk(1,1) + b2 * f
    dk(2,1) = dk(2,1) - a2 * f

    dk(1,2) = dk(1,2) - a2 * f
    dk(2,2) = dk(2,2) - b2 * f
  enddo


  do m = 1, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dk(1,2) = dk(1,2) + b2 * f
    dk(2,2) = dk(2,2) - a2 * f

    dk(1,1) = dk(1,1) + a2 * f
    dk(2,1) = dk(2,1) + b2 * f
  enddo
endif

end subroutine ab_multipole_kick

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine elec_multipole_field (a, b, n, coord, Ex, Ey, dE, compute_dE)
!
! Subroutine to put in the field due to an electric_multipole.
!
! Input:
!   a       -- Real(rp): Multipole skew component.
!   b       -- Real(rp): Multipole normal component.
!   n       -- Real(rp): Multipole order.
!   coord   -- Coord_struct:
!
! Output:
!   Ex          -- Real(rp): X field component
!   Ey          -- Real(rp): Y field component.
!   dE(2,2)     -- Real(rp), optional: Field derivatives: dfield(x,y)/d(x,y).
!   compute_dE  -- logical, optional: If False, do not compute the field derivatives 
!                     even if dE is present. Default is True.
!-

subroutine elec_multipole_field (a, b, n, coord, Ex, Ey, dE, compute_dE)

type (coord_struct)  coord

real(rp) a, b, x, y
real(rp), optional :: dE(2,2)
real(rp) Ex, Ey, f

integer n, m, n1
logical, optional :: compute_dE
logical compute

! Init

Ex = 0
Ey = 0
compute = (present(dE) .and. logic_option(.true., compute_dE))

if (compute) dE = 0

! simple case

if (a == 0 .and. b == 0) return

! normal case
! Note that c_multi can be + or -

x = coord%vec(1)
y = coord%vec(3)

do m = 0, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  Ex = Ex - b * f
  Ey = Ey - a * f
enddo

do m = 1, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  Ex = Ex + a * f
  Ey = Ey - b * f
enddo

! dE calc

if (compute) then

  n1 = n - 1
  
  do m = 0, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dE(1,1) = dE(1,1) - b * f
    dE(2,1) = dE(2,1) - a * f

    dE(1,2) = dE(1,2) - a * f
    dE(2,2) = dE(2,2) + b * f
  enddo


  do m = 1, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dE(1,2) = dE(1,2) - b * f
    dE(2,2) = dE(2,2) - a * f

    dE(1,1) = dE(1,1) + a * f
    dE(2,1) = dE(2,1) - b * f
  enddo

endif

end subroutine elec_multipole_field

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine pointer_to_ele_multipole (ele, a_pole, b_pole, pole_type)
!
! Routine to point to the appropriate magnetic or electric poles in an element.
!
! Input:
!   ele          -- Ele_struct: Lattice element.
!   pole_type    -- integer, optional: Type of multipole. magnetic$ (default) or electric$.
!
! Output:
!   a_pole(:)    -- real(rp), pointer: Pointer to skew electric or magnetic poles.
!   b_pole(:)    -- real(rp), pointer: Pointer to normal electric or magnetic poles.
!-

subroutine pointer_to_ele_multipole (ele, a_pole, b_pole, pole_type)

type (ele_struct), target :: ele

real(rp), pointer :: a_pole(:), b_pole(:)
integer, optional :: pole_type

!

select case (integer_option(magnetic$, pole_type))
case (magnetic$)
  a_pole => ele%a_pole
  b_pole => ele%b_pole

case (electric$)
  if (ele%key == multipole$) then
    a_pole => null()   ! Multipoles do not have electric fields.
    b_pole => null()
  else
    a_pole => ele%a_pole_elec
    b_pole => ele%b_pole_elec
  endif
end select

end subroutine pointer_to_ele_multipole

end module
