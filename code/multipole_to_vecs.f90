!+
! Subroutine multipole_to_vecs (ele, particle, knl, tilt)
!
! Subroutine to put the multipole components (strength and tilt)
! into 2 vectors along with the appropriate scaling.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele      -- Ele_struct: Multipole element.
!   particle -- Integer: Particle species (+1 or -1).
!
! Output:
!   knl(0:n_pole_maxx)  -- Real: Vector of strengths, MAD units.
!   tilt(0:n_pole_maxx) -- Real: Vector of tilts.
!-

subroutine multipole_to_vecs (ele, particle, knl, tilt)

  use bmad_struct
  implicit none
  type (ele_struct)  ele

  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx), signn, a_n, b_n
  real value(n_attrib_maxx), a(0:n_pole_maxx), b(0:n_pole_maxx)

  integer n, particle, n_fact

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on) then
    knl = 0
    tilt = 0
    return
  endif

! multipole
                    
  if (ele%key == multipole$) then
    knl  = ele%value(ix1_m$:ix2_m$-1:2)
    tilt = ele%value(ix1_m$+1:ix2_m$:2) + ele%value(tilt$)
    return
  endif

! ab_multiple, etc...

  call multipole_ab_scale (ele, particle, a, b)
  call multipole_ab_to_kt (a, b, knl, tilt)
  tilt = tilt + ele%value(tilt$)
                       
end subroutine
