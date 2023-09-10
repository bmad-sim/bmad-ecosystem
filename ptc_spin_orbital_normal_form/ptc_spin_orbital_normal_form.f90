!+
! Program ptc_spin_orbital_normal_form
!
! Program to construct the one turn spin/orbit map using two different methods:
!   1) By symplectic integration using an explicit Hamiltonian and 
!   2) By using PTC.
!
! A writeup of this program is in:
!    ptc_spin_orbital_normal_form.pdf
!-

program ptc_spin_orbital_normal_form

use pointer_lattice

implicit none

integer n1, n2, n3, i, j, je(6), ier
integer nd2, npara

real(dp) thin, x(6), total_time

type(real_8) y(6), y_bmad(6), z 
type(normalform) normal
type(layout), pointer :: psr
type(c_normal_form) cn
type(fibre), pointer :: p, f1, f2
type(internal_state), target:: state
type(real_8) vm, phi, b(3), e(3)
type(work) w
integer npe, neq
 
type(teapotp), pointer ::el
type(c_damap) ma, id
type(damap) id0
type(probe) xs0
type(probe_8) xs
type(c_vector_field) c_h
integer ord
type(real_8) dE, az
type(c_taylor) phase(3), ht, spin_tune
type(taylor) pot, t, h
type(real_8) dp1, n, dl_ds, gamma, beta
type(real_8) b_par(3), be, B_perp(3), E_cross_V(3), efd(3), v_unit(3), om(3)
 
! Init
use_quaternion=.true.
call ptc_ini_no_append
 
state = time0+nocavity0+spin0+totalpath0
call read_lattice_append(M_U, "g-2.flat")

! Obliterating all fringe fields (Vertical nonlinear focussing in bends) to facilitate comparison

psr =>m_u%start
p =>psr%start

do j = 1, psr%n
  p%mag%p%KILL_ENT_FRINGE = .true.
  p%mag%p%KILL_exi_FRINGE = .true.
  p%magp%p%KILL_ENT_FRINGE = .true.
  p%magp%p%KILL_exi_FRINGE = .true.
  p%mag%p%KILL_ENT_spin = .true.
  p%mag%p%KILL_exi_spin = .true.
  p%magp%p%KILL_ENT_spin = .true.
  p%magp%p%KILL_exi_spin = .true.
  p =>p%next
enddo

! Calculate number of steps for tracking for each element in the lattice

thin = 0.001
lielib_print(14) = 0  ! No printing of info messages.
call thin_lens_resplit(psr, THIN, lim = [3, 6])

! Init FPP/PTC to bea ready for tracking.

ord = 3   ! Order of DA calculation.
call init(state, ord, 0)   ! Init fpp/ptc

! Init some varialbes

call alloc(y);
call alloc(phi); 
call alloc(vm);
call alloc(b); call alloc(e)
call alloc(dE, az)
call alloc(id)
call alloc(ma)
call alloc(cn)
call alloc(phase)
call alloc(h, pot, t)
call alloc(ht, spin_tune)
call alloc(c_h)
call alloc(dp1, n)
call alloc(b_par);call alloc(be);call alloc(B_perp)
call alloc(E_cross_V);call alloc(efd);call alloc(v_unit)
call alloc(DL_DS, gamma, beta)
call alloc(om)
call alloc(xs)

! Init y to zero orbit

y(1) = 1.d0.mono.1
y(2) = 1.d0.mono.2
y(3) = 1.d0.mono.3
y(4) = 1.d0.mono.4
y(5) = 1.d0.mono.5
y(6) = 1.d0.mono.6

! Loop over all lattice elements to get one-turn spin/orbit map using an explicit Hamiltonian

p =>psr%start
w = p  ! Computes reference energy, etc.

ma = 1    ! Initialize map to be unit map.
total_time = 0

do j = 1, psr%n

  if (p%mag%kind /= kind10) then   ! Ignore markers
    p =>p%next
    cycle
  endif

  ! The analysis here assumes that the lattice only has "teapot" style elements which include
  ! sbends and quadrupoles but does not include truerbend elements

  el =>p%magp%tp10   ! Element must be "teapot" style

  ! The electric potential phi and the vm are computed by PTC (not FPP)
  ! GETELECTRIC will correct the fields and potentials if in a bend to keep Maxwell's equations satisfied.
  ! e(3) = normalized electric field
  ! phi = normalized electric potential
  ! b(3) = normalized magnetic field
  ! vm = normalized magnetic scalar potential

  call GETELECTRIC(EL, E, phi, B, VM, y, .false.)

  ! compute the az = (1+h*x) * a_s vector potential. get_az is a contained routine below.
  
  call get_az(vm, az, el%p%b0)

  ! h = hamiltonian = -(1 + g * x) * sqrt(1+2*E/beta0+E^2 - px^2 - py^2) - az
  ! Note: All bend elements are assumed to be in the horizontal plane.

  dE = y(5) - phi  ! = (E - E0) / P0
  h = -(1.0_dp+y(1)*el%p%b0) * sqrt(1.0_dp+2*dE/w%beta0+dE**2-y(2)**2-y(4)**2) - az  
 
  ! ht is -L * H, L = length. The associated map is exp(-:L*H:)Id
  ! The complex package has NO Poisson Bracket operator
  ! so convert the c_taylor ht into a c_vector_field c_h
  ! The c_vector_field has spin built into it.

  ht = -h * p%mag%tp10%l
  c_h = getvectorfield(ht)

  ! For the BMT equation need the components of E & B with respect to the direction of propagation.

  DP1 = sqrt(1.0_dp+2.0_dp*dE/w%BETA0+dE**2)
  N = sqrt(DP1**2-y(2)**2-y(4)**2)

  ! V_unit is the unit vector in the direction of propagation
  V_unit(1) = y(2)/DP1
  V_unit(2) = y(4)/DP1
  V_unit(3) = N/DP1

  ! B_par is the parallel component of B
  ! B_perp is the perpendicular component of B

  be = b(1)*v_unit(1)+b(2)*v_unit(2)+b(3)*v_unit(3)  ! Magnitude along axis parallel to propagation

  do i = 1, 3
     B_par(i) = be*v_unit(i)
  enddo

  do i = 1, 3
    B_perp(i) = B(i)-B_par(i)
  enddo

  ! E_cross_V = (direction of ray) cross (Electric field) = V_unit x E

  E_cross_V(1) = -E(2)*V_unit(3) + E(3)*V_unit(2)      
  E_cross_V(2) = -E(3)*V_unit(1) + E(1)*V_unit(3)
  E_cross_V(3) = -E(1)*V_unit(2) + E(2)*V_unit(1)

  ! This is proportional to dcT/ds (time as function of s derivative)

  DL_DS = 1.0_dp/sqrt(1.0_dp+2.0_dp*dE/w%BETA0+dE**2-y(2)**2-y(4)**2)*(1.0_dp+el%p%b0*y(1))

  ! In a sector bend the reference frame rotates continuously by el%p%b0: it must be included in the BMT equation.
  OM(2) = el%p%b0

  GAMMA = w%BETA0/w%GAMMA0I * (1.0_dp/w%BETA0 + dE)
  beta = sqrt(1.0_dp+2.0_dp*dE/w%BETA0+dE**2)/(1.0_dp/w%BETA0 + dE)  

  ! Hamiltonian with spin: c_h = F . grad + OM . L  
  ! but so far c_h%om = 0  %om contains the spin operator.

  ! This is the magnetic part of the BMT equation dOmega/ds

  om(1) = -DL_DS*( (1.0_dp+p%AG*GAMMA)*B_perp(1) + (1.0_dp+p%AG)*B_par(1) )
  om(2) = -DL_DS*( (1.0_dp+p%AG*GAMMA)*B_perp(2) + (1.0_dp+p%AG)*B_par(2) )+OM(2)
  om(3) = -DL_DS*( (1.0_dp+p%AG*GAMMA)*B_perp(3) + (1.0_dp+p%AG)*B_par(3) )

  ! This is the electric part of the BMT equation

  do I = 1, 3
     OM(I) = OM(I) + DL_DS*beta*gamma*(p%AG+1.0_dp/(1.0_dp+GAMMA))*E_cross_V(I)
  enddo
 
  ! c_h will now contain the spin effect

  c_h%q%x(0)=0.0_dp
  do i = 1, 3
    t = om(i)
    c_h%q%x(i) = t * p%mag%tp10%l
  enddo

  ! Here I am "COSY-INFINITIZING" (exponanting the Lee operator) to the get the one-turn map with spin.
  ! Note: The "*" operator in ma = id * ma ignores constant terms. 
  ! [The operation ma = id .o. ma would not ignore constant terms]
  ! Since the zeroth order terms in the map gets ignored, and we need to keep track of 
  ! the conatant part of the total time separately.

  id = 1
  id = exp(c_h, id)
  ma = id * ma
  total_time = total_time + (id%v(6) .sub. "0")
  p => p%next
enddo

! Put the constant part of the total time back in the map

ma%v(6) = ma%v(6) - (ma%v(6) .sub. '0') + total_time

! The map is normalized to check chromaticities

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'Output from call to c_normal:'
call c_normal (ma, cn, dospin = .true., phase = phase, nu_spin = spin_tune)

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'Tunes:'
write(6, *) cn%tune(1:3)


print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'phase(1) map. Independent variables are (x+I*px, x-I*px, y+I*py, y-I*py, E, time)'
print '(a)', 'Remember: The map is complex so each term has real and imaginary parts.'
print '(a)', 'The real part of the "000000" term is the tune of the horizontal like mode.'
print '(a)', 'The real part of the "000010" term is the chromaticity of the horizontal like mode.'

call print(phase(1), 6)

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'phase(2) map. Independent variables are (x+I*px, x-I*px, y+I*py, y-I*py, E, time)'
print '(a)', 'Remember: The map is complex so each term has real and imaginary parts.'
print '(a)', 'The real part of the "000000" term is the tune of the vertical like mode.'
print '(a)', 'The real part of the "000010" term is the chromaticity of the vertical like mode.'

call print(phase(2), 6)

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'phase(3) map.'

phase(3) = phase(3) + (ma%v(6).sub.'0')  ! Add in total time constant
call print(phase(3), 6)

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'spin_tune map. Independent variables are (x+I*px, x-I*px, y+I*py, y-I*py, E, time)'

call print(spin_tune, 6)

!-----------------------------------------------
! Get the spin orbital map using PTC.

id = 1
x = 0.d0
xs0 = x
xs = xs0+id

call track_probe (psr, xs, state, fibre1 = 1)

id = xs

! I substract the two maps and print them
! Since we use a second order integrator for spin, I do not expect prefection on the spin.
! It should be perfect on the orbital part.

ma = ma - id
ma = ma .cut. (ord-1)

print '(a)'
print '(a)', '!-------------------------------------------------------'
print '(a)', 'Differnece between direct Hamiltonian tracking and PTC tracking.'
print '(a)', 'The first 6 Taylor series are the orbital part and the last 9 are the spin part.'

call print(ma, 6)

!----------------------------------------------------------------------
contains

subroutine get_az(vm, az, h)
implicit none
TYPE(real_8) vm, az
TYPE(taylor)v(2), a, vmt
type(c_taylor) cc
integer, allocatable :: je(:)
integer j
complex(dp) vvc
real(dp) vv, h
allocate(je(c_%nv))


call alloc(v)
call alloc(a, vmt)
call alloc(cc)
az = 0.0_dp
vmt = vm
v(1) = (vmt .d. 3) * (1.d0 + h*(1.d0 .mono. 1))
v(2) = -(vmt .d. 1) * (1.d0 + h*(1.d0 .mono. 1))

j = 1

do while(.true.) 
  cc = v(1)
  call  c_cycle(cc, j, vvc , je); if (j == 0) exit;
  vv = vvc
  az = az+morph(((vv/(je(1)+je(3)+1)) .mono. je)*(1.0_dp .mono. 1))
enddo

j = 1

do while(.true.) 
  cc = v(2)
  call  c_cycle(cc, j, vvc , je); if (j == 0) exit;
  vv = vvc
  az = az+morph(((vv/(je(1)+je(3)+1)) .mono. je)*(1.0_dp .mono. 3))
enddo

call kill(v)
call kill(a, vmt)
call kill(cc)
deallocate(je)

end subroutine get_az

end program ptc_spin_orbital_normal_form
