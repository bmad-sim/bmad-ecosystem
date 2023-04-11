!+
! Program to calculate the Mais-Ripken coupled Twiss parameters 
! (as opposed to the Edwards-Teng formalism used by Bmad). 
!
! Henry Lovelace III 
! Program Based on Scott's tech note.
!-

program coupling_convert_et2rip

use bmad
use twiss_and_track_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele_n
type (ele_struct), target :: ele0
type (coord_struct), allocatable :: orbit(:)
type (coord_struct), target :: parti, partf, orb
type (branch_struct), pointer :: branch
integer :: ix_branch = 0
integer:: n, scount, track_state
integer:: twiss_file = 12, open_file = 13


real(rp) b11, b12, b21, b22, a11, a12, a21, a22, step
real(rp) kx, ky, Ax, Ay, u_plus, u_minus
logical err_flag

real(rp) total_length, length, s, s_out, u
complex v1, v2, num_v1, num_v2, denom_v1, denom_v2

character(120) lat_file_name, infile_name, twiss_name

namelist / input_params / lat_file_name, step, ix_branch

! Read in parameters

infile_name = 'coupling.init'
twiss_name = 'twiss.txt'
open(open_file, file = infile_name)
read (open_file, nml = input_params)
close (open_file)

!

call bmad_parser (lat_file_name, lat)
call twiss_and_track(lat, orbit)

branch => lat%branch(ix_branch)       

scount = int(lat%param%total_length / step)
print *, 'Number of points: ', scount

branch => lat%branch(ix_branch)

open(twiss_file, file = twiss_name)
write (twiss_file, '(a12,20a15)' )  "#N","NAME","S[m]","BetaX1[cm]","AlphaX1","BetaY1[cm]","AlphaY1","Nu1/(2*PI)","BetaX2[cm]","AlphaX2","BetaY2[cm]","AlphaY2","Nu2/(2*PI)","U","DspX[cm]","DspXp","DspY[cm]","DspYp","Q1","Q2","M56[cm]"

do n = 0, scount
  s = n*step + branch%ele(0)%s
  call twiss_and_track_at_s(lat, s, ele0, orbit, orb, ix_branch, err_flag, .false., .false.)
  !Matrix elements
  b11 = ele0%gamma_c**2 * ele0%a%beta
  a11 = ele0%gamma_c**2 * ele0%a%alpha 
  b22 = ele0%gamma_c**2 * ele0%b%beta
  a22 = ele0%gamma_c**2 * ele0%b%alpha
  b21 = ele0%c_mat(2,2)**2*ele0%a%beta + 2 * ele0%c_mat(1,2)*ele0%c_mat(2,2)*ele0%gamma_c * ele0%a%alpha + ele0%c_mat(1,2)**2*ele0%a%gamma
  a21 = (ele0%c_mat(1,1)*ele0%c_mat(2,2)+ele0%c_mat(1,2)*ele0%c_mat(2,1))*ele0%a%alpha + ele0%c_mat(2,1)*ele0%c_mat(2,2)*ele0%a%beta+ele0%c_mat(1,1)*ele0%c_mat(1,2)*ele0%a%gamma
  b12 = ele0%c_mat(1,1)**2*ele0%b%beta - 2 * ele0%c_mat(1,1)*ele0%c_mat(1,2)*ele0%gamma_c * ele0%b%alpha + ele0%c_mat(1,2)**2*ele0%b%gamma
  a12 = (ele0%c_mat(1,1)*ele0%c_mat(2,2)+ele0%c_mat(1,2)*ele0%c_mat(2,1))*ele0%b%alpha - ele0%c_mat(1,1)*ele0%c_mat(2,1)*ele0%b%beta - ele0%c_mat(1,2)*ele0%c_mat(2,2)*ele0%b%gamma
  

  kx = sqrt(b21/b11)
  ky = sqrt(b12/b22)
  Ax = kx * a11 - a21 / kx
  Ay = ky * a22 - a12 / ky
  u_plus = ((-kx**2*ky**2)+sqrt(kx**2*ky**2*(1+(Ax**2-Ay**2)/(kx**2-ky**2) * (1-kx**2*ky**2))))/(1-kx**2*ky**2)
  u_minus = ((-kx**2*ky**2)-sqrt(kx**2*ky**2*(1+(Ax**2-Ay**2)/(kx**2-ky**2) * (1-kx**2*ky**2))))/(1-kx**2*ky**2)  
  if (u_plus <= 1) then
    u = u_plus
  else
    u = u_minus 
  endif
  
  if (kx == 0 .and. ky == 0) then
    u = 0
  endif
  num_v1 = cmplx(Ax, (kx*(1-u)+u/kx))
  denom_v1 = cmplx(Ay, -(ky*(1-u)+u/ky))
  num_v2 = cmplx(Ax, (kx*(1-u)-u/kx))
  denom_v2 = cmplx(Ay, (ky*(1-u)-u/ky))
  

  v1 = (num_v1/denom_v1)
  v2 = (num_v2/denom_v2)
  write(twiss_file, '(i10, 2x, a15, 20f15.6)') ele0%ix_ele, ele0%name, s, b11*100, a11, b21*100, a21, v1%re, b12*100, a12, b22*100, a22, v2%re, u, ele0%a%eta*100, ele0%a%etap, ele0%b%eta*100, ele0%b%etap,ele0%a%phi/twopi, ele0%b%phi/twopi, ele0%mat6(5,6)*100

end do

close(twiss_file)
end program
