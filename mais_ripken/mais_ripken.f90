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
write (twiss_file, '(a12,20a15)' )  "#N", "NAME", "S[m]", "BetaX1[m]", "AlphaX1", "BetaY1[m]", "AlphaY1", "Nu1/(2*PI)", &
                                                          "BetaX2[m]", "AlphaX2", "BetaY2[m]", "AlphaY2", "Nu2/(2*PI)", &
                                                          "DspX[m]", "DspXp", "DspY[m]", "DspYp", "U", "Q1", "Q2", "M56[m]"

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
       
  u = 1 - ele0%gamma_c**2
  v1 = atan2(-ele0%c_mat(1,2)/sqrt(ele0%a%beta), -ele0%c_mat(2,2)*sqrt(ele0%a%beta) - ele0%c_mat(1,2)*ele0%a%alpha/sqrt(ele0%a%beta))
  v2 = atan2(-ele0%c_mat(1,2)/sqrt(ele0%b%beta), ele0%c_mat(1,1)*sqrt(ele0%b%beta) - ele0%c_mat(1,2)*ele0%b%alpha/sqrt(ele0%b%beta))
 
  write(twiss_file, '(i10, 2x, a15, f15.6, 3(f15.4, f15.6, f15.4, 2f15.6), 2f15.6, f15.4)') ele0%ix_ele, ele0%name, s, &
                                        b11, a11, b21, a21, real(v1,rp)/twopi, &
                                        b12, a12, b22, a22, real(v2,rp)/twopi, &
                                        ele0%a%eta, ele0%a%etap, ele0%b%eta, ele0%b%etap, u, ele0%a%phi/twopi, ele0%b%phi/twopi, ele0%mat6(5,6)
end do

close(twiss_file)
end program
