program spin_matching

use ptc_layout_mod
use f95_lapack

implicit none

type (lat_struct), target :: lat
type (spin_matching_struct), allocatable, target :: match_info(:)
type (spin_matching_struct), pointer :: sm
type (branch_struct), pointer :: branch

real(rp) s_1turn(3,3), s_ele(3,3), f, dn_dpz(2), sh_mat(6,6), m8(8,8)
complex(rp) cvec(6), ccvec(6,1)
complex(rp) eval(6), evec(6,6), dd(2,2), gv(2,1), w_vec(6,2), q_vec(8), cmat6(6,6)

integer i, p, k, kk, ix, ix_lat_branch, ix_ele, ix_ele2
integer stat, ipiv(2), ipiv6(6)
logical err

character(200) bmad_lat, param_file

namelist / params / bmad_lat, ix_lat_branch, ix_ele, ix_ele2

!

param_file = 'spin_matching.init'
if (command_argument_count() == 1) call get_command_argument(1, param_file)

ix_lat_branch = 0
ix_ele2 = -1
open (1, file = param_file)
read (1, nml = params)
close (1)

!

call bmad_parser (bmad_lat, lat)
call lat_to_ptc_layout (lat)
branch => lat%branch(ix_lat_branch)

!

call ptc_spin_matching_calc (branch, match_info)


if (ix_ele2 == -1) ix_ele2 = ix_ele

do i = ix_ele, ix_ele2
  print *, '!============================================================'
  sm => match_info(i)
  print '(a, i6, 3x, a, 5x, a)', 'At element: ', i, trim(branch%ele(i)%name), trim(key_name(branch%ele(i)%key))
  print '(a, 6f12.6)', 'l_axis:   ', sm%axis%l
  print '(a, 6f12.6)', 'n0:       ', sm%axis%n0
  print '(a, 6f12.6)', 'm_axis:   ', sm%axis%m
  print '(a, 6f12.6)', 'alpha:    ', sm%alpha
  print '(a, 6f12.6)', 'beta:     ', sm%beta
  print '(a, 6f12.6)', 'orb0:     ', sm%orb0
  print '(a, 6f12.6)', 'dn_dpz:   ', sm%dn_dpz
  print '(a, 6f12.6)', 'm_1turn:    '
  call mat_type (sm%m_1turn)
  print '(a, 6f12.6)', 'm_ele:    '
  call mat_type (sm%m_ele)

  s_1turn = quat_to_w_mat(sm%sq_1turn)
  s_ele = quat_to_w_mat(sm%sq_ele)

  print '(a, 6f12.6)', 's_1turn:    '
  call mat_type (s_1turn)
  print '(a, 6f12.6)', 's_ele:    '
  call mat_type (s_ele)

  print '(a, 3f12.6)', 'Axis dot products:', dot_product(sm%axis%l, sm%axis%n0), dot_product(sm%axis%n0, sm%axis%m), dot_product(sm%axis%m, sm%axis%l)
  print '(a, 3f12.6)', 'n.dn:             ', dot_product(sm%axis%n0, sm%dn_dpz)
  print '(a, 3f12.6)', 'det(D) 1turn, ele:', determinant(sm%m_1turn(7:8,7:8)), determinant(sm%m_ele(7:8,7:8))
  print '(a, 3f12.6)', 'det(S) 1turn, ele:', determinant(s_1turn), determinant(s_ele)
  print '(a, 3f12.6, 5x, 3f12.6)', 'n-to-n 1turn, ele:', matmul(s_1turn, sm%axis%n0) - sm%axis%n0, matmul(s_ele, match_info(i-1)%axis%n0) - sm%axis%n0

  call mat_make_unit(m8)
  do ix = i+1, branch%n_ele_track
    m8 = matmul(match_info(ix)%m_ele, m8)
  enddo
  do ix = 1, i
    m8 = matmul(match_info(ix)%m_ele, m8)
  enddo

  print '(a, 6f12.6)', 'm_1turn:    '
  call mat_type (m8)

  !---------------------------------------------------
  !cycle

  call mat_eigen(sm%m_1turn(1:6,1:6), eval, evec, err)
  print *
  do k = 1, 6
    print '(2f10.5, 6x, 6(2f10.5, 2x))', eval(k), evec(k,:)
  enddo

  print *
  do kk = 1, 3
    k = 2*kk - 1
    f = 2 * aimag(evec(k,1) * conjg(evec(k,2)) + evec(k,3) * conjg(evec(k,4)) + evec(k,5) * conjg(evec(k,6)))
    if (f < 0) then
      evec(k:k+1,:) = evec(k+1:k:-1, :)
      eval(k:k+1) = [eval(k+1), eval(k)]
    endif

    evec(k,:)   = evec(k,:) / sqrt(abs(f))
    evec(k+1,:) = evec(k+1,:) / sqrt(abs(f))

    print '(2f10.5, 6x, 6(2f10.5, 2x))', eval(k), evec(k,:)
    print '(2f10.5, 6x, 6(2f10.5, 2x))', eval(k+1), evec(k+1,:)
  enddo

  print *
  sh_mat = 0
  sh_mat(1,2) = -1
  sh_mat(3,4) = -1
  sh_mat(5,6) = -1
  sh_mat(2,1) =  1
  sh_mat(4,3) =  1
  sh_mat(6,5) =  1
  do kk = 1, 3
    k = 2*kk - 1
    print '(3(2f10.5, 4x))', sum(conjg(evec(k,:)) * matmul(sh_mat, evec(k,:)))
  enddo

  print *
  do k = 1, 6
    dd = sm%m_1turn(7:8,7:8)
    dd(1,1) = dd(1,1) - eval(k)
    dd(2,2) = dd(2,2) - eval(k)
    gv(:,1) = -matmul(sm%m_1turn(7:8,1:6), evec(k,1:6))
    call zgesv_f95 (dd, gv, ipiv, stat)
    w_vec(k,:) = gv(:,1)
    q_vec = [evec(k,:), w_vec(k,:)]
    print '(16f10.5)', matmul(sm%m_1turn, q_vec) - eval(k) * q_vec
  enddo

  dn_dpz = 0
  do kk = 1, 3
    k = 2*kk - 1
    dn_dpz = dn_dpz - 2 * aimag(conjg(evec(k,5)) * w_vec(k,:))
  enddo

  print *
  print '(a, 6f12.6)', 'dn_dpz:', dn_dpz
  print '(a, 6f12.6)', 'dn_dpz:', dot_product(sm%axis%l, sm%dn_dpz), dot_product(sm%axis%m, sm%dn_dpz)

  cmat6 = transpose(evec)
  ccvec(:,1) = [0, 0, 0, 0, 0, 1]
  call zgesv_f95(cmat6, ccvec, ipiv6, stat)
  print *
  print '(6(2f9.5, 4x))', matmul(transpose(evec), ccvec(:,1))
  print '(6(2f9.5, 4x))', ccvec(:,1)
  print '(a, 6f12.6)', 'dn_dpz:', matmul(transpose(w_vec), ccvec(:,1))

enddo

end program
