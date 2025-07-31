program multipass_test

use bmad

implicit none

type (lat_struct), target :: lat, lat2
type (ele_struct), pointer :: ele, ele2, lord
type (ele_pointer_struct), allocatable :: eles(:)
type (rf_stair_step_struct), pointer :: step

integer ii, ie, is, n_loc

!

open (1, file = 'output.now')

! Space_charge_method test

call bmad_parser ('multipass_and_superimpose.bmad', lat)
call write_bmad_lattice_file ('lat_out.bmad', lat)
call bmad_parser ('lat_out.bmad', lat2)

do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  write (1, '(a, i0, 3a, 2x, a)') '"scm-', ie, trim(ele%name), '" STR  ', &
      quote(space_charge_method_name(ele%space_charge_method)), &
      quote(space_charge_method_name(lat2%ele(ie)%space_charge_method))
  write (1, '(a, i0, a, es22.14)') '"ref_time-', ie, '" ABS 1E-16', ele%ref_time
enddo

call lat_ele_locator('CAVITY1', lat, eles, n_loc)
lord => eles(1)%ele


do is = 0, 2
  if (is == 0) then
    ele => lord
  else
    ele => pointer_to_slave(lord, is)
  endif

  do ii = 0, ubound(ele%rf%steps, 1)
    step => ele%rf%steps(ii)
    write (1, '(5a, 6es22.14)') '"', trim(ele%name), '-a-', int_str(ii), '" REL 1E-8', step%E_tot0, &
                                                                   step%E_tot1, step%p0c, step%dp0c, step%dE_amp
    write (1, '(5a, 6es22.14)') '"', trim(ele%name), '-b-', int_str(ii), '" REL 1E-8', step%scale, step%time, step%dt_rf, step%s
  enddo
enddo

call lat_ele_locator('END', lat, eles, n_loc)
ele => eles(1)%ele
write (1, '(a, 2es22.14)') '"END-energy" REL 4E-9', ele%value(p0c$), ele%value(E_tot$)



! Forking with a branch element

call bmad_parser ('branch_fork.bmad', lat)
call write_bmad_lattice_file ('lat_out1.bmad', lat)
call bmad_parser ('lat1.bmad', lat)

ele => lat%branch(1)%ele(2)
write (1, '(3a)')       '"BF-01"  STR  "', trim(ele%name), '"'
write (1, '(a, 2i3)')   '"BF-02"  ABS   0 ', nint(ele%value(ix_to_branch$)), nint(ele%value(ix_to_element$))
write (1, '(a, f10.4)') '"BF-03"  REL  1e-12', ele%value(p0c$)

! Multipass and superimpose

call bmad_parser ('multipass_and_superimpose.bmad', lat)

write (1, '(3a)') '"MS-01"  STR  "', trim(lat%ele(01)%name), '"'
write (1, '(3a)') '"MS-02"  STR  "', trim(lat%ele(02)%name), '"'
write (1, '(3a)') '"MS-03"  STR  "', trim(lat%ele(03)%name), '"'
write (1, '(3a)') '"MS-04"  STR  "', trim(lat%ele(04)%name), '"'
write (1, '(3a)') '"MS-05"  STR  "', trim(lat%ele(05)%name), '"'
write (1, '(3a)') '"MS-06"  STR  "', trim(lat%ele(06)%name), '"'
write (1, '(3a)') '"MS-07"  STR  "', trim(lat%ele(07)%name), '"'
write (1, '(3a)') '"MS-08"  STR  "', trim(lat%ele(08)%name), '"'
write (1, '(3a)') '"MS-09"  STR  "', trim(lat%ele(09)%name), '"'
write (1, '(3a)') '"MS-10"  STR  "', trim(lat%ele(10)%name), '"'
write (1, '(3a)') '"MS-11"  STR  "', trim(lat%ele(11)%name), '"'
write (1, '(3a)') '"MS-12"  STR  "', trim(lat%ele(12)%name), '"'
write (1, '(3a)') '"MS-13"  STR  "', trim(lat%ele(13)%name), '"'

! fiducial and flexible patch

call bmad_parser ('patch.bmad', lat)
call write_bmad_lattice_file ('lat_out2.bmad', lat)
call bmad_parser ('lat2.bmad', lat)

write (1, '(a, f12.6)')  '"P-0S" ABS 0', lat%branch(1)%ele(0)%s
write (1, '(a, es14.6)') '"P-0T" ABS 0', lat%branch(1)%ele(0)%ref_time
ele => lat%branch(1)%ele(6)
write (1, '(a, 3f12.6)') '"P-6MI" ABS 0', ele%value(z_offset$), ele%value(x_pitch$), ele%value(y_pitch$)
ele => lat%branch(1)%ele(7)
write (1, '(a, 3f12.6)') '"P-7FR" ABS 0', ele%floor%r
write (1, '(a, 3f12.6)') '"P-7FA" ABS 0', ele%floor%theta, ele%floor%phi, ele%floor%psi

! And close

close (1)

end program
