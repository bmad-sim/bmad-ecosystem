module bmad2ac_mod

use bmad

contains

!--------------------------------------------------------------------------

subroutine bmad2ac(file_name, erl, ele_name,xmap,zmap,tmap)

implicit none

type(lat_struct), target :: erl
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)

integer i, j, n, n_loc
real(dp) x0, z0, th0, c0, s0, xloc, zloc, xmap, zmap, tmap

character(*) file_name, ele_name
character(80) file_name1, file_name2
character(40) name

!

file_name1 = file_name

call file_suffixer (file_name1, file_name1, '.lat_list', .true.)
open (1, file = file_name1)
call file_suffixer (file_name1, file_name2, '.lat_list_no_offset', .true.)
open (2, file = file_name2)

call lat_ele_locator(ele_name, erl, eles, n_loc)
if (n_loc == 0) then
  print *, 'ERROR: CANNOT FIND ELEMENT: ', ele_name
  print *, '       WILL USE BEGINNING OF LATTICE INSTEAD!'
  allocate (eles(1))
  eles(1)%ele => erl%ele(0)
endif

x0=eles(1)%ele%floor%r(1) ; z0=eles(1)%ele%floor%r(3) ; th0=eles(1)%ele%floor%theta
c0 = cos(tmap-th0) ; s0 = sin(tmap-th0)

do n = 0, ubound(erl%branch, 1)
  do i = 0, erl%branch(n)%n_ele_track

    if (i == 0 .and. n == 0) cycle
    ele => erl%branch(n)%ele(i)

    xloc = xmap + c0*(ele%floor%r(1)-x0) +  s0*(ele%floor%r(3)-z0)
    zloc = zmap - s0*(ele%floor%r(1)-x0) +  c0*(ele%floor%r(3)-z0)

    name = ele%name
    if (i == 0) name = 'BRANCH'

    write (1, '(2es18.10, 2x, a)') xloc, zloc, name
    write (2, '(i6, 2es18.10, 2x, a)') i, ele%floor%r(1), ele%floor%r(3), name
  end do
enddo

print *, 'Created: ', trim(file_name1)
print *, 'Created: ', trim(file_name2)
close (1)
close (2)

! list of elements by name

call file_suffixer (file_name1, file_name1, '.name_list', .true.)
open (1, file = file_name1)

do n = 0, ubound(erl%branch, 1)
  ele_loop: do i = 1, erl%branch(n)%n_ele_track

    ele => erl%branch(n)%ele(i)
    do j = 1, i-1
      if (erl%branch(n)%ele(j)%name == ele%name) cycle ele_loop  
    end do

    if (ele%key == sbend$) then
      write (1, '(a, 2x, a, 6es18.10)') ele%name, key_name(ele%key), ele%value(l$), ele%value(x1_limit$), &
               ele%value(angle$), ele%value(e1$), ele%value(e2$)
    else if (ele%key == wiggler$) then
      write (1, '(a, 2x, a, 6es18.10)') ele%name, key_name(ele%key), ele%value(l$), ele%value(x1_limit$)
    else if (ele%type == "BPM") then
      write (1, '(a30, a17, 11x, 2es18.10)') ele%name, "BPMON", ele%value(l$), ele%value(x1_limit$)
    else
      write (1, '(a, 2x, a, 6es18.10)') ele%name, key_name(ele%key), ele%value(l$), ele%value(x1_limit$)
    endif

    if (ele%value(l$) < 0) then
      print *, 'WARNING: Element has negative drift length!'
      print *, '         Element: ', trim(ele%name)
      print *, '         Length:', ele%value(l$)
    endif

  end do ele_loop
enddo

print *, 'Created: ', trim(file_name1)
close (1)

end subroutine bmad2ac

end module

!------------------------------------------------------------------------------------------

program bmad2autocad

use bmad2ac_mod

implicit none

type(lat_struct) :: erl

!

call bmad_parser('lat.bmad', erl)

call bmad2ac('lat.bmad', erl, 'IP_L3',&
             0.8915840379E+06*12/39.3700787D+0,0.8496528770E+06*12/39.3700787D+0,&
             0.2281071935E+02*atan(1.D+0)/45.D+0)

end program
