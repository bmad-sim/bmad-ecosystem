program bend_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) ele0

real(rp) dparam(3), vals(4), f, dg, dangle
integer ie, ip, iz, ix_param(3), np, nargs
character(40) str
character(200) lat_file

namelist / params / dg, dangle

!

lat_file = 'bend_test.bmad'
nargs = command_argument_count()
if (nargs > 0) then
  call get_command_argument(1, lat_file)
endif

!

open (1, file = 'output.now', recl = 200)

bmad_com%auto_bookkeeper = .false.
call bmad_parser(lat_file, lat)

open (2, file = lat_file)
read (2, nml = params)
close (2)

!

f = lat%ele(0)%value(p0c$) / c_light
dparam = [dg, dangle, dg * f]
ix_param = [g$, angle$, b_field$]


do ie = 1, lat%n_ele_track
  ele0 = lat%ele(ie)
  if (ie > 1) write (1, *)

  do ip = 1, size(ix_param)
    ele => lat%ele(ie)
    ele = ele0
    np = ix_param(ip)
    ele%value(np) = ele%value(np) + dparam(ip)
    call set_flags_for_changed_attribute(ele, ele%value(np))
    call lattice_bookkeeper(lat)
    str = '"' // int_str(ie) // '-' // trim(ele%name) // '-' // trim(attribute_name(ele, np))

    vals(1:3) = [ele%value(g$), ele%value(angle$), ele%value(b_field$) / f] 
    write (1, '(2a, 4f16.10)') trim(str), '  g-ang-field" ABS 1E-10', vals(1:3)

    vals = [ele%value(l$), ele%value(l_rectangle$), ele%value(e1$), ele%value(e2$)] 
    write (1, '(2a, 4f16.10)') trim(str), '  l-lr-e1/2"   ABS 1E-10', vals
  enddo
enddo

close (1)

end program
