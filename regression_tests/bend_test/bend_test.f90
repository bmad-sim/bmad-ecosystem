program bend_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) ele0

real(rp) dparam(5), pa(4), f
integer ie, ip, iz, ix_param(5), np

!

open (1, file = 'output.now', recl = 200)

bmad_com%auto_bookkeeper = .false.
call bmad_parser('bend_test.bmad', lat)

f = lat%ele(0)%value(p0c$) / c_light
dparam = [0.01_rp, 0.1_rp, 0.1_rp * f, 0.1_rp, 0.1_rp]
ix_param = [g$, angle$, b_field$, l$, l_rectangle$]


do ie = 1, lat%n_ele_track
  ele0 = lat%ele(ie)

  do ip = 1, size(ix_param)
    ele => lat%ele(ie)
    ele = ele0
    np = ix_param(ip)
    ele%value(np) = ele%value(np) + dparam(ip)
    call set_flags_for_changed_attribute(ele, ele%value(np))
    call lattice_bookkeeper(lat)

    pa(1:3) = [ele%value(ix_param(1)), ele%value(ix_param(2)), ele%value(ix_param(3)) / f] 
    write (1, '(a, i0, 3a, 4f16.10)') '"', ie, '-p1-', trim(attribute_name(ele, np)), '" ABS 1E-10', pa(1:3)

    pa = [ele%value(ix_param(4)), ele%value(ix_param(5)), ele%value(e1$), ele%value(e2$)] 
    write (1, '(a, i0, 3a, 4f16.10)') '"', ie, '-p2-', trim(attribute_name(ele, np)), '" ABS 1E-10', pa
  enddo
enddo

close (1)

end program
