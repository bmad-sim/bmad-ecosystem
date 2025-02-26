program compare_tracking_methods_text

use bmad

implicit none

type (lat_struct), target :: lat

character(38) :: final_str
integer, parameter :: n_methods = ubound(tracking_method_name, 1)
integer :: i, j
 
type (coord_struct) end_orb

call bmad_parser ('compare_tracking_methods_text.bmad', lat)

bmad_com%auto_bookkeeper = .false.

open (1, file = 'compare_tracking_methods_text.out')

DO i = 1, lat%n_ele_max - 1
   DO j = 1, n_methods
      if(.not. valid_tracking_method(lat%ele(i), lat%param%particle, j) .or. j == custom$) cycle
      if(lat%ele(i)%key == elseparator$ .and. (j == runge_kutta$ .or. j == time_runge_kutta$)) cycle
      lat%ele(i)%tracking_method = j
      call set_flags_for_changed_attribute(lat%ele(i),lat%ele(i)%tracking_method)
      call lattice_bookkeeper(lat)
      if (j == linear$) then
        lat%ele(i)%tracking_method = symp_lie_ptc$
        if(lat%ele(i)%key == beambeam$) lat%ele(i)%tracking_method = bmad_standard$
        call set_flags_for_changed_attribute(lat%ele(i),lat%ele(i)%tracking_method)
        call lattice_bookkeeper(lat)
        call make_mat6 (lat%ele(i), lat%param, lat%particle_start)
        lat%ele(i)%tracking_method = j
        call set_flags_for_changed_attribute(lat%ele(i),lat%ele(i)%tracking_method)
        call lattice_bookkeeper(lat)
      endif
      call init_coord (lat%particle_start, lat%particle_start, ele = lat%ele(i), element_end = upstream_end$)
      call track1 (lat%particle_start, lat%ele(i), lat%param, end_orb)
      final_str = '"' // trim(lat%ele(i)%name) // ':' // trim(tracking_method_name(j)) // '"' 
      write (1,'(a,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)',advance='no') final_str, end_orb%vec(1), end_orb%vec(2), end_orb%vec(3), end_orb%vec(4), end_orb%vec(5), end_orb%vec(6), end_orb%t
      write (1,*)
   END DO
   write (1,*)
END DO

close(1)

end program
