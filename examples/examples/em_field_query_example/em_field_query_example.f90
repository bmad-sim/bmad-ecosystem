!+ 
! Program em_field_query_example
!
! Example program to calculate the E and B fields at any point in the lattice
!
! Input (command line)
!   bmad lattice file 
!   coordinate file.dat
!		Format:
!		n (integer)
!		t_1 x_1 y_1 s_1 (reals)
!		t_2 x_2 y_2 s_2
!	    . . .
!		t_n x_n y_n s_n
!
! Output
!   coordinate file.out
!	Format:
!	n
!	t_1 x_1 y_1 s_1 Ex Ey Es Bx By Bs
!   . . .
!-

program em_field_query_example

use bmad

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele 
type (coord_struct), allocatable, target:: coord(:)
type (coord_struct), pointer :: orbit
type (em_field_struct) field
real(rp) rawvec(4)
real(rp) :: s_rel, t_rel
character(100) lat_filename
character(100) coords_filename
character(100) outfile_filename
integer, parameter ::  coordsfile = 1, outfile = 2
integer :: i, n_pts, ix_ele
logical :: local_ref_frame


!------------------------------------------
!Get lattice file name
call getarg(1, lat_filename)
call getarg(2, coords_filename)
print *, "Using lattice: ", lat_filename
print *, "Using coordinates in: ", coords_filename

!Prepare output file name
call file_suffixer (coords_filename, outfile_filename, '.out', .true.)

!Parse Lattice
call bmad_parser (lat_filename, lat)

!Open coordinate file and read the number of points
open(coordsfile, file = coords_filename)
read(coordsfile, '(i8)') n_pts

!Allocate coord(n_pts)
if (.not. allocated (coord)) then
  allocate(coord(1:n_pts))
endif

!Set coord(:) to contain these points
do i = 1, n_pts
	read (coordsfile, *) rawvec
	ix_ele =  element_at_s (lat, rawvec(4), .false.)
    call init_coord(coord(i), [rawvec(2),  0.0_rp,  rawvec(3),  0.0_rp,  0.0_rp,  0.0_rp], &
                    lat%ele(ix_ele), inside$, particle=lat%param%particle)
	coord(i)%t = rawvec(1)
	coord(i)%s = rawvec(4)
end do
close(coordsfile)


open(outfile, file = outfile_filename)
write (outfile, '(i8)') n_pts
do i = 1, n_pts
	orbit => coord(i)
	ix_ele =  element_at_s (lat, orbit%s, .false.)
	ele => lat%ele(ix_ele)
	s_rel = orbit%s - (ele%s - ele%value(L$) ) 
	t_rel = orbit%t - (ele%ref_time - ele%value(delta_ref_time$) ) 
	! calculate the field. 
	! Note that only orbit%vec(1) = x and orbit%vec(3) = y are used in em_field_calc,
	!	and they coincide in both coordinate systems,
	!	so we can use the 'normal' routine:
	call em_field_calc (ele, lat%param, s_rel, orbit, local_ref_frame, field, .false., rf_time = t_rel)
	write (outfile, '(10es18.10)')	orbit%t, orbit%vec(1), orbit%vec(3), orbit%s, &
									field%E, field%B

end do
close(outfile)
print *, "Written: ", outfile_filename

end program
