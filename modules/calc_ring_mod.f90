module calc_ring_mod

contains

subroutine calc_ring(ring,dims,co,err,mat_err)
  use bmad
  implicit none
  type(lat_struct) ring
  integer dims
  type(coord_struct), allocatable :: co(:)
  logical err
  logical, optional :: mat_err

  integer status
  integer i

  err = .false.
  if(present(mat_err)) mat_err = .false.
  call closed_orbit_calc(ring,co,dims,err_flag=err)
  if(present(mat_err)) then
    if(err) mat_err = .true.
  endif

  if(.not. err) then
    call lat_make_mat6(ring, -1, co, err_flag=err)
    if(present(mat_err)) then
      if(err) mat_err = .true.
    endif
    if(.not. err) then
      call twiss_at_start(ring,status)
      if(status .eq. ok$) then
        err = .false.
        call twiss_propagate_all(ring, err_flag=err)
      else
        err = .true.
      endif
    endif
  endif
end subroutine

end module
