!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine check_end (point_, ix, string_in)

  use sr_struct
  use sr_interface

  implicit none

  type (outline_pt_struct) point_(*)
  integer ix
  character*(*) string_in

!

  if (abs (point_(ix)%x - 0.045) > .0001) then
      type *, 'WARNING: ', string_in, ' ON OUTLINE NOT AT X = 4.5 CM: ', &
                                                                point_(ix)%name
  endif

  return

end subroutine
