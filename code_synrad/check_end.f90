subroutine check_end (point, ix, string_in)

  use synrad_struct
  use synrad_interface, except => check_end

  implicit none

  type (outline_pt_struct) point(:)
  integer ix
  character(*) string_in

!

  if (abs (point(ix)%x - 0.045) > .0001) then
      print *, 'WARNING: ', string_in, ' ON OUTLINE NOT AT X = 4.5 CM: ', &
                                                                point(ix)%name
  endif

  return

end subroutine
