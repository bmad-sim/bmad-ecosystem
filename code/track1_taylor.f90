!+
! Subroutine track1_taylor (start, ele, param, end)
!
! Subroutine to track through an element using the elements taylor series.
!
! Moudules needed:
!   use bmad
!
! Input:
!   start      -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- Param_struct: Beam parameters.
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
! Output:
!   end        -- Coord_struct: Ending coords.
!-

subroutine track1_taylor (start, ele, param, end)

  use ptc_interface_mod

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (param_struct), intent(inout) :: param
  type (ele_struct), intent(inout) :: ele

!

  if (.not. associated(ele%taylor(1)%term)) then
    if (bmad_status%type_out) then
      print *, 'WARNING FROM TRACK1_TAYLOR: TAYLOR SERIES NOT PRESENT FOR: ', &
                                                                      ele%name
      print *, '        I WILL MAKE A TAYLOR SERIES AROUND THE GIVEN ORBIT...'
    endif
    call ele_to_taylor(ele, param, start)
  endif

  call track_taylor (start%vec, ele%taylor, end%vec)

end subroutine
