!+
! Function logic_str(logic) result (str)
!
! Routine to return a string representation (T/F) of a logical.
!
! Input:
!   logic    -- logical
!
! Output:
!   str      -- character(1), String representation.


function logic_str(logic) result (str)

implicit none

logical logic
character(1) :: str

select case (logic)
case (.true.);  str = 'T'
case (.false.); str = 'F'
end select

end function logic_str

