!+
! Function find_location(arr, value) result (ix_match)
!
! Routine to find the first location in array arr(:) with given value.
!
! Input:
!   arr(:)      -- real(rp), logical, or integer
!   value       -- real(rp), logical, or integer. Must match type of arr(:).
!
! Output:
!   ix_match    -- integer: Index of match. Zero if no match found.
!-

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function find_location_real(arr, value) result (ix_match)

use precision_def

implicit none

real(rp) arr(:), value
integer ix_match

!

do ix_match = 1, size(arr)
  if (arr(ix_match) == value) return
enddo

ix_match = 0

end function find_location_real

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function find_location_int(arr, value) result (ix_match)

implicit none

integer arr(:), value
integer ix_match

!

do ix_match = 1, size(arr)
  if (arr(ix_match) == value) return
enddo

ix_match = 0

end function find_location_int

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function find_location_str(arr, value) result (ix_match)

implicit none

character(*) arr(:), value
integer ix_match

!

do ix_match = 1, size(arr)
  if (arr(ix_match) == value) return
enddo

ix_match = 0

end function find_location_str

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function find_location_logic(arr, value) result (ix_match)

implicit none

logical arr(:), value
integer ix_match

!

do ix_match = 1, size(arr)
  if (arr(ix_match) .eqv. value) return
enddo

ix_match = 0

end function find_location_logic

