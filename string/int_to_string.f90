!	call int_to_string(buffer,jnum ,setup_data)
! converts jnum integers, contained in the array BUFFER into strings of
! characters. Each SETUP_DATA string is a maximum of 1024 characters long.

	subroutine int_to_string(buffer,jnum ,setup_data)
	implicit none
!
        integer*4::  buffer(*)
        integer*4::  jnum
        character(1024)::  setup_data(*)
        integer*4,save:: j, k, count
	j=1
	k=1
	count = 1
	do while (count.le.jnum )  !Loop until all data has been converted.
	 if (k .gt. 1024)then    !Once the first string is filled,
	  j = j+1                !go to next string in array
	  k = 1                  !and reset the character index
	 endif
         setup_data(j)(k:k) = char( buffer(count) )
         k = k + 1
	 count = count + 1       !Keep track of total data.
        enddo
	return
	end
