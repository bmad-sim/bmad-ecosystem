!+             
! Function attribute_index (ele, name)
!
! Function to return the index of a attribute for a given BMAD element type
! and the name of the attribute.
!
! Modules Needed:
!   use bmad
!
! Input:
!     ele  -- Ele_struct: ATTRIBUTE_INDEX will restrict the name search to 
!             valid attributes of the given element. Note: If 
!             ELE%KEY = OVERLAY then the entire name table will be searched.
!     name -- Character*16: Attribute name. Must be uppercase
!
! Output:
!     attribute_index -- Integer: Index of the attribute. If the attribute name
!                            is not appropriate then 0 will be returned.
!
! Example:
!     ele%key = sbend$
!     ix = attribute_index (ele, 'K1')
! Result:
!     ix -> k1$
!-

#include "CESR_platform.inc"

function attribute_index (ele, name) result (at_index)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele, ele0

  integer i, j, k, key, num
  integer attrib_num(n_key)
  integer attrib_ix(n_key, n_attrib_special_maxx)
  integer at_index

  character*(*) name
  character*16 name16
  character*16 attrib_name_array(n_key, n_attrib_special_maxx)
      
  logical init_needed / .true. /

!-----------------------------------------------------------------------
! Init: We make a short list to compare against to make things go faster

  if (init_needed) then

    do i = 1, n_key
      ele0%key = i
      num = 0
      do j = 1, n_attrib_special_maxx
        if (attribute_name(ele0, j) /= null_name) then
          num = num + 1
          attrib_name_array(i, num) = attribute_name(ele0, j)
          attrib_ix(i, num) = j
        endif
      enddo
      attrib_num(i) = num
    enddo
    init_needed = .false.

  endif

!-----------------------------------------------------------------------
! search for name

  name16 = name          ! make sure we have 16 characters
  key = ele%key
  at_index = 0           ! match not found

  if (key == overlay$) then
    do k = 1, n_key
      do i = 1, attrib_num(k)
        if (attrib_name_array(k, i) == name16) then
          at_index = attrib_ix(k, i)
          return
        endif
      enddo
    enddo

  elseif (key > 0 .and. key <= n_key) then
    do i = 1, attrib_num(key)
      if (attrib_name_array(key, i) == name16) then
        at_index = attrib_ix(key, i)
        return
      endif
    enddo      

    if (key == rfcavity$ .and. name16 == 'LAG') at_index = phi0$

  else
    print *, 'ERROR IN ATTRIBUTE_INDEX: BAD KEY', key
    call err_exit
  endif

end function
