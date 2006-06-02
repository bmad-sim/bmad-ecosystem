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
!   ele  -- Ele_struct: ATTRIBUTE_INDEX will restrict the name search to 
!             valid attributes of the given element. Note: If 
!             ELE%KEY = OVERLAY then the entire name table will be searched.
!   name -- Character(40): Attribute name. Must be uppercase.
!
! Output:
!   attribute_index -- Integer: Index of the attribute. If the attribute name
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
  use bmad_interface, except => attribute_index

  implicit none

  type (ele_struct) ele, ele0

  integer i, j, k, key, num, ilen, n_abbrev, ix_abbrev
  integer attrib_num(n_key)
  integer attrib_ix(n_key, n_attrib_special_maxx)
  integer at_index

  character(*) name
  character(40) name40
  character(40) attrib_name_array(n_key, n_attrib_special_maxx)
      
  logical :: init_needed = .true.

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

  name40 = name          ! make sure we have 40 characters
  key = ele%key
  at_index = 0           ! match not found

  ilen = len_trim(name)
  if (ilen == 0) return
  n_abbrev = 0            ! number of abbreviation matches

! Overlays search all types of elements

  if (key == overlay$) then
    do k = 1, n_key
      do i = 1, attrib_num(k)
        if (attrib_name_array(k, i) == name40) then
          at_index = attrib_ix(k, i)
          return
        endif
        if (attrib_name_array(k, i)(1:ilen) == name40(1:ilen)) then
          n_abbrev = n_abbrev + 1
          ix_abbrev = attrib_ix(k, i)
        endif 
      enddo
    enddo

    if (name40 == 'CURRENT') then
      at_index = current$
      return
    endif

! else only search this type of element

  elseif (key > 0 .and. key <= n_key) then
    do i = 1, attrib_num(key)
      if (attrib_name_array(key, i) == name40) then
        at_index = attrib_ix(key, i)
        return
      endif
      if (attrib_name_array(key, i)(1:ilen) == name40(1:ilen)) then
        n_abbrev = n_abbrev + 1
        ix_abbrev = attrib_ix(key, i)
      endif 
    enddo      

    if (key == rfcavity$ .and. name40 == 'LAG') then
      at_index = phi0$
      return
    endif

! error

  else
    print *, 'ERROR IN ATTRIBUTE_INDEX: BAD KEY', key
    call err_exit
  endif

! If there is one unique abbreviation then use it.

  if (n_abbrev == 1) at_index = ix_abbrev

end function
