!+
! Subroutine tao_load_data_array (s)
!
! Routine to take data from the model lattice and model orbit
! and put that into the s%u(:)%data(:) arrays.
!
! Input:
!   s  -- Super_universe_struct:
!-

subroutine tao_load_data_array(s)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_universe_struct), pointer :: u
type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: data

real(rp) cbar(2,2), f
integer i, j, k, m, ix
character(20) :: r_name = 'tao_load_data_array'
character(16) prefix, postfix
logical found

! loop over all universes 

do i = 1, size(s%u)

  u => s%u(i)

! loop over all data structs

  do j = 1, size(u%data)
  
    data => u%data(j)
    if (.not. data%exists) cycle
    call tao_hook_load_data_array (s, data, found)
    if (found) cycle
    ele => u%model%ele_(data%ix_ele)

! loop over all data 

    ix = index(data%class, ':')
    if (ix == 0) then
      prefix = data%class
      postfix = ' '
    else
      prefix = data%class(:ix-1)
      postfix = data%class(ix+1:)
    endif

    select case (prefix)

    case ('orbit')
      if (postfix == 'x') then
        data%model_value = u%model_orb(data%ix_ele)%vec(1)
      elseif (postfix == 'y') then
        data%model_value = u%model_orb(data%ix_ele)%vec(3)
      else
        call name_error
      endif

    case ('phase')
      if (postfix == 'x') then
        data%model_value = ele%x%phi
      elseif (postfix == 'y') then
        data%model_value = ele%y%phi
      else
        call name_error
      endif

    case ('beta')
      if (postfix == 'x') then
        data%model_value = ele%x%beta
      elseif (postfix == 'y') then
        data%model_value = ele%y%beta
      else
        call name_error
      endif

    case ('coupling')
      call c_to_cbar (ele, cbar)
      if (postfix == '11') then
        f = sqrt(ele%x%beta/ele%y%beta) / ele%gamma_c
        data%model_value = cbar(1,1) * f
        data%conversion_factor = 1 / f
      elseif (postfix == '12a') then
        f = sqrt(ele%y%beta/ele%x%beta) / ele%gamma_c
        data%model_value = cbar(1,2) * f
        data%conversion_factor = 1 / f
      elseif (postfix == '12b') then
        f = sqrt(ele%x%beta/ele%y%beta) / ele%gamma_c
        data%model_value = cbar(1,2) * f
        data%conversion_factor = 1 / f
      elseif (postfix == '22') then
        f = sqrt(ele%y%beta/ele%x%beta) / ele%gamma_c
        data%model_value = cbar(2,2) * f
        data%conversion_factor = 1 / f
      else
        call name_error
      endif

    case ('cbar')
      call c_to_cbar (ele, cbar)
      if (postfix == '11') then
        data%model_value = cbar(1,1)
      elseif (postfix == '12') then
        data%model_value = cbar(1,2)
      elseif (postfix == '21') then
        data%model_value = cbar(1,2)
      elseif (postfix == '22') then
        data%model_value = cbar(2,2)
      else
        call name_error
      endif

    case ('eta')
      if (postfix == 'x') then
        data%model_value = ele%x%beta
      elseif (postfix == 'y') then
        data%model_value = ele%y%beta
      else
        call name_error
      endif

    case default
      call out_io (s_error$, r_name, 'UNKNOWN DATA TYPE: ' // data%class)
      return

    end select

  enddo

enddo

!----------------------------------------------------------------------------
contains

subroutine name_error

  call out_io (s_error$, r_name, 'BAD DATA TYPE: ' // u%data%class)

end subroutine

end subroutine

