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
type (tao_d1_data_struct), pointer :: d1
type (ele_struct), pointer :: ele

real(rp) cbar(2,2), f
integer i, j, k, m
character(20) :: r_name = 'tao_load_data_array'
logical found

! loop over all universes 

do i = 1, size(s%u)

  u => s%u(i)

! loop over all d2_data structs

! do nothing if no d2_data
  if (.not. associated (u%d2_data)) return

  do j = 1, size(u%d2_data)
  
    call tao_hook_load_data_array (s, found)
    if (found) cycle

! loop over all data 

    do k = 1, size(u%d2_data(j)%d1)
      d1 => u%d2_data(j)%d1(k)

      select case (u%d2_data(j)%class)

      case ('orbit')
        if (d1%sub_class == 'x') then
          where (d1%d%exists) d1%d(:)%model_value = u%model_orb(d1%d(:)%ix_ele)%vec(1)
        elseif (d1%sub_class == 'y') then
          where (d1%d%exists) d1%d(:)%model_value = u%model_orb(d1%d(:)%ix_ele)%vec(3)
        else
          call sub_class_error
        endif

      case ('phase')
        if (d1%sub_class == 'x') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%x%phi
        elseif (d1%sub_class == 'y') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%y%phi
        else
          call sub_class_error
        endif

      case ('beta')
        if (d1%sub_class == 'x') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%x%beta
        elseif (d1%sub_class == 'y') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%y%beta
        else
          call sub_class_error
        endif

      case ('coupling')
        do m = lbound(d1%d, 1), ubound(d1%d, 1)
          if (.not. d1%d(m)%exists) cycle
          ele => u%model%ele_(d1%d(m)%ix_ele)
          call c_to_cbar (ele, cbar)
          if (d1%sub_class == '11') then
            f = sqrt(ele%x%beta/ele%y%beta) / ele%gamma_c
            d1%d(m)%model_value = cbar(1,1) * f
            d1%d(m)%conversion_factor = 1 / f
          elseif (d1%sub_class == '12a') then
            f = sqrt(ele%y%beta/ele%x%beta) / ele%gamma_c
            d1%d(m)%model_value = cbar(1,2) * f
            d1%d(m)%conversion_factor = 1 / f
          elseif (d1%sub_class == '12b') then
            f = sqrt(ele%x%beta/ele%y%beta) / ele%gamma_c
            d1%d(m)%model_value = cbar(1,2) * f
            d1%d(m)%conversion_factor = 1 / f
          elseif (d1%sub_class == '22') then
            f = sqrt(ele%y%beta/ele%x%beta) / ele%gamma_c
            d1%d(m)%model_value = cbar(2,2) * f
            d1%d(m)%conversion_factor = 1 / f
          else
            call sub_class_error
          endif
        enddo

      case ('cbar')
        do m = lbound(d1%d, 1), ubound(d1%d, 1)
          if (.not. d1%d(m)%exists) cycle
          ele => u%model%ele_(d1%d(m)%ix_ele)
          call c_to_cbar (ele, cbar)
          if (d1%sub_class == '11') then
            d1%d(m)%model_value = cbar(1,1)
          elseif (d1%sub_class == '12') then
            d1%d(m)%model_value = cbar(1,2)
          elseif (d1%sub_class == '21') then
            d1%d(m)%model_value = cbar(1,2)
          elseif (d1%sub_class == '22') then
            d1%d(m)%model_value = cbar(2,2)
          else
            call sub_class_error
          endif
        enddo

      case ('eta')
        if (d1%sub_class == 'x') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%x%beta
        elseif (d1%sub_class == 'y') then
          where (d1%d%exists) d1%d(:)%model_value = u%model%ele_(d1%d(:)%ix_ele)%y%beta
        else
          call sub_class_error
        endif

      case default
        call out_io (s_error$, r_name, &
                                  'UNKNOWN DATA TYPE: ' // u%d2_data(j)%class)
        return

      end select

    enddo

  enddo

enddo

!----------------------------------------------------------------------------
contains

subroutine sub_class_error

  call out_io (s_error$, r_name, 'BAD DATA SUB_TYPE: ' // d1%sub_class, &
                                      'FOR DATA TYPE: ' // u%d2_data(j)%class)

end subroutine

end subroutine

