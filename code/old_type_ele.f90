!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_twiss, 
!                                                           type_control, ring)
!
! Subroutine to type out information on an element. 
! See also the subroutine type2_ele.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ELE              -- Ele_struct: Element
!   TYPE_ZERO_ATTRIB -- Logical: If true then type all attributes even if the
!                          attribute value is 0.
!   TYPE_MAT6        -- Integer:
!                          TYPE_MAT6 = 0   => Do not type ELE%MAT6
!                          TYPE_MAT6 = 4   => Type 4X4 XY submatrix
!                          TYPE_MAT6 = 6   => Type full 6x6 matrix
!   TYPE_TWISS       -- Logical: If true then type the twiss parameters
!                          at the end of the element.
!   TYPE_CONTROL     -- Logical: If true then type control status.
!   RING             -- Ring_struct: Needed for control typeout.
!
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine type_ele (ele, type_zero_attrib, type_mat6, type_twiss,  &
                                                       type_control, ring)

  use bmad_struct
  implicit none
  type (ele_struct)  ele
  type (ring_struct)  ring

  integer i, j, n, type_mat6, ix, iv, ic, ct
  real coef
  character*16 a_name, attribute_name, name, null_type / ' ' /

  logical type_twiss, type_control, type_zero_attrib

!

  ct = ele%control_type

  type *, 'Element Name: ', ele%name
  if (ele%type /= null_type) type *, 'Element Type: "', ele%type, '"'
  if (ele%key <= 0) then
    type *, 'Key: UNKNOWN!', ele%key
  else
    type *, 'Key: ', key_name(ele%key)
    if (.not. type_zero_attrib)  &
            type *, 'Attribute values [Only non-zero values shown]:'
    if (ct == overlay_lord$) then
      i = ele%ix_value
      ix = ring%control_(ele%ix1_slave)%ix_slave
      name = attribute_name(ring%ele_(ix)%key, i)
      type '(5x, 2a, 1pe12.3)', name, ' = ', ele%value(i)
    else
      do i = 1, attrib_maxx
        if (attribute_name(ele%key, i) /= null_name) then
          if (ele%value(i) /= 0 .or. type_zero_attrib)  &
                         type '(5x, 2a, 1pe12.3)',  &
                         attribute_name(ele%key, i), ' = ', ele%value(i)
        endif
      enddo
    endif
  endif
  if (ele%is_off) type *, '*** Note: Element is turned OFF ***'
  type '(1x, a, f13.3)', 'S: ', ele%s


  if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
    type *
    type *, 'ERROR: N_SLAVE DOES NOT MATCH IX2_SLAVE-IX1_SLAVE+1 !!!'
    type *, ele%n_slave, ele%ix2_slave, ele%ix1_slave
  endif

  if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
    type *
    type *, 'ERROR: N_LORD DOES NOT MATCH Ic1_LORD-Ic2_LORD+1 !!!'
    type *, ele%n_lord, ele%ic1_lord, ele%ic2_lord
  endif

  if (type_control) then
    type *
    if (ct <= 0) then
      type *, 'Control_type: UNKNOWN!', ct
    else
      type *, 'Control_type: ', control_name(ct)

      if (ele%n_slave /= 0) then
        type '(1x, a, i4)', 'Slaves: Number:', ele%n_slave
        type *, '    Name            Ring_index  Attribute        Coefficient'
        do i = ele%ix1_slave, ele%ix2_slave
          j = ring%control_(i)%ix_slave
          coef = ring%control_(i)%coef
          if (ct == super_lord$ .or. &
                          ring%ele_(j)%control_type == scale_slave$) then
            a_name = '   ----'
          else
            iv = ring%control_(i)%ix_attrib
            a_name = attribute_name(ring%ele_(j)%key, iv)
          endif
          type '(5x, a, i10, 2x, a16, 1pe12.3)', ring%ele_(j)%name, j,  &
                                                            a_name, coef
        enddo
      endif

      if (ele%n_lord /= 0) then
        type '(1x, a, i4)', 'Lords: Number:', ele%n_lord
        type *, '    Name            Ring_index  Attribute        Coefficient'
        do i = ele%ic1_lord, ele%ic2_lord
          ic = ring%ic_(i)
          j = ring%control_(ic)%ix_lord
          coef = ring%control_(ic)%coef
          ct = ele%control_type
          if (ct == super_slave$ .or. ct == scale_slave$) then
            a_name = '   ----'
          else
            iv = ring%control_(ic)%ix_attrib
            a_name = attribute_name(ele%key, iv)
          endif
          type '(5x, a, i10, 2x, a16, 1pe12.3)', ring%ele_(j)%name,  &
                                                        j, a_name, coef
        enddo
      endif

    endif
  endif

! twiss output

  if (type_twiss) then
    type *
    type '(10x, a)',  &
            'Beta     Alpha     Gamma       Phi        Eta       Etap'
    type '(a, f11.4, 2f10.3, f10.4, 2f11.4)', ' X:', ele%x%beta,  &
                 ele%x%alpha, ele%x%gamma, ele%x%phi, ele%x%eta, ele%x%etap
    type '(a, f11.4, 2f10.3, f10.4, 2f11.4)', ' Y:', ele%y%beta,  &
                 ele%y%alpha, ele%y%gamma, ele%y%phi, ele%y%eta, ele%y%etap
  endif

! MAT6 output

  n = type_mat6

  if (n /= 0) type *

  if (any(abs(ele%mat6(1:n,1:n)) >= 1000)) then
    do i = 1, n
      type '(1p6e11.3)', (ele%mat6(i, j), j = 1, n)
    enddo
  else
    do i = 1, n
      type '(6f10.5)', (ele%mat6(i, j), j = 1, n)
    enddo
  endif

  return
  end

