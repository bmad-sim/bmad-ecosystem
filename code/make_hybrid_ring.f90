!+
! Subroutine MAKE_HYBRID_RING (RING_IN, USE_ELE,
!                                           REMOVE_MARKERS, RING_OUT, IX_OUT)
!
! Subroutine to concatinate together elements to make a hybrid ring
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING_IN        -- Ring_struct: Input ring.
!     USE_ELE(n_ele_maxx) 
!                    -- Logical array: USE_ELE(I) = .true. indicates an element
!                            that is NOT to be concatenated.
!     REMOVE_MARKERS -- Logical: If .true. then marker elements in USE_ELE
!                            are removed from RING_OUT. In this case IX_OUT
!                            points to the element before the marker
!     RING_IN
!       %PARAM%Z_DECOUPLED -- Logical: If .true. then the suboutine
!                       assumes that the RF is off and the Z coordinate is
!                       decoupled from X and Y. In this case the routine can
!                       multiply 4x4 matrices instead of the full 6x6.
!     RING_OUT
!       %PARAM%SYMMETRY   -- Integer: See BMAD_STRUCT for logical values.
!
! Output:
!     RING_OUT  -- Ring_struct. Input ring.
!     IX_OUT(n_ele_maxx) 
!               -- Integer array. Ix_out(i) is the index for ring_in%ele_(i) 
!                  of the corresponding element in ring_out%ele(). 
!                  ix_out(i) set to 0 if ring_in%ele_(i) is concatenated.
!
! Note: You need to have made the ring_in%ele_()%mat6 matrices before you
!       call this routine.
!-

!$Id$
!$Log$
!Revision 1.4  2002/01/08 21:44:39  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.3  2001/11/29 19:39:53  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.2  2001/09/27 18:31:53  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine make_hybrid_ring (r_in, use_ele, remove_markers, r_out, ix_out)

  use bmad_struct
  use bmad_interface

  implicit none


  type (ring_struct)  r_in, r_out

  real e_vec(4)

  integer j_in, i_out, ix_out(:), i
  integer n_ele, j, ix, ic, o_key

  logical init_ele_needed, remove_markers, use_ele(:), out_symmetry

! Init

  i_out = 0                        ! index for current out ring
  r_out%ele_(0) = r_in%ele_(0)     !
  init_ele_needed = .true.         ! we need to init out ring element

  if (r_out%param%symmetry /= ew_antisymmetry$) then       ! normal
    n_ele = r_in%n_ele_ring
    if (n_ele == 0) then
      type *, 'ERROR IN MAKE_HYBRID_RING: RING_IN.N_ELE_RING = 0!'
      call err_exit
    endif
  else                                                 ! use only 1/2 ring
    n_ele = r_in%n_ele_symm
    if (n_ele == 0) then
      type *, 'ERROR IN MAKE_HYBRID_RING: RING_IN.N_ELE_SYMM = 0!'
      call err_exit
    endif
  endif

! loop over all in ring elements

  do j_in = 1, n_ele

! if a match...

    if (use_ele(j_in)) then

! if current out-element is a hybrid then calculate dispersion part of mat6

      if (i_out /= 0) then
        if (r_out%ele_(i_out)%key == hybrid$ .and. r_in%param%z_decoupled)  &
                       call mat6_dispersion (r_out%ele_(i_out)%mat6, e_vec)
      endif

! on to the next out-element which is a simple element

      if (remove_markers .and. r_in%ele_(j_in)%key == marker$) then
        ix_out(j_in) = i_out
      else
        i_out = i_out + 1                     ! starting next element
        r_out%ele_(i_out) = r_in%ele_(j_in)   ! single element
        ix_out(j_in) = i_out
      endif

      init_ele_needed = .true.                ! need to init next ele

! here if no match found...
! If this is the first element after a matched element then just transfer in
! to out. Else modify the out MAT6 transfer matrix using the in MAT6 matrix.

    else

      ix_out(j_in) = 0          ! point to nothing

      if (init_ele_needed) then
        i_out = i_out + 1                       ! starting next element
        r_out%ele_(i_out) = r_in%ele_(j_in)
        r_out%ele_(i_out)%control_type = free$
        r_out%ele_(i_out)%n_slave = 0
        r_out%ele_(i_out)%ix1_slave = 0
        r_out%ele_(i_out)%ix2_slave = -1
        r_out%ele_(i_out)%n_lord = 0
        r_out%ele_(i_out)%ic1_lord = 0
        r_out%ele_(i_out)%ic2_lord = -1

        if (r_in%param%z_decoupled) then
          e_vec(1:4) = r_out%ele_(i_out)%mat6(1:4, 6)
        endif
        init_ele_needed = .false.

      else
        if (r_in%ele_(j_in)%key == marker$) cycle

        if (r_in%param%z_decoupled) then
          e_vec = matmul(r_in%ele_(j_in)%mat6(1:4,1:4), e_vec)
          e_vec(1:4) = e_vec(1:4) + r_in%ele_(j_in)%mat6(1:4, 6)
          r_out%ele_(i_out)%mat6(1:4,1:4) = &
                                     matmul(r_in%ele_(j_in)%mat6(1:4,1:4), &
                                            r_out%ele_(i_out)%mat6(1:4,1:4))
        else
          r_out%ele_(i_out)%mat6 = matmul(r_in%ele_(j_in)%mat6, &
                                               r_out%ele_(i_out)%mat6)
        endif

        r_out%ele_(i_out)%s = r_in%ele_(j_in)%s
        r_out%ele_(i_out)%value(l$) = r_out%ele_(i_out)%value(l$) +  &
                                             r_in%ele_(j_in)%value(l$)
        r_out%ele_(i_out)%value(hkick$) = r_out%ele_(i_out)%value(hkick$) +  &
                                             r_in%ele_(j_in)%value(hkick$)
        r_out%ele_(i_out)%value(vkick$) = r_out%ele_(i_out)%value(vkick$) +  &
                                             r_in%ele_(j_in)%value(vkick$)
        r_out%ele_(i_out)%value(x_limit$) = r_in%ele_(j_in)%value(x_limit$)
        r_out%ele_(i_out)%value(y_limit$) = r_in%ele_(j_in)%value(y_limit$)

        o_key = r_out%ele_(i_out)%key 
        if (r_in%ele_(j_in)%key == drift$ .and. &
                          (o_key == drift$ .or. o_key == marker$)) then
          r_out%ele_(i_out)%name = 'DRIFT_HYBRID' 
          r_out%ele_(i_out)%key = drift$
        else
          r_out%ele_(i_out)%name = 'HYBRID'
          r_out%ele_(i_out)%key = hybrid$
          if (r_in%ele_(j_in)%coupled) r_out%ele_(i_out)%coupled = .true.
        endif

        r_out%ele_(i_out)%x       = r_in%ele_(j_in)%x
        r_out%ele_(i_out)%y       = r_in%ele_(j_in)%y
        r_out%ele_(i_out)%c_mat   = r_in%ele_(j_in)%c_mat
        r_out%ele_(i_out)%gamma_c = r_in%ele_(j_in)%gamma_c

      endif

    endif ! use_ele

  enddo

! end cleanup

  if (r_out%ele_(i_out)%key == hybrid$ .and. r_in%param%z_decoupled)  &
                          call mat6_dispersion (r_out%ele_(i_out)%mat6, e_vec)

  out_symmetry = r_out%param%symmetry   ! save symmetry
  r_out%parameters = r_in%parameters
  r_out%n_ele_ring = i_out
  r_out%param%symmetry = out_symmetry

  if (r_out%param%symmetry /= ew_antisymmetry$) then       ! normal
    r_out%n_ele_use = i_out
    if (r_in%n_ele_symm /= 0) then
      r_out%n_ele_symm = ix_out(r_in%n_ele_symm)
    else
      r_out%n_ele_symm = 0
    endif
  else
    r_out%n_ele_symm = i_out
    r_out%n_ele_use = i_out
  endif

! put control elements in

  do j_in = r_in%n_ele_ring+1, r_in%n_ele_max
    if (use_ele(j_in)) then
      i_out = i_out + 1
      ix_out(j_in) = i_out
      r_out%ele_(i_out) = r_in%ele_(j_in)
    endif
  enddo
  r_out%n_ele_max = i_out

! update control pointers
! if symmetry is EW_ANTISYMMETRY$ then ignore fact that controllers may
! control non-existant elements in the east.

  r_out%control_(:)%ix_lord = 0
  r_out%control_(:)%ix_slave = 0
  r_out%ic_(:) = r_in%ic_(:)

  do i_out = 1, r_out%n_ele_max

    do j = r_out%ele_(i_out)%ix1_slave, r_out%ele_(i_out)%ix2_slave
      ix = r_in%control_(j)%ix_slave
      if (ix_out(ix) == 0) then
        if (r_out%param%symmetry /= ew_antisymmetry$ .or.  &
                         ix <= n_ele .or. ix > r_in%n_ele_ring) then
          type '(a, /, 2a, /, 2a)',  &
      ' WARNING IN MAKE_HYBRID_RING: SLAVE ELEMENT LIST NOT COMPLETE.',  &
      '         FOR LORD ELEMENT: ', r_out%ele_(i_out)%name,  &
      '         CANNOT FIND: ', r_in%ele_(ix)%name
        endif
        r_out%control_(j)%ix_lord = i_out
        r_out%control_(j)%ix_slave = n_ele_maxx     ! point to dummy ele
        r_out%control_(j)%ix_attrib = -1
        r_out%control_(j)%coef = 0
      else
        r_out%control_(j)%ix_lord = i_out
        r_out%control_(j)%ix_slave = ix_out(ix)
        r_out%control_(j)%ix_attrib = r_in%control_(j)%ix_attrib
        r_out%control_(j)%coef = r_in%control_(j)%coef
      endif
    enddo

    do j = r_out%ele_(i_out)%ic1_lord, r_out%ele_(i_out)%ic2_lord
      ic = r_out%ic_(j)
      ix = r_in%control_(ic)%ix_lord
      if (ix == 0) then
        type *, 'WARNING IN MAKE_HYBRID_RING: LORD ELEMENT ',  &
                    'LIST NOT COMPLETE FOR: ', r_out%ele_(i_out)%name
      endif
    enddo

  enddo

! end

  if (r_out%n_ele_ring == 0) then
    type *, 'ERROR IN MAKE_HYBRID_RING: OUTPUT RING HAS 0 ELEMENTS!'
    call err_exit
  endif

  call check_ring_controls (r_out, .false.)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mat6_dispersion (mat6, e_vec)
!
! Subroutine to put the dispersion into ELE.MAT6 given the eta vector E_VEC
!
! Input:
!   E_VEC(4) -- Real: eta vector
!
! Output:
!   mat6(6,6) -- Real: Matrix with 4x4 x-y submatrix already made.
!-



subroutine mat6_dispersion (mat6, e_vec)

  use bmad_struct
  use bmad_interface

  implicit none

  real, intent(inout) :: mat6(6,6)
  real, intent(in) :: e_vec(*)

  real e2_vec(4)

  integer i

!

  mat6(1:4, 6) = e_vec(1:4)

  e2_vec(1) = -e_vec(2)
  e2_vec(2) =  e_vec(1)
  e2_vec(3) = -e_vec(4)
  e2_vec(4) =  e_vec(3)

  mat6(5,1:4) = matmul (e2_vec, mat6(1:4,1:4))

end subroutine
