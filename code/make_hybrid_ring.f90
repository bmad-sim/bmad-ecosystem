!+
! Subroutine make_hybrid_ring (ring_in, keep_ele,
!                    remove_markers, ring_out, ix_out, use_taylor, orb0_)
!
! Subroutine to concatinate together the elements in a ring to make
! a ring with fewer elements. This is used to speed up computation times.
! The concatinated elements in the new ring are known as hybrid elements.
!
! Note: For hybrid elements ring_out%ele_(i)%tracking_method and 
! ring_out%ele_(i)%mat6_calc_method are set as follows:
!
!   use_taylor    tracking_method     mat6_calc_method
!   ----------    ---------------     ----------------  
!   False         linear$             none$
!   True          taylor$             taylor$
!
! Note: For use_taylor = .false. You need to have made the 
! ring_in%ele_()%mat6 matrices before you call this routine.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring_in        -- Ring_struct: Input ring.
!   keep_ele(n_ele_maxx) 
!                  -- Logical array: keep_ele(I) = True indicates an element
!                        that is NOT to be concatenated. That is, there will be
!                        a corresponding element in ring_out.
!   remove_markers -- Logical: If .true. then marker elements in keep_ele
!                        are removed from RING_OUT. In this case IX_OUT
!                        points to the element before the marker
!   ring_out       -- Ring_struct: Ring with hybrid elements.
!     %param%symmetry -- Integer: See bmad_struct for logical values.
!   use_taylor     -- Logical, optional: If present and True then the
!                        hybrid elements will have a taylor series 
!                        instead of a simple linear matrix. If an element to
!                        be concatenated has a taylor series then this taylor
!                        series will be concatenated with the other elements
!                        in the hybrid element. 
!   orb0_(0:n_ele_maxx) 
!                  -- Coord_struct, optional: Central orbit for taylor stuff.
!
! Output:
!   ring_out  -- Ring_struct. Ring with hybrid elements.
!   ix_out(n_ele_maxx) 
!             -- Integer array. Ix_out(i) is the index for ring_in%ele_(i) 
!                of the corresponding element in ring_out%ele(). 
!                ix_out(i) set to 0 if ring_in%ele_(i) is concatenated.
!-

!$Id$
!$Log$
!Revision 1.8  2002/08/05 20:04:16  dcs
!Updated Documentation.
!
!Revision 1.7  2002/07/16 20:44:01  dcs
!*** empty log message ***
!
!Revision 1.6  2002/06/13 14:54:26  dcs
!Interfaced with FPP/PTC
!
!Revision 1.5  2002/02/23 20:32:18  dcs
!Double/Single Real toggle added
!
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


subroutine make_hybrid_ring (r_in, keep_ele, remove_markers, &
                                       r_out, ix_out, use_taylor, orb0_)

  use accelerator

  implicit none

  type (ring_struct), target :: r_in, r_out
  type (coord_struct), optional, volatile :: orb0_(0:)
  type (coord_struct) c0, c1
  type (ele_struct), pointer :: ele_in, ele_out
  type (real_8) y8(6)

  real(rdef) e_vec(4)

  integer j_in, i_out, ix_out(:), i
  integer n_ele, j, ix, ic, o_key

  logical init_hybrid_needed, remove_markers, keep_ele(:), out_symmetry
  logical z_decoupled, do_taylor
  logical, optional :: use_taylor

! Init

  if (present(use_taylor)) then
    do_taylor = use_taylor
  else
    do_taylor = .false.
  endif

  if (all(r_in%ele_(1:r_in%n_ele_ring)%mat6(6,5) == 0) .and. &
                                                  .not. do_taylor) then
    z_decoupled = .true.
  else
    z_decoupled = .false.
  endif

  i_out = 0                        ! index for current out ring
  r_out%ele_(0) = r_in%ele_(0)     !
  init_hybrid_needed = .true.         ! we need to init out ring element

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

    ele_in => r_in%ele_(j_in)

! if a match...

    if (keep_ele(j_in)) then

! if current out-element is a hybrid then calculate dispersion part of mat6

      if (i_out /= 0) then
        if (ele_out%key == hybrid$ .and. z_decoupled) &
                        call mat6_dispersion (ele_out%mat6, e_vec)
      endif

! on to the next out-element which is a simple element

      if (remove_markers .and. ele_in%key == marker$) then
        ix_out(j_in) = i_out
      else
        i_out = i_out + 1                     ! starting next element
        ele_out => r_out%ele_(i_out)
        ele_out = ele_in   ! single element
        ix_out(j_in) = i_out
      endif

      init_hybrid_needed = .true.                ! need to init next ele

! here if no match found...
! If this is the first element after a matched element then just transfer in
! to out. Else modify the out MAT6 transfer matrix using the in MAT6 matrix.

    else

      ix_out(j_in) = 0          ! point to nothing

      if (init_hybrid_needed) then
        i_out = i_out + 1                       ! starting next element
        ele_out => r_out%ele_(i_out)
        ele_out = ele_in
        ele_out%control_type = free$
        ele_out%n_slave = 0
        ele_out%ix1_slave = 0
        ele_out%ix2_slave = -1
        ele_out%n_lord = 0
        ele_out%ic1_lord = 0
        ele_out%ic2_lord = -1
        ele_out%tracking_method = linear$
        ele_out%mat6_calc_method = none$

        if (present (orb0_)) then
          c0 = orb0_(j_in)
        else
          c0%vec = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        endif

        if (z_decoupled) then
          e_vec(1:4) = ele_out%mat6(1:4, 6)
        elseif (do_taylor) then
          ele_out%tracking_method = taylor$
          ele_out%mat6_calc_method = taylor$
          if (.not. associated(ele_out%taylor(1)%term)) then ! construct taylor
            call ele_to_taylor (ele_out, c0, r_in%param)
          endif
        endif

        init_hybrid_needed = .false.

      else
        if (ele_in%key == marker$) cycle

        if (z_decoupled) then
          e_vec = matmul(ele_in%mat6(1:4,1:4), e_vec)
          e_vec(1:4) = e_vec(1:4) + ele_in%mat6(1:4, 6)
          ele_out%mat6(1:4,1:4) = matmul(ele_in%mat6(1:4,1:4), &
                                                      ele_out%mat6(1:4,1:4))
        elseif (do_taylor) then
          if (associated(ele_in%taylor(1)%term)) then
            call concat_taylor (ele_out%taylor, ele_in%taylor, ele_out%taylor)
          else
            call taylor_propagate1 (ele_out%taylor, ele_in, r_in%param)
          endif
        else
          ele_out%mat6 = matmul(ele_in%mat6, ele_out%mat6)
        endif

        ele_out%s = ele_in%s
        ele_out%value(l$) = ele_out%value(l$) + ele_in%value(l$)
        ele_out%value(hkick$) = ele_out%value(hkick$) + ele_in%value(hkick$)
        ele_out%value(vkick$) = ele_out%value(vkick$) + ele_in%value(vkick$)
        ele_out%value(x_limit$) = ele_in%value(x_limit$)
        ele_out%value(y_limit$) = ele_in%value(y_limit$)

        o_key = ele_out%key 
        if (ele_in%key == drift$ .and. &
                            (o_key == drift$ .or. o_key == marker$)) then
          ele_out%name = 'DRIFT_HYBRID' 
          ele_out%key = drift$
        else
          ele_out%name = 'HYBRID'
          ele_out%key = hybrid$
        endif

        ele_out%x       = ele_in%x
        ele_out%y       = ele_in%y
        ele_out%c_mat   = ele_in%c_mat
        ele_out%gamma_c = ele_in%gamma_c

      endif

      if (ele_out%key == hybrid$ .and. .not. do_taylor) then
        if (present(orb0_)) then
          ele_out%vec0 = orb0_(j_in)%vec - matmul(ele_out%mat6, c0%vec)
        else
          ele_out%vec0 = 0
        endif
      endif

    endif ! keep_ele

  enddo

! end cleanup

  if (ele_out%key == hybrid$ .and. z_decoupled)  &
                          call mat6_dispersion (ele_out%mat6, e_vec)

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
    ele_in => r_in%ele_(j_in)    
    if (keep_ele(j_in)) then
      i_out = i_out + 1
      ele_out => r_out%ele_(i_out)
      ix_out(j_in) = i_out
      ele_out = ele_in
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
    ele_out => r_out%ele_(i_out)

    do j = ele_out%ix1_slave, ele_out%ix2_slave
      ix = r_in%control_(j)%ix_slave
      if (ix_out(ix) == 0) then
        if (r_out%param%symmetry /= ew_antisymmetry$ .or.  &
                         ix <= n_ele .or. ix > r_in%n_ele_ring) then
          type '(a, /, 2a, /, 2a)',  &
      ' WARNING IN MAKE_HYBRID_RING: SLAVE ELEMENT LIST NOT COMPLETE.',  &
      '         FOR LORD ELEMENT: ', ele_out%name,  &
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

    do j = ele_out%ic1_lord, ele_out%ic2_lord
      ic = r_out%ic_(j)
      ix = r_in%control_(ic)%ix_lord
      if (ix == 0) then
        type *, 'WARNING IN MAKE_HYBRID_RING: LORD ELEMENT ',  &
                    'LIST NOT COMPLETE FOR: ', ele_out%name
      endif
    enddo

  enddo

! end

  if (r_out%n_ele_ring == 0) then
    type *, 'ERROR IN MAKE_HYBRID_RING: OUTPUT RING HAS 0 ELEMENTS!'
    call err_exit
  endif

  r_out%ele_init = r_in%ele_init

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
!   E_VEC(4) -- Real(rdef): eta vector
!
! Output:
!   mat6(6,6) -- Real(rdef): Matrix with 4x4 x-y submatrix already made.
!-



subroutine mat6_dispersion (mat6, e_vec)

  use bmad

  implicit none

  real(rdef), intent(inout) :: mat6(6,6)
  real(rdef), intent(in) :: e_vec(:)

  real(rdef) e2_vec(4)

  integer i

!

  mat6(1:4, 6) = e_vec(1:4)

  e2_vec(1) = -e_vec(2)
  e2_vec(2) =  e_vec(1)
  e2_vec(3) = -e_vec(4)
  e2_vec(4) =  e_vec(3)

  mat6(5,1:4) = matmul (e2_vec, mat6(1:4,1:4))

end subroutine
