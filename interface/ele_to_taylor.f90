!+
! Subroutine ele_to_taylor (ele, orb0, param)
!
! Subroutine to make a taylor map for an element. 
! The order of the map is set by set_ptc
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Element_struct: 
!     %integration_order  -- Order for the symplectic integrator: 2, 4, or 6.
!     %num_steps          -- Number of integrater steps.
!   orb0  -- Coord_struct, optional: Starting coords around which the Taylor series 
!              is evaluated.
!   param -- Param_struct, optional: 
!     %energy -- Needed for wigglers.
!
! Output:
!   ele -- Element_struct:
!     %taylor(6)  -- Taylor maps.
!-

subroutine ele_to_taylor (ele, orb0, param)

  use accelerator
  
  implicit none
  
  type (ele_struct), intent(inout) :: ele
  type (coord_struct), optional, intent(in) :: orb0
  type (param_struct), optional, intent(in) :: param

  type (fibre), pointer, save :: a_fibre
  type (real_8) y(6), y2(6)
  type (universal_taylor), save :: u_taylor(6)

  real(dp) x(6)
  
  integer i
  
  logical :: init_needed = .true.

! Init

  if (init_needed) then
    allocate (a_fibre)
    do i = 1, 6
      u_taylor(i) = 0  ! nullify
    enddo
    init_needed = .false.
  endif

  if (bmad_com%taylor_order_ptc == 0) then
    call set_ptc (taylor_order = bmad_com%taylor_order)
  endif

! Track with offset

  call alloc_fibre (a_fibre)
  call ele_to_fibre (ele, a_fibre, param)
 
  if (present(orb0)) then
    ele%taylor(:)%ref = orb0%vec
    call vec_bmad_to_ptc (orb0%vec, x)
  else
    ele%taylor(:)%ref = 0
    x = 0
  endif

  call real_8_init(y)
  y = x
  call track (a_fibre, y, default, +1)

! take out the offset

  call real_8_init(y2)
  y2 = -x
  call concat_real_8 (y, y2, y)

! convert to bmad_taylor  

  do i = 1, 6
    u_taylor(i) = y(i)%t
  enddo
  
  call universal_to_bmad_taylor (u_taylor, ele%taylor)
  ele%taylor_order = bmad_com%taylor_order_ptc

  call kill(a_fibre)
  call kill(y)
  call kill(y2)

  if (associated (ele%gen_field)) call kill_gen_field (ele%gen_field)

end subroutine
