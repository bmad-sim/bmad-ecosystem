!+
! Subroutine multipole_init (ele, who, zero)
!
! Subroutine to allocate memory for the ele%a_pole and ele%b_pole magnetic multipoles
! and/or ele%a_pole_elec and ele%b_ple_elec electric multipoles
!
! Input:
!   who  -- Integer: electric$, magnetic$, or all$
!   zero -- Logical, optional: If present and True then zero the arrays
!             even if they already exist when this routine is called. 
!             Default is False which means that if the arrays already 
!             exist then this routine will do nothing.
!
! Output:
!   ele -- Ele_struct: Element holding the multipoles.
!     %a_pole(0:n_pole_maxx) -- Multipole An array 
!     %b_pole(0:n_pole_maxx) -- Multipole Bn array
!-

subroutine multipole_init (ele, who, zero)

use bmad_interface, except_dummy => multipole_init
  
implicit none

type (ele_struct) ele
integer who
logical, optional :: zero

!

if (who == magnetic$ .or. who == all$) then
  if (allocated(ele%multipole_cache)) ele%multipole_cache%ix_pole_mag_max = invalid$

  ! If %a_pole and %b_pole already exist then zero them if zero argument present and True.

  if (associated (ele%a_pole)) then
    if (logic_option(.false., zero)) then
      ele%a_pole = 0
      ele%b_pole = 0
      ele%multipole_cache = multipole_cache_struct()
    endif

  ! If memory not allocated then allocate and zero.
  else
    allocate (ele%a_pole(0:n_pole_maxx), ele%b_pole(0:n_pole_maxx))
    if (.not. allocated(ele%multipole_cache)) allocate(ele%multipole_cache)
    ele%a_pole = 0
    ele%b_pole = 0
    ele%multipole_cache = multipole_cache_struct()
  endif
endif

!

if (who == electric$ .or. who == all$) then
  if (allocated(ele%multipole_cache)) ele%multipole_cache%ix_pole_elec_max = invalid$

  ! If %a_pole_elec and %b_pole_elec already exist then zero them if zero argument present and True.

  if (associated (ele%a_pole_elec)) then
    if (logic_option(.false., zero)) then
      ele%a_pole_elec = 0
      ele%b_pole_elec = 0
      ele%multipole_cache = multipole_cache_struct()
    endif

  ! If memory not allocated then allocate and zero.
  else
    allocate (ele%a_pole_elec(0:n_pole_maxx), ele%b_pole_elec(0:n_pole_maxx))
    ele%a_pole_elec = 0
    ele%b_pole_elec = 0
    ele%multipole_cache = multipole_cache_struct()
  endif
endif

end subroutine multipole_init

