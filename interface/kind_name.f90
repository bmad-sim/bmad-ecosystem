!+
! function kind_name (this_kind)
!
! function to return the name of a PTC kind.
!
! Input:
!   this_kind -- Integer: 
!
! Output:
!   kind_name -- Character*20: 
!-

function kind_name (this_kind)

  use accelerator_struct

  implicit none

  integer this_kind
  character*20 kind_name

!

  select case (this_kind)
  case (kind0); kind_name  = 'KIND0'
  case (kind1); kind_name  = 'DRIFT1'
  case (kind2); kind_name  = 'DKD2 (Gen Element)' 
  case (kind3); kind_name  = 'KICKT3 (Thin Ele)'
  case (kind4); kind_name  = 'CAV4 (RF Cavity)'
  case (kind5); kind_name  = 'SOL5 (Solenoid)'
  case (kind6); kind_name  = 'KTK (Slow Thick)'
  case (kind7); kind_name  = 'TKTF (Fast Thick)'
  case (kind8); kind_name  = 'NSMI (Normal SMI)'
  case (kind9); kind_name  = 'SSMI (Skew SMI)'
  case (kind10); kind_name = 'TEAPOT (Sector Bend)'
  case (kindfitted); kind_name = 'FITTED KIND'
  case (kinduser1); kind_name = 'USER1 KIND'
  case (kinduser2); kind_name = 'USER2 KIND'
  case default; write (kind_name, '(a, i5)') 'UNKNOWN KIND!', this_kind 
  end select

end function
  

