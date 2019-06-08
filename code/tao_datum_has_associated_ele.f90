!+
! Function tao_datum_has_associated_ele (data_type) result (has_associated_ele)
!
! Routine to determine if a given element type has an associated lattice element.
!
! Input:
!   data_type     -- character(*): Type of data.
!
! Output:
!   has_associated_ele -- integer:
!                           no$          -- Does not have an associated lattice element.
!                           yes$         -- Has an associated lattice element.
!                           provisional$ -- Does not have an associated lattice element for circular lbranches.
!                           maybe$       -- Is possible to have an associated lattice element or not.
!-

function tao_datum_has_associated_ele (data_type) result (has_associated_ele)

use sim_utils

implicit none

character(*) data_type
integer has_associated_ele

!

if (data_type == 'unstable.orbit' .or. data_type(1:7) == 'normal.' .or. data_type(1:5) == 'srdt.' .or. &
                      data_type(1:18) == 'spin.polarization_' .or. data_type == 'spin.depolarization_rate') then
  has_associated_ele = no$

elseif (data_type == 'chrom.a' .or. data_type == 'chrom.b' .or. &
             data_type(1:12) == 'chrom.dtune.' .or. data_type(1:5) == 'damp.' .or. &
             data_type(1:17) == 'multi_turn_orbit.' .or. data_type(1:5) == 'tune.' .or. &
             data_type(1:13) == 'unstable.ring' .or. index(data_type, 'emit.') /= 0) then
  has_associated_ele = provisional$

elseif ((data_type(1:11) == 'expression:' .and. index(data_type, 'ele::#[') == 0) .or. data_type(1:8) == 'rad_int.') then
  has_associated_ele = maybe$

else
  has_associated_ele = yes$
endif

end function tao_datum_has_associated_ele
