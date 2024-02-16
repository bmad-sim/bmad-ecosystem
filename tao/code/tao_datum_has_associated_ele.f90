!+
! Function tao_datum_has_associated_ele (data_type, branch_geometry) result (has_associated_ele)
!
! Routine to determine if a given element type has an associated lattice element.
!
! Input:
!   data_type        -- character(*): Type of data.
!   branch_geometry  -- integer, optional: Geometry of the associated lattice branch. open$ or closed$.
!
! Output:
!   has_associated_ele -- integer:
!                           no$          -- Must not have an associated lattice element.
!                           yes$         -- Must have an associated lattice element.
!                           maybe$       -- Is possible to have an associated lattice element or not.
!                           provisional$ -- Does not have an associated lattice element for closed geometies.
!                                            Only used if the branch_geometry argument is not present.
!-

function tao_datum_has_associated_ele (data_type, branch_geometry) result (has_associated_ele)

use bmad_struct

implicit none

character(*) data_type
character(len(data_type)) short_type
integer has_associated_ele
integer, optional :: branch_geometry
integer ix

!

ix = index(data_type, '.')
if (ix == 0) then
  short_type = data_type
else
  short_type = data_type(:ix-1)
endif

if ((short_type == 'unstable' .and. data_type /= 'unstable.orbit') .or. short_type == 'normal' .or. &
              data_type(1:18) == 'spin.polarization_' .or. data_type == 'spin.depolarization_rate' .or. &
              data_type == 'chrom.a' .or. data_type == 'chrom.b' .or. short_type == 'chrom_ptc' .or. &
              data_type(1:12) == 'chrom.dtune.' .or. short_type == 'slip_factor_ptc' .or. &
              short_type == 'srdt' .or. short_type == 'damp' .or. short_type == 'tune' .or. &
              short_type == 'momentum_compaction_ptc') then
  has_associated_ele = no$

elseif (data_type == 'emit.a' .or. data_type == 'norm_emit.a' .or. data_type == 'emit.b' .or. &
        data_type == 'norm_emit.b') then
  if (present(branch_geometry)) then
    if (branch_geometry == closed$) then
      has_associated_ele = no$
    else
      has_associated_ele = yes$
    endif
  else
    has_associated_ele = provisional$
  endif

elseif (data_type == 'unstable.orbit' .or. short_type == 'dynamic_aperture') then
  has_associated_ele = maybe$
  if (integer_option(int_garbage$, branch_geometry) == closed$) has_associated_ele = no$

elseif ((data_type(1:11) == 'expression:' .and. index(data_type, 'ele::#[') == 0) .or. data_type(1:8) == 'rad_int.') then
  has_associated_ele = maybe$

else
  has_associated_ele = yes$
endif

end function tao_datum_has_associated_ele
