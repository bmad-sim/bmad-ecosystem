!+
!-

!$Id$
!$Log$
!Revision 1.4  2002/01/08 21:45:23  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.3  2001/10/05 18:24:28  rwh24
!Increased size of filename variables to 200 chars
!
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

module local_bmad_struct

! This is for bmad_parser and bmad_parser2

  use bmad_struct

! structure for a decleared variable

  type parser_var_struct
    character*16 name              ! variable name
    real value                     ! variable value
  end type                      

! structure for holding the contents of lines and lists (sequences)

  type seq_ele_struct
    character*16 name              ! name of element, subline, or sublist
    integer type                   ! LINE$, LIST$, ELEMENT$
    integer ix_array               ! if an element: pointer to ELE_ array
                                   ! if a list: pointer to SEQ_ array
    integer ix_arg                 ! for replacement lines
    logical reflect                ! reflection of subline?
  end type

  type replacement_arg_struct
    character*16 dummy_name
    character*16 actual_name
  end type

! Head structure for lines and lists (sequences)

  type seq_struct
    character*16 name                  ! name of sequence
    integer type                       ! LINE$, REPLACEMENT_LINE$ or LIST$
    integer ix                         ! current index of element in %ELE
    integer indexx                     ! alphabetical order sorted index
    type (seq_ele_struct), pointer :: ele(:)
    type (replacement_arg_struct), pointer :: arg(:)
  end type

! stack structure

  type seq_stack_struct
    integer ix_seq                ! index to seq_(:) array
    integer ix_ele                ! index to seq%ele(:) array
    integer reflect               ! reflection sequence?
  end type

!-----------------------------------------------------------
! structure for holding the control names and pointers for
! superimpose and overlay elements

  type parser_ele_struct
    character*16 ref_name
    character*16, pointer :: name_(:)
    character*16, pointer :: attrib_name_(:)
    type (control_struct), pointer :: cs_(:)
    integer ix_count
    integer ele_pt, ref_pt
    logical common_lord
    real s
    integer indexx
  end type

  type parser_ring_struct
    type (parser_ele_struct) ele(0:n_ele_maxx)
  end type

! component_struct

  type component_struct
    character*16 name
    integer ix_ele
  end type

  integer ivar_maxx
  parameter (ivar_maxx =  500)

!

  integer, parameter :: line$ = 1, list$ = 2, element$ = 3
  integer, parameter :: replacement_line$ = 4

  integer begin$, center$, end$
  parameter (begin$ = -1)
  parameter (center$ = 0)
  parameter (end$ = 1)

  integer def$, redef$
  parameter (def$ = 1)
  parameter (redef$ = 2)

!

  type pcom_struct
    integer i_line, f_unit, n_files
    character*160 parse_line
    character*16 parser_name
    character*72 current_file_name, file_name_(20), debug_line
    logical parser_debug, no_digested, error_flag
    integer iseq_tot, ivar_tot, ivar_init
  end type

! common stuff

  type (pcom_struct)  pcom
  type (parser_var_struct)  var_(ivar_maxx)

  character*16 :: blank = ' '

!--------------------------------------------------
! this is for multi_turn_tracking_to_mat

  type (coord_struct), pointer :: multi_turn_func_common(:)  

end module
