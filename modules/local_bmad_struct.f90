module local_bmad_struct

! This is for bmad_parser and bmad_parser2

  use bmad_struct

! structure for a decleared variable

  type parser_var_struct
    character*16 name              ! variable name
    real value                     ! variable value
  end type                      

! Head structure for lines and lists (sequences)

  type seq_struct
    character*16 name              ! name of sequence
    integer type                   ! either LINE$ or LIST$
    integer ix1                    ! start index to SEQ_ELE array
    integer ix2                    ! stop index to SEQ_ELE array
    integer ix                     ! current index of element in SEQ_ELE
    integer indexx                 ! alphabetical order sorted index
  end type

! structure for holding the contents of lines and lists (sequences)

  type seq_ele_struct
    character*16 name              ! name of element, subline, or sublist
    integer type                   ! LINE$, LIST$, or ELEMENT$
    integer ix_array               ! if an element: pointer to ELE_ array
                                   ! if a list: pointer to SEQ_ array
    logical reflection             ! reflection of subline?
  end type

! structure for holding the control names and pointers for
! superimpose and overlay elements

  type parser_ele_struct
    character*16 name(100)
    character*16 attrib_name(100)
    type (control_struct)  cs(100)
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

  integer line$, list$, element$
  parameter (line$    = 1)
  parameter (list$    = 2)
  parameter (element$ = 3)

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
    character*200 current_file_name, file_name_(20)
    logical parser_debug, no_digested, error_flag
    integer iseq_tot, ivar_tot, ivar_init
  end type

! common stuff

  type (pcom_struct)  pcom
  type (parser_var_struct)  var_(ivar_maxx)

  character*16 :: blank = ' '

!--------------------------------------------------

!$Id$
!$Log$
!Revision 1.3  2001/10/05 18:24:28  rwh24
!Increased size of filename variables to 200 chars
!
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

! this is for multi_turn_tracking_to_mat

  type (coord_struct), pointer :: multi_turn_func_common(:)  

end module
