!+
!-

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:45:23  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

module local_bmad_interface

  use bmad_interface

  interface
    subroutine get_next_word (word, ix_word, delim_list, &
                                       delim, delim_found, upper_case_word)
      use local_bmad_struct
      implicit none
      integer ix_word
      character*(*) word
      character*(*) delim_list
      character*(*) delim
      logical delim_found, upper_case_word
    end subroutine
  end interface
 
  interface
    subroutine file_stack (how, file_name, finished)
      use local_bmad_struct
      implicit none
      character*(*) how
      character*(*) file_name
      logical finished
    end subroutine
  end interface
 
  interface
    subroutine load_parse_line (how, ix_cmd, file_end)
      use local_bmad_struct
      implicit none
      integer ix_cmd
      character*(*) how
      logical file_end
    end subroutine
  end interface
 
  interface
    subroutine evaluate_value (var, ring, final_delim, final_delim_found, err_flag)
      use local_bmad_struct
      implicit none
      type (parser_var_struct) var
      type (ring_struct) ring
      character final_delim*1
      logical final_delim_found
      logical err_flag
    end subroutine
  end interface
 
  interface
    subroutine increment_pointer (ix, reflect)
      implicit none
      integer ix
      integer reflect
    end subroutine
  end interface
 
  interface
    subroutine pushit (stack, i_lev, value)
      implicit none
      integer stack(:)
      integer i_lev
      integer value
    end subroutine
  end interface
 
  interface
    subroutine word_to_value (word, ring, value)
      use local_bmad_struct
      implicit none
      type (ring_struct) ring
      real value
      character*(*) word
    end subroutine
  end interface
 
  interface
    subroutine get_attribute (how, ele, ring, pring, &
                                                 delim, delim_found, err_flag)
      use local_bmad_struct
      implicit none
      type (ring_struct) ring
      type (parser_ring_struct) pring
      type (ele_struct) ele
      integer how
      character delim*1
      logical delim_found
      logical err_flag
    end subroutine
  end interface
 
  interface
    subroutine type_get (ele, ix_type, delim, delim_found)
      use local_bmad_struct
      implicit none
      type (ele_struct) ele
      integer ix_type
      character delim*1
      logical delim_found
    end subroutine
  end interface
 
  interface
    subroutine get_overlay_group_names (ele, ring, pring, delim, delim_found)
      use local_bmad_struct
      implicit none
      type (ele_struct) ele
      type (parser_ring_struct) pring
      type (ring_struct) ring
      character delim*1
      logical delim_found
    end subroutine
  end interface
 
  interface
    subroutine verify_valid_name (name, ix_name)
      implicit none
      integer ix_name
      character*(*) name
    end subroutine
  end interface
 
  interface
    subroutine error_exit (what1, what2)
      use local_bmad_struct
      implicit none
      character*(*) what1
      character*(*) what2
    end subroutine
  end interface
 
  interface
    subroutine warning (what1, what2)
      use local_bmad_struct
      implicit none
      character*(*) what1
      character*(*), optional :: what2
    end subroutine
  end interface
 
  interface
    subroutine init_bmad_parser_common
      use local_bmad_struct
    end subroutine
  end interface
 
  interface
    subroutine compute_super_lord_s (ele, ring, pring)
      use local_bmad_struct
      implicit none
      type (ring_struct) ring
      type (ele_struct) ele
      type (parser_ring_struct) pring
    end subroutine
  end interface
 
  interface
    subroutine add_all_superimpose (ring, ele_in, pele)
      use local_bmad_struct
      implicit none
      type (ring_struct) ring
      type (ele_struct) ele_in
      type (parser_ele_struct) pele
    end subroutine
  end interface
 
  interface
    subroutine compute2_super_lord_s (ring, i_ref, ele, pele)
      use local_bmad_struct
      implicit none
      type (ring_struct) ring
      type (ele_struct) ele
      type (parser_ele_struct) pele
      integer i_ref
    end subroutine
  end interface
 
  interface
    subroutine find_indexx (name, names, an_indexx, n_max, ix_match)
      implicit none
      integer n_max
      integer ix_match
      integer an_indexx(:)
      character*16 name
      character*16 names(:)
    end subroutine
  end interface
 
!-------------------------------------------------------------------------
contains

!+
! Subroutine seq_expand1 (seq_, ix_seq, top_level)
!
! Subroutine to expand a sequence
!-

recursive subroutine seq_expand1 (seq_, ix_seq, top_level)

  use local_bmad_struct

  implicit none

  type (seq_struct), target :: seq_(:)
  type (seq_struct), pointer :: seq
  type (seq_ele_struct) s_ele(n_ele_maxx)

  integer ix_seq, ix_ele, ix_word, ix, count, i, j, n, ios, i_rl
  integer, save :: ix_internal = 0

  character*20 word
  character delim*1, str*16, name*16
              
  logical delim_found, top_level, replacement_line_here

! first thing should be a "("

  seq => seq_(ix_seq)
  ix_ele = 1

  call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

  if (delim /= '(') call warning  &
        ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
  if (ix_word /= 0)  call warning  &
        ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
        'FOR LINE: ' // seq%name)

! now parse list proper

  do 

    call get_next_word (word, ix_word, ':=(,)', delim, delim_found, .true.)

    ix = index(word, '*')          ! E.g. word = '-3*LINE'
    if (ix /= 0) then
      read (word(:ix-1), *, iostat = ios) count
      if (ios /= 0) then
        call warning ('BAD COUNT: ' // word, ' ')
        return
      endif
      s_ele(ix_ele)%name = word(ix+1:)
      if (count < 0) then
        s_ele(ix_ele)%reflect = .true.
      else
        s_ele(ix_ele)%reflect = .false.
      endif
      count = abs(count)
      ix_word = ix_word - ix
    elseif (word(1:1) == '-') then
      s_ele(ix_ele)%reflect = .true.
      count = 1
      s_ele(ix_ele)%name = word(2:)
      ix_word = ix_word - 1
    else
      s_ele(ix_ele)%reflect = .false.
      count = 1
      s_ele(ix_ele)%name = word
    endif

! check for a subline or replacement line
! if there is one then save as an internal sequence

    name = s_ele(ix_ele)%name
    if (name /= ' ') call verify_valid_name (name, ix_word)

    replacement_line_here = .false.
    if (delim == '(') then ! subline or replacement line
      if (name /= ' ') replacement_line_here = .true.
      ix_internal = ix_internal + 1
      write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name
      s_ele(ix_ele)%name = str
      ix_seq = ix_seq + 1
      seq_(ix_seq)%name = str
      seq_(ix_seq)%type = line$
      pcom%parse_line = '(' // pcom%parse_line 
      call seq_expand1 (seq_, ix_seq, .false.)
      call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
      if (word /= ' ') call warning &
                ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
                 word, 'IN THE SEQUENCE: ' // seq%name)
    endif

! if a replacement line then switch the real list for the actual args.

    if (replacement_line_here) then  ! replacement line 
      do i_rl = 1, ix_seq
        if (i_rl == ix_seq) then
          call warning ('CANNOT FIND REPLACEMENT LINE DEFINITION FOR: ' &
                          // name, 'WHICH APPEARS IN LINE: ' // seq%name)
          return
        endif
        if (seq_(i_rl)%name == name) exit
      enddo

      if (seq_(i_rl)%type /= replacement_line$) then
        call warning (trim(name) // ' IS USED AS A REPLACEMENT LINE IN: ' // &
              seq%name, 'BUT IT IS NOT A REPLACEMENT LINE')
        return
      endif

      if (size(seq_(i_rl)%arg) /= size(seq_(ix_seq)%ele)) then
        call warning ('NUMBER OF ARGUMENTS IN REPLACEMENT LINE: ' &
                // seq_(i_rl)%name, 'DOES NOT MATCH USE IN LINE: ' // seq%name)
        return
      endif

      seq_(i_rl)%arg(:)%actual_name = seq_(ix_seq)%ele(:)%name
      n = size(seq_(i_rl)%ele)
      deallocate (seq_(ix_seq)%ele)
      allocate (seq_(ix_seq)%ele(n))
      do j = 1, n
        ix = seq_(i_rl)%ele(j)%ix_arg
        if (ix > 0) then
          seq_(ix_seq)%ele(j)%name = seq_(i_rl)%arg(ix)%actual_name
        else
          seq_(ix_seq)%ele(j) = seq_(i_rl)%ele(j)
        endif
      enddo 

   endif

! if a replacement line then look for element in argument list

    s_ele(ix_ele)%ix_arg = 0
    if (seq%type == replacement_line$) then
      do i = 1, size(seq%arg)
        if (seq%arg(i)%dummy_name == s_ele(ix_ele)%name) then
          s_ele(ix_ele)%ix_arg = i
          exit
        endif
      enddo
    endif

! if multiple repetition count then expand

    do i = 1, count-1
      s_ele(ix_ele+i) = s_ele(ix_ele)
    enddo

    ix_ele = ix_ele + count

    if (delim == ')') exit

    if (delim /= ',') then
      call warning  &
               ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
      exit
    endif
           
  enddo

! make sure there is nothing else if at top level

  if (top_level) then
    call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)

    if (delim_found .or. ix_word /= 0) call warning  &
          ('EXTRA CHARACTERS AFTER CLOSING ")"',  'FOR LINE: ' // seq%name)
  endif

! transfer

  ix_ele = ix_ele - 1
  allocate (seq%ele(ix_ele))
  seq%ele(:) = s_ele(1:ix_ele)

  seq%ix = 1   ! Init. Used for replacement list index

end subroutine

end module
