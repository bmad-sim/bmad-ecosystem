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
      integer stack(*)
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
    subroutine seq_expand1 (b_line, seq_ele)
      use local_bmad_struct
      implicit none
      type (seq_struct) b_line
      type (seq_ele_struct) seq_ele(n_ele_maxx)
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
      character*(*) what2
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
      integer an_indexx(*)
      character*16 name
      character*16 names(*)
    end subroutine
  end interface
 


end module
