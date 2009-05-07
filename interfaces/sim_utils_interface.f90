#include "CESR_platform.inc"

module sim_utils_interface

interface
  subroutine abs_sort (array, index, n)
    use precision_def
    implicit none
    real(rp) array(*)
    integer index(*)
    integer n
  end subroutine
end interface
 
interface
  subroutine bbi_kick (x, y, r, kx, ky)
    use precision_def
    implicit none
    real(rp) x
    real(rp) y
    real(rp) r
    real(rp) kx
    real(rp) ky
  end subroutine
end interface
 
interface
  subroutine bracket_index (s_arr, i_min, i_max, s, ix)
    use precision_def
    implicit none
    integer i_min, i_max
    real(rp) s_arr(i_min:), s
    integer ix
  end subroutine
end interface
 
interface
  subroutine calc_file_number (file_name, num_in, num_out, err_flag)
    implicit none
    character(*) file_name
    integer num_in
    integer num_out
    logical err_flag
  end subroutine
end interface
 
interface
  subroutine change_file_number (file_name, change)
    implicit none
    character(*) file_name
    integer change
  end subroutine
end interface
 
interface
  subroutine complex_error_function (wr, wi, zr, zi)
    use precision_def
    implicit none
    real(rp) wr
    real(rp) wi
    real(rp) zr
    real(rp) zi
  end subroutine
end interface
 
interface
  function determinant (mat) result (det)
    use precision_def
    implicit none
    real(rp) mat(:,:)
    real(rp) det
  end function
end interface
 
interface
  subroutine doubleup_quotes (str_in, str_out, quote)
    character(*) str_in
    character(*) str_out
    character quote*1
  end subroutine
end interface
 
interface
  function even (num) result (is_even)
    implicit none
    integer num
    logical is_even
  end function
end interface
 
interface
  subroutine fff_sub(line, error)
    implicit none
    character(*) line
    logical error
  end subroutine
end interface
 
interface
  subroutine file_directorizer (in_file, out_file, directory, add_switch)
    implicit none
    logical add_switch
    character(*) in_file
    character(*) out_file
    character(*) directory
  end subroutine
end interface
 
interface
  subroutine file_get (string, dflt_file_name, file_name)
    implicit none
    character(*) file_name
    character(*) dflt_file_name
    character(*) string
  end subroutine
end interface
 
interface
  subroutine file_get_open (string, dflt_file_name, &
                                    file_name, file_unit, readonly)
    implicit none
    character(*) file_name
    character(*) dflt_file_name
    character(*) string
    integer file_unit
    logical readonly
  end subroutine
end interface
 
interface
  subroutine file_suffixer (in_file_name, out_file_name, suffix, add_switch)
    implicit none
    logical add_switch
    character(*) in_file_name
    character(*) out_file_name
    character(*) suffix
  end subroutine
end interface
 
interface
  subroutine type_this_file(filename)
    implicit none
    character(*) filename
  end subroutine
end interface
 
interface
  subroutine get_file_number (file_name, cnum_in, num_out, err_flag)
    implicit none
    character(*) file_name
    character(*) cnum_in
    integer num_out
    logical err_flag
  end subroutine
end interface
 
interface
  subroutine get_next_number (filein, cnum, digits)
    implicit none
    character(*) filein
    character(*) cnum
    integer digits
  end subroutine
end interface

interface
  subroutine if_error (idelim, icmd, error_string, line_number, end_check)
    implicit none
    integer idelim
    integer line_number
    integer end_check
    integer icmd
    character(*) error_string
  end subroutine
end interface
 
interface
  subroutine increment_file_number (file_name, digits, number, cnumber)
    implicit none
    character(*) file_name
    character(*) cnumber
    integer digits
    integer number
  end subroutine
end interface
 
interface
  function index_nocase(string1, string2) result (indx)
    implicit none
    integer indx
    character(*) string1
    character(*) string2
  end function
end interface
 
interface
  function integer_read(error_message) result (int_read)
    implicit none
    integer int_read
    character(*) error_message
  end function
end interface
 
interface
  function inverse (funct, y, x1, x2, tol) result (x)
    use precision_def
    implicit none
    real(rp) x 
    real(rp) y, x1, x2, tol
    interface
      function funct(x)
        use precision_def
        real(rp) funct, x
      end function funct
     end interface
  end function
end interface
 
interface
  function inverse_prob (val) result (prob)
    use precision_def
    implicit none
    real(rp) prob 
    real(rp) val
  end function
end interface
 
interface
  subroutine ion_kick (x, y, x_kicker, y_kicker, s_kicker)
    use precision_def
    real(rp) x
    real(rp) y
    real(rp) x_kicker
    real(rp) y_kicker
    real(rp) s_kicker
  end subroutine
end interface
 
interface
  subroutine ion_kick_2d (x, y, x_kicker, y_kicker)
    use precision_def
    real(rp) x
    real(rp) y
    real(rp) x_kicker
    real(rp) y_kicker
  end subroutine
end interface
 
interface
  function is_integer (string) result (valid)
    character(*) string
    logical valid
  end function
end interface

interface
  function is_logical (string, ignore) result (valid)
    character(*) string
    logical, optional :: ignore
    logical valid
  end function
end interface

interface
  function is_real (string, ignore) result (valid)
    character(*) string
    logical, optional :: ignore
    logical valid
  end function
end interface

interface
  subroutine linear_fit (x, y, n_data, a, b, sig_a, sig_b)
    use precision_def
    implicit none
    integer n_data
    real(rp) x(*)
    real(rp) y(*)
    real(rp) a
    real(rp) b
    real(rp) sig_a
    real(rp) sig_b
  end subroutine
end interface
 
interface
  subroutine location_decode(string, array, ix_min, num, names, exact_case)
    implicit none
    integer num
    integer ix_min
    character(*) string
    logical array(ix_min:)
    character(*), optional :: names(ix_min:)
    logical, optional :: exact_case
  end subroutine
end interface

interface
  subroutine make_legal_comment (comment_in, comment_out)
    implicit none
    character(*) comment_in
    character(*) comment_out
  end subroutine
end interface
 
interface
  function match_wild (string, template) result (is_match)
    implicit none
    logical is_match
    character(*) string
    character(*) template
  end function
end interface
 
interface
  subroutine mat_det(mat, det)
    use precision_def
    implicit none
    real(rp) mat(:,:)
    real(rp) det
  end subroutine
end interface
 
interface
  subroutine mat_rotation (mat, angle, bet_1, bet_2, alph_1, alph_2)
    use precision_def
    implicit none
    real(rp) angle
    real(rp) bet_1, bet_2
    real(rp) alph_1, alph_2
    real(rp) mat(:,:)
  end subroutine
end interface
 
interface
  subroutine mat_symp_conj(mat1, mat2)
    use precision_def
    implicit none
    real(rp) mat1(:,:)
    real(rp) mat2(:,:)
  end subroutine
end interface
 
interface
  subroutine mat_type (mat, nunit, header)
    use precision_def
    implicit none
    real(rp) mat(:,:)
    integer, optional, intent(in) :: nunit
    character(*), optional, intent(in) :: header
  end subroutine
end interface
 
interface
  subroutine mat_make_unit (mat)
    use precision_def
    real(rp) mat(:,:)
  end subroutine
end interface
 
interface
  subroutine node_put (node, n1, n2, val_in, cmd_only, val_out, bad_set)
    implicit none
    integer n1
    integer n2
    integer val_in(*)
    integer val_out(*)
    character(12) node
    logical cmd_only
    logical bad_set
  end subroutine
end interface
 
interface
  function odd (num) result (is_odd)
    implicit none
    integer num
    logical is_odd
  end function
end interface
 
interface
  function probability_funct(x) result (prob)
    use precision_def
    implicit none
    real(rp) prob
    real(rp) x
  end function
end interface
 
interface
  subroutine query_string (query_str, upcase, return_str, ix, ios)
    implicit none
    character(*) return_str
    character(*) query_str
    integer ix
    integer ios
    logical upcase
  end subroutine
end interface
 
interface
  function real_read(error_message)
    use precision_def
    implicit none
    real(rp) real_read
    character(*) error_message
  end function
end interface
 
interface
  subroutine run_timer(command, time)
    use precision_def
    implicit none
    real(rp), optional :: time
    character(*) command
  end subroutine
end interface
 
interface
  subroutine skip_header (ix_unit, error_flag)
    implicit none
    integer ix_unit
    logical error_flag
  end subroutine
end interface
 
interface
  subroutine spline_akima (spline, ok)
    use sim_utils_struct
    implicit none
    type (spline_struct) :: spline(:)
    logical ok
  end subroutine
end interface

interface
  subroutine spline_evaluate (spline, x, ok, y, dy)
    use sim_utils_struct
    type (spline_struct), target :: spline(:)
    real(rp) :: x
    real(rp), optional :: y
    real(rp), optional :: dy
    logical ok
  end subroutine
end interface
 
interface
  subroutine spawn_command (command)
    implicit none
    character(*) command
  end subroutine
end interface
 
interface
  subroutine string_to_int (line, default, value, err_flag)
    implicit none
    integer default
    integer value
    character(*) line
    logical err_flag
  end subroutine
end interface
 
interface
  subroutine string_to_real (line, default, value, err_flag)
    use precision_def
    implicit none
    real(rp) default
    real(rp) value
    character(*) line
    logical err_flag
  end subroutine
end interface
 
interface
  subroutine string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next)
    implicit none
    character(*) in_str
    character(*) out_str
    character(*) delimitors
    character(1) delim
    integer ix_word
    integer ix_next
  end subroutine
end interface
 
interface
  subroutine test_tune_tracker_lock (tracker_locked)
    implicit none
    logical tracker_locked(2)
  end subroutine
end interface
 
interface
  subroutine downcase_string(string)
    implicit none
    character(*) string
  end subroutine
end interface
 
interface
  subroutine upcase_string(string)
    implicit none
    character(*) string
  end subroutine
end interface
 
interface
  function word_len (wording) result (wlen)
    implicit none
    integer wlen 
    character(*) wording
  end function
end interface
 
interface
  subroutine word_read (in_str, delim_list, word, &
                                       ix_word, delim, delim_found, out_str)
    implicit none
    character(*) in_str
    character(*) out_str
    character(*) word
    character(*) delim_list
    character(1) delim
    integer ix_word
    logical delim_found
  end subroutine
end interface
 
interface
  subroutine str_substitute (string, str_match, str_replace, do_trim)
    implicit none
    character(*) string
    character(*), optional :: str_match, str_replace
    logical, optional :: do_trim
  end subroutine
end interface

interface
  subroutine date_and_time_stamp (string, numeric_month)
    implicit none
    character(*) string
    logical, optional :: numeric_month
  end subroutine
end interface
 
interface
  subroutine find_file (file_in, found, file_out, dirs)
    character(*) file_in, dirs(:), file_out
    logical found
  end subroutine
end interface

interface                   
  subroutine get_file_time_stamp (file, time_stamp)
    implicit none
    character(*) file, time_stamp
  end subroutine
end interface

interface                   
   function lunget ()
     implicit none
     integer lunget
   end function lunget
end interface

interface
   function match_reg(str, pat) result (is_match)
     implicit none
     logical is_match
     character(*) str, pat
   end function match_reg
end interface

interface
   subroutine milli_sleep (milli_sec)
    implicit none
    integer milli_sec
  end subroutine
end interface    

interface
  subroutine ps2gif (ps_file, gif_file, kill_ps_file)
    implicit none
    character(*) ps_file, gif_file
    logical, optional :: kill_ps_file
  end subroutine
end interface    

interface
   recursive function str_match_wild(str, pat) result (a_match)
     implicit none
     character(*) pat, str
     logical a_match
   end function str_match_wild
end interface

interface
   subroutine str_upcase(dst, src)
     implicit none
     character(*) dst, src
   end subroutine str_upcase
end interface

interface
   subroutine str_downcase(dst, src)
     implicit none
     character(*) dst, src
   end subroutine str_downcase
end interface

interface
   subroutine system_command (line)
     implicit none
     character(*) line
   end subroutine
end interface

interface
   subroutine string_trim (in_string, out_string, word_len)
     implicit none
     character(*) in_string, out_string
     integer word_len
   end subroutine string_trim
end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
