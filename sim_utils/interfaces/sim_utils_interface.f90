module sim_utils_interface

use sim_utils_struct

interface set_parameter
  function set_parameter_real (param_val, set_val) result (save_val)
    import
    implicit none
    real(rp) param_val, set_val, save_val
  end function
  function set_parameter_int (param_val, set_val) result (save_val)
    import
    implicit none
    integer param_val, set_val, save_val
  end function
  function set_parameter_logic (param_val, set_val) result (save_val)
    import
    implicit none
    logical param_val, set_val, save_val
  end function
end interface

!--------------------------------------------------------------------------------

interface

elemental function cosc(x, nd) result (y)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) y
  integer, optional, intent(in) :: nd
end function

elemental function asinc(x, nd) result (y)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) y
  integer, optional, intent(in) :: nd
end function

function assert_equal (int_arr, err_str) result (ival)
  import
  implicit none
  integer, intent(in) :: int_arr(:)
  integer ival
  character(*) err_str
end function

function bracket_index (s, s_arr, i_min, dr, restrict) result (ix)
  import
  implicit none
  integer i_min, i_max
  real(rp) s_arr(i_min:), s
  real(rp), optional :: dr
  logical, optional :: restrict
  integer ix
end function

function bracket_index_int (s, int_arr, i_min, restrict) result (ix)
  import
  implicit none
  integer i_min, i_max
  integer int_arr(i_min:), s
  logical, optional :: restrict
  integer ix
end function

function bracket_index2 (s, ix0, s_arr, i_min, dr, restrict) result (ix)
  import
  implicit none
  integer i_min, i_max
  real(rp) s_arr(i_min:), s
  real(rp), optional :: dr
  integer ix0, ix
  logical, optional :: restrict
end function

subroutine calc_file_number (file_name, num_in, num_out, err_flag)
  implicit none
  character(*) file_name
  integer num_in
  integer num_out
  logical err_flag
end subroutine

subroutine change_file_number (file_name, change)
  implicit none
  character(*) file_name
  integer change
end subroutine

elemental function cos_one(angle) result (cos1)
  import
  implicit none
  real(rp), intent(in) :: angle
  real(rp) cos1
end function

subroutine cplx_mat_inverse(mat, mat_inv, ok, print_err)
  import
  implicit none
  complex(rp) :: mat(:,:)
  complex(rp) :: mat_inv(:,:)
  logical, optional :: ok, print_err
end subroutine

subroutine cplx_mat_make_unit (mat)
  import
  implicit none
  complex(rp) mat(:,:)
end subroutine

subroutine complex_error_function (wr, wi, zr, zi)
  import
  implicit none
  real(rp) wr
  real(rp) wi
  real(rp) zr
  real(rp) zi
end subroutine

function cross_product (a, b) result (c)
  import
  implicit none
  real(rp) a(:), b(:)
  real(rp) c(3)
end function

subroutine date_and_time_stamp (string, numeric_month, include_zone)
  implicit none
  character(*) string
  logical, optional :: numeric_month, include_zone
end subroutine

function determinant (mat) result (det)
  import
  implicit none
  real(rp) mat(:,:)
  real(rp) det
end function

subroutine display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res)
  import
  implicit none
  real(rp) x_size, y_size, x_res, y_res
  integer ix_screen
end subroutine

function dJ_bessel(m, arg) result (dj_bes)
  import
  implicit none
  integer m
  real(rp) arg, dj_bes
end function

function djb_hash (str, old_hash) result (hash)
  implicit none
  character(*) :: str
  integer :: hash
  integer, optional :: old_hash
end function

function djb_str_hash (in_str) result (hash_str)
  implicit none
  character(*) :: in_str
  character(6) :: hash_str
end function

subroutine doubleup_quotes (str_in, str_out, quote)
  character(*) str_in
  character(*) str_out
  character quote*1
end subroutine

subroutine downcase_string(string)
  implicit none
  character(*) string
end subroutine

function downcase(str_in) result (str_out)
  implicit none
  character(*) str_in
  character(len(str_in)) str_out
end function

subroutine err_exit(err_str)
  implicit none
  character(*), optional :: err_str
end subroutine

function factorial(n) result (fact)
  import
  implicit none
  real(rp) fact
  integer n
end function

subroutine fff_sub(line, error)
  implicit none
  character(*) line
  logical error
end subroutine

subroutine fft_1d (arr, isign)
  import
  complex(rp) arr(:)
  integer isign
end subroutine

subroutine file_directorizer (in_file, out_file, directory, add_switch)
  implicit none
  logical add_switch
  character(*) in_file
  character(*) out_file
  character(*) directory
end subroutine

subroutine file_get (string, dflt_file_name, file_name)
  implicit none
  character(*) file_name
  character(*) dflt_file_name
  character(*) string
end subroutine

subroutine file_get_open (string, dflt_file_name, file_name, file_unit, readonly)
  implicit none
  character(*) file_name
  character(*) dflt_file_name
  character(*) string
  integer file_unit
  logical readonly
end subroutine

subroutine file_suffixer (in_file_name, out_file_name, suffix, add_switch)
  implicit none
  logical add_switch
  character(*) in_file_name
  character(*) out_file_name
  character(*) suffix
end subroutine

function gen_complete_elliptic(kc, p, c, s, err_tol) result (value)
  import
  implicit none
  real(rp) kc, p, c, s, value
  real(rp), optional :: err_tol
end function

subroutine get_file_number (file_name, cnum_in, num_out, err_flag)
  implicit none
  character(*) file_name
  character(*) cnum_in
  integer num_out
  logical err_flag
end subroutine

subroutine get_next_number (filein, cnum, digits)
  implicit none
  character(*) filein
  character(*) cnum
  integer digits
end subroutine

subroutine get_file_time_stamp (file, time_stamp)
  implicit none
  character(*) file, time_stamp
end subroutine

function I_bessel(m, arg) result (i_bes)
  import
  implicit none
  integer m
  real(rp) arg, i_bes
end function

function I_bessel_extended(m, arg) result (i_bes)
  import
  implicit none
  integer m
  real(rp) arg
  complex(rp) i_bes
end function

subroutine increment_file_number (file_name, digits, number, cnumber)
  implicit none
  character(*) file_name
  character(*) cnumber
  integer digits
  integer number
end subroutine

function index_nocase(string1, string2) result (indx)
  implicit none
  integer indx
  character(*) string1
  character(*) string2
end function

function int_str(int, width) result (str)
  implicit none
  integer int
  integer, optional :: width
  character(:), allocatable :: str
end function

function inverse (funct, y, x1, x2, tol) result (x)
  import
  implicit none
  real(rp) x 
  real(rp) y, x1, x2, tol
  interface
    function funct(x)
      import
      real(rp) funct, x
    end function funct
  end interface
end function
 
function inverse_prob (val) result (prob)
  import
  implicit none
  real(rp) prob 
  real(rp) val
end function

function is_alphabetic (string, valid_chars) result (is_alpha)
  implicit none
  character(*) string
  character(*), optional :: valid_chars
  logical is_alpha
end function

function is_decreasing_sequence (array, strict) result (is_decreasing)
  import
  implicit none
  real(rp) array(:)
  logical, optional :: strict
  logical is_decreasing
end function

function is_increasing_sequence (array, strict) result (is_increasing)
  import
  implicit none
  real(rp) array(:)
  logical, optional :: strict
  logical is_increasing
end function

function is_integer (string, int, delims, ix_word) result (valid)
  implicit none
  character(*) string
  integer, optional :: int, ix_word
  character(*), optional :: delims
  logical valid
end function

function is_logical (string, ignore) result (valid)
  implicit none
  character(*) string
  logical, optional :: ignore
  logical valid
end function

function is_real (string, ignore, real_num) result (valid)
  import
  implicit none
  character(*) string
  real(rp), optional :: real_num
  logical, optional :: ignore
  logical valid
end function

! isatty is a standard POSIX C routine.

function isatty(ii) result (t_type) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int
  integer(c_int), value :: ii
  integer(c_int) :: t_type
end function

function J_bessel(m, arg) result (j_bes)
  import
  implicit none
  integer m
  real(rp) arg, j_bes
end function

subroutine linear_fit (x, y, n_data, a, b, sig_a, sig_b)
  import
  implicit none
  integer n_data
  real(rp) x(:)
  real(rp) y(:)
  real(rp) a
  real(rp) b
  real(rp) sig_a
  real(rp) sig_b
end subroutine

subroutine linear_fit_2d (x, y, z, coef)
  import
  implicit none
  real(rp) x(:), y(:), z(:)
  real(rp) coef(3)
end subroutine

subroutine location_decode(string, array, ix_min, num, names, exact_case, can_abbreviate, print_err)
  implicit none
  integer num
  integer ix_min
  character(*) string
  logical array(ix_min:)
  character(*), optional :: names(ix_min:)
  logical, optional :: exact_case, can_abbreviate, print_err
end subroutine

function logic_str(logic) result (str)
  implicit none
  logical logic
  character(1) :: str
end function

function lunget ()
  implicit none
  integer lunget
end function lunget

function match_reg(str, pat) result (is_match)
  implicit none
  logical is_match
  character(*) str, pat
end function match_reg

subroutine milli_sleep (milli_sec)
  implicit none
  integer milli_sec
end subroutine

subroutine make_legal_comment (comment_in, comment_out)
  implicit none
  character(*) comment_in
  character(*) comment_out
end subroutine

function match_wild (string, template) result (is_match)
  implicit none
  logical is_match
  character(*) string
  character(*) template
end function

subroutine mat_pseudoinverse(A,Ap,svd_condition,print_err,ok)
  import
  real(rp) A(:,:)
  real(rp) Ap(:,:)
  real(rp), optional :: svd_condition
  logical, optional :: print_err
  logical, optional :: ok
end subroutine

subroutine mat_eigen (mat, eigen_val, eigen_vec, error, print_err)
  import
  implicit none
  real(rp) mat(:,:), mat2(size(mat,1), size(mat,2))
  complex(rp) eigen_val(:), eigen_vec(:,:)
  logical error
  logical, optional :: print_err
end subroutine

subroutine mat_inverse (mat, mat_inv, ok, print_err)
  import
  real(rp) :: mat(:,:)
  real(rp) :: mat_inv(:,:)
  logical, optional :: ok, print_err
end subroutine

subroutine mat_rotation (mat, angle, bet_1, bet_2, alph_1, alph_2)
  import
  implicit none
  real(rp) angle
  real(rp) bet_1, bet_2
  real(rp) alph_1, alph_2
  real(rp) mat(:,:)
end subroutine

function mat_scale_p0 (mat_in, p0_ratio, invert) result (mat_out)
  import
  real(rp), intent(in) :: mat_in(:,:)
  real(rp), optional :: p0_ratio
  logical, optional :: invert
  real(rp) :: mat_out(size(mat_in, 1), size(mat_in, 2))
end function

function mat_symp_conj(mat) result (mat_conj)
  import
  implicit none
  real(rp) mat(:,:)
  real(rp) mat_conj(size(mat, 1), size(mat, 2))
end function

function mat_symp_conj_i(mat) result (mat_conj)
  import
  implicit none
  complex(rp) mat(:,:)
  complex(rp) mat_conj(size(mat, 1), size(mat, 2))
end function

function mat_symp_error (mat, p0_ratio, err_mat) result (error)
  import
  real(rp), intent(in) :: mat(:,:)
  real(rp), optional :: p0_ratio
  real(rp), optional :: err_mat(:,:)
  real(rp) error
end function

subroutine mat_symplectify (mat_in, mat_symp, p0_ratio, r_root)
  import
  real(rp), intent(in)  :: mat_in(:,:)
  real(rp), intent(out) :: mat_symp(:,:)
  real(rp), intent(in), optional :: p0_ratio, r_root
end subroutine

subroutine mat_type (mat, nunit, header, num_form, lines, n_lines)
  import
  implicit none
  real(rp) mat(:,:)
  integer, optional :: nunit, n_lines
  character(*), optional :: header, num_form, lines(:)
end subroutine

subroutine mat_make_unit (mat)
  import
  real(rp) mat(:,:)
end subroutine

subroutine match_word (string, names, ix, exact_case, can_abbreviate, matched_name)
  import
  implicit none
  character(*) string, names(:)
  character(*), optional :: matched_name
  integer ix
  logical, optional :: exact_case, can_abbreviate
end subroutine

function max_nonzero (lbnd, array1, array2) result (ix_max)
  import
  implicit none
  integer lbnd, ix_max
  real(rp) array1(lbnd:)
  real(rp), optional :: array2(lbnd:)
end function

function n_choose_k(n, k) result (nck)
  import
  implicit none
  real(rp) nck
  integer n, k
end function

subroutine n_spline_create (deriv0, deriv1, x1, n_spline)
  import
  implicit none
  real(rp) deriv0(0:), deriv1(0:)
  real(rp) x1
  real(rp) :: n_spline(0:)
end subroutine

subroutine nametable_add (nametable, name, ix_name)
  import
  implicit none
  type (nametable_struct), target :: nametable
  character(*) name
  integer ix_name
end subroutine

function nametable_bracket_indexx (nametable, name, n_match) result (ix_max)
  import
  implicit none
  type (nametable_struct) nametable
  character(*) name
  integer ix_max
  integer, optional :: n_match
end function

subroutine nametable_change1 (nametable, name, ix_name)
  import
  implicit none
  type (nametable_struct), target :: nametable
  character(*) name
  integer ix_name
end subroutine

subroutine nametable_init (nametable, n_min, n_max)
  import
  implicit none
  type (nametable_struct), target :: nametable
  integer, optional :: n_min, n_max
end subroutine

subroutine nametable_remove (nametable, ix_name)
  import
  implicit none
  type (nametable_struct), target :: nametable
  integer ix_name
end subroutine

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

function ordinal_str(n) result (str)
  implicit none
  integer n
  character(:), allocatable :: str
end function

function outer_product (a, b) result (c)
  import
  implicit none
  real(rp) a(:), b(:)
  real(rp) c(size(a), size(b))
end function

subroutine parse_fortran_format (format_str, n_repeat, power, descrip, width, digits)
  implicit none
  integer n_repeat, power, width, digits
  character(*) format_str, descrip
end subroutine

subroutine pointer_to_locations(string, array, num, ix_min, ix_max, names, exact_case, print_err)
  implicit none
  integer, allocatable :: array(:)
  integer num, ix_min, ix_max
  logical, optional :: exact_case, print_err
  character(*) string
  character(*), optional :: names(:)
end subroutine

function poly_eval(poly, x, diff_coef) result (y)
  import
  implicit none
  real(rp) poly(0:), x, y
  logical, optional :: diff_coef
end function

function probability_funct(x) result (prob)
  import
  implicit none
  real(rp) prob
  real(rp) x
end function

subroutine query_string (query_str, upcase, return_str, ix, ios)
  implicit none
  character(*) return_str
  character(*) query_str
  integer ix
  integer ios
  logical upcase
end subroutine

function quote(str) result (q_str)
  character(*) str
  character(:), allocatable :: q_str
end function

function quoten(str, delim) result (q_str)
  character(*) str(:)
  character(*), optional :: delim
  character(:), allocatable :: q_str
end function

function real_to_string (real_num, width, n_signif, n_decimal) result (str)
  import
  implicit none
  real(rp) real_num
  integer width
  integer, optional :: n_signif, n_decimal
  character(width) str
end function

function reals_to_string (real_arr, width, n_blank, n_signif, n_decimal) result (str)
  import
  implicit none
  real(rp) real_arr(:)
  integer width, n_blank
  integer, optional :: n_signif, n_decimal
  character(width*size(real_arr)+n_blank*(size(real_arr)-1)) str
end function

function reals_to_table_row (real_arr, width, n_decimal, n_blank) result (str)
  import
  implicit none
  real(rp) real_arr(:)
  integer width, n_decimal
  integer, optional :: n_blank
  character(width*size(real_arr)) str
end function

function all_pointer_to_string (a_ptr, err) result (str)
  import
  implicit none
  type (all_pointer_struct) a_ptr
  logical, optional :: err
  character(24) str
end function

function unquote (str_in) result (str_out)
  implicit none
  character(*) str_in
  character(len(str_in)) str_out
end function

function real_num_fortran_format (number, width, n_blanks) result (fmt_str)
  import
  implicit none
  real(rp) number
  integer, optional :: n_blanks
  integer width
  character(9) fmt_str
end function

subroutine str_set(str_out, str_in)
  implicit none
  character(:), allocatable :: str_out
  character(*) str_in
end subroutine

subroutine svd_fit (A, b, cutoff, x_min, chisq, w_vec, v_mat)
  import
  implicit none
  real(rp) A(:,:), b(:), cutoff, x_min(:)
  real(rp), optional :: chisq, w_vec(:), v_mat(:,:)
end subroutine

function real_path (path_in, path_out) result (is_ok)
  implicit none
  character(*) path_in, path_out
  logical is_ok
end function

function real_str(r_num, n_signif, n_decimal) result (str)
  import
  implicit none
  real(rp) r_num
  integer n_signif
  integer, optional :: n_decimal
  character(:), allocatable :: str
end function

function rms_value(val_arr, good_val, ave_val) result (rms_val)
  import
  implicit none
  real(rp) val_arr(:), rms_val
  logical, optional :: good_val(:)
  real(rp), optional :: ave_val
end function

function rot_2d (vec_in, angle) result (vec_out)
  import
  implicit none
  real(rp) vec_in(2), angle, vec_out(2)
end function

subroutine run_timer(command, time)
  import
  implicit none
  real(rp), optional :: time
  character(*) command
end subroutine

elemental function sinc(x, nd) result (y)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) y
  integer, optional, intent(in) :: nd
end function

elemental function sincc(x, nd) result (y)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) y
  integer, optional, intent(in) :: nd
end function

elemental function sinhx_x(x, nd) result (y)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) y
  integer, optional, intent(in) :: nd
end function

subroutine skip_header (ix_unit, error_flag)
  implicit none
  integer ix_unit
  logical error_flag
end subroutine

elemental function sqrt_one(x, nd) result (ds1)
  import
  real(rp), intent(in) :: x
  real(rp) ds1
  integer, optional, intent(in) :: nd
end function

function str_first_in_set(line, set, ignore_clauses) result (ix_match)
  implicit none
  character(*) line
  character(*) set
  logical, optional :: ignore_clauses
  integer ix_match
end function

function str_first_not_in_set(line, set) result (ix_match)
  implicit none
  character(*) line
  character(*) set
  integer ix_match
end function

function str_last_in_set(line, set) result (ix_match)
  implicit none
  character(*) line
  character(*) set
  integer ix_match
end function

function str_last_not_in_set(line, set) result (ix_match)
  implicit none
  character(*) line
  character(*) set
  integer ix_match
end function

function string_to_int (line, default, err_flag, err_print_flag) result (value)
  implicit none
  integer default
  integer value
  character(*) line
  logical err_flag
  logical, optional :: err_print_flag
end function

function string_to_real (line, default, err_flag, err_print_flag) result (value)
  import
  implicit none
  real(rp) default
  real(rp) value
  character(*) line
  logical err_flag
  logical, optional :: err_print_flag
end function

subroutine string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next)
  implicit none
  character(*) in_str
  character(*) out_str
  character(*) delimitors
  character(1) delim
  integer ix_word
  integer ix_next
end subroutine

function substr(var_str, n1, n2) result (sub_str)
  implicit none
  character(:), allocatable :: var_str
  integer n1, n2
  character(n2-n1+1) sub_str
end function

subroutine test_tune_tracker_lock (tracker_locked)
  implicit none
  logical tracker_locked(2)
end subroutine

function to_str(num, max_signif) result (string)
  import
  implicit none
  real(rp) num
  integer, optional :: max_signif
  character(:), allocatable :: string
end function

subroutine type_this_file(filename)
  implicit none
  character(*) filename
end subroutine

function upcase(str_in) result (str_out)
  implicit none
  character(*) str_in
  character(len(str_in)) str_out
end function

subroutine upcase_string(string)
  implicit none
  character(*) string
end subroutine

function word_len (wording) result (wlen)
  implicit none
  integer wlen 
  character(*) wording
end function

subroutine word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str, ignore_interior)
  implicit none
  character(*) in_str, out_str
  character(*) word
  character(*) delim_list, delim
  integer ix_word
  logical delim_found
  logical, optional :: ignore_interior
end subroutine

subroutine str_substitute (string, str_match, str_replace, do_trim, ignore_escaped)
  implicit none
  character(*) string
  character(*), optional :: str_match, str_replace
  logical, optional :: do_trim, ignore_escaped
end subroutine

recursive function str_match_wild(str, pat) result (a_match)
  implicit none
  character(*) pat, str
  logical a_match
end function str_match_wild

subroutine str_upcase(dst, src)
  implicit none
  character(*) dst, src
end subroutine str_upcase

subroutine str_downcase(dst, src)
  implicit none
  character(*) dst, src
end subroutine str_downcase

subroutine system_command (line, err_flag)
  implicit none
  character(*) line
  logical, optional :: err_flag
end subroutine

subroutine string_trim (in_string, out_string, word_len)
  implicit none
  character(*) in_string, out_string
  integer word_len
end subroutine string_trim

function virtual_memory_usage() result (usage)
  implicit none
  integer usage
end function

end interface

!

interface find_location
  function find_location_real(arr, value) result (ix_match)
    import
    real(rp) arr(:), value
    integer ix_match
  end function

  function find_location_int(arr, value) result (ix_match)
    integer arr(:), value
    integer ix_match
  end function

  function find_location_logic(arr, value) result (ix_match)
    logical arr(:), value
    integer ix_match
  end function

  function find_location_str(arr, value) result (ix_match)
    character(*) arr(:), value
    integer ix_match
  end function
end interface

! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
