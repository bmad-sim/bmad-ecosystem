module sim_utils_interface

use precision_def

interface
  subroutine bbi_kick (x, y, r, kx, ky)
    import
    implicit none
    real(rp) x
    real(rp) y
    real(rp) r
    real(rp) kx
    real(rp) ky
  end subroutine

  subroutine bracket_index (s_arr, i_min, i_max, s, ix)
    import
    implicit none
    integer i_min, i_max
    real(rp) s_arr(i_min:), s
    integer ix
  end subroutine

  subroutine bracket_index2 (s_arr, i_min, i_max, s, ix0, ix)
    import
    implicit none
    integer i_min, i_max
    real(rp) s_arr(i_min:), s
    integer ix0, ix
  end subroutine

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

  subroutine cplx_lubksb(a,indx,b)
    use nrtype
    import
    complex(dp), dimension(:,:), intent(in) :: a
    integer(I4B), dimension(:), intent(in) :: indx
    complex(dp), dimension(:), intent(inout) :: b
  end subroutine

  subroutine cplx_ludcmp(a,indx,d)
    use nrtype
    import
    complex(dp), dimension(:,:), intent(inout) :: a
    integer(I4B), dimension(:), intent(out) :: indx
    real(dp), intent(out) :: d
  end subroutine

  subroutine cplx_mat_inverse(mat_r, mat_i, inv_r, inv_i, ok, print_err)
    import
    real(rp) :: mat_r(:,:)
    real(rp) :: mat_i(:,:)
    real(rp) :: inv_r(:,:)
    real(rp) :: inv_i(:,:)
    logical, optional :: ok, print_err
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

  subroutine date_and_time_stamp (string, numeric_month)
    implicit none
    character(*) string
    logical, optional :: numeric_month
  end subroutine

  function determinant (mat) result (det)
    import
    implicit none
    real(rp) mat(:,:)
    real(rp) det
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

  function even (num) result (is_even)
    implicit none
    integer num
    logical is_even
  end function

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

  function is_integer (string) result (valid)
    character(*) string
    logical valid
  end function

  function is_logical (string, ignore) result (valid)
    character(*) string
    logical, optional :: ignore
    logical valid
  end function

  function is_real (string, ignore) result (valid)
    character(*) string
    logical, optional :: ignore
    logical valid
  end function

  subroutine linear_fit (x, y, n_data, a, b, sig_a, sig_b)
    import
    implicit none
    integer n_data
    real(rp) x(*)
    real(rp) y(*)
    real(rp) a
    real(rp) b
    real(rp) sig_a
    real(rp) sig_b
  end subroutine

  subroutine location_decode(string, array, ix_min, num, names, exact_case)
    implicit none
    integer num
    integer ix_min
    character(*) string
    logical array(ix_min:)
    character(*), optional :: names(ix_min:)
    logical, optional :: exact_case
  end subroutine

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

  subroutine mat_type (mat, nunit, header, num_form)
    import
    implicit none
    real(rp) mat(:,:)
    integer, optional :: nunit
    character(*), optional :: header, num_form
  end subroutine

  subroutine mat_make_unit (mat)
    import
    real(rp) mat(:,:)
  end subroutine

  function max_nonzero (lbnd, array1, array2) result (ix_max)
    import
    integer lbnd, ix_max
    real(rp) array1(lbnd:)
    real(rp), optional :: array2(lbnd:)
  end function

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

  function odd (num) result (is_odd)
    implicit none
    integer num
    logical is_odd
  end function

  function probability_funct(x) result (prob)
    import
    implicit none
    real(rp) prob
    real(rp) x
  end function

  subroutine ps2gif (ps_file, gif_file, kill_ps_file)
    implicit none
    character(*) ps_file, gif_file
    logical, optional :: kill_ps_file
  end subroutine

  subroutine query_string (query_str, upcase, return_str, ix, ios)
    implicit none
    character(*) return_str
    character(*) query_str
    integer ix
    integer ios
    logical upcase
  end subroutine

  function real_to_string (num, fmt) result (str)
    import
    implicit none
    real(rp) num
    character(*), optional :: fmt
    character(24) str
  end function

  function remove_quotes (str_in) result (str_out)
    implicit none
    character(*) str_in
    character(len(str_in)) str_out
  end function

  subroutine run_timer(command, time)
    import
    implicit none
    real(rp), optional :: time
    character(*) command
  end subroutine

  subroutine skip_header (ix_unit, error_flag)
    implicit none
    integer ix_unit
    logical error_flag
  end subroutine

  function str_find_first_in_set(line, set) result (ix_match)
    implicit none
    character(*) line
    character(*) set
    integer ix_match
  end function

  function str_find_first_not_in_set(line, set) result (ix_match)
    implicit none
    character(*) line
    character(*) set
    integer ix_match
  end function

  function str_find_last_in_set(line, set) result (ix_match)
    implicit none
    character(*) line
    character(*) set
    integer ix_match
  end function

  function str_find_last_not_in_set(line, set) result (ix_match)
    implicit none
    character(*) line
    character(*) set
    integer ix_match
  end function

  subroutine string_to_int (line, default, value, err_flag)
    implicit none
    integer default
    integer value
    character(*) line
    logical err_flag
  end subroutine

  subroutine string_to_real (line, default, value, err_flag)
    import
    implicit none
    real(rp) default
    real(rp) value
    character(*) line
    logical err_flag
  end subroutine

  subroutine string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next)
    implicit none
    character(*) in_str
    character(*) out_str
    character(*) delimitors
    character(1) delim
    integer ix_word
    integer ix_next
  end subroutine

  subroutine test_tune_tracker_lock (tracker_locked)
    implicit none
    logical tracker_locked(2)
  end subroutine

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

  subroutine word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str)
    implicit none
    character(*) in_str, out_str
    character(*) word
    character(*) delim_list, delim
    integer ix_word
    logical delim_found
  end subroutine

  subroutine str_substitute (string, str_match, str_replace, do_trim)
    implicit none
    character(*) string
    character(*), optional :: str_match, str_replace
    logical, optional :: do_trim
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

  subroutine system_command (line)
    implicit none
    character(*) line
  end subroutine

  subroutine string_trim (in_string, out_string, word_len)
    implicit none
    character(*) in_string, out_string
    integer word_len
  end subroutine string_trim

end interface

! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
