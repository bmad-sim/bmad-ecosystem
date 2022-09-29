! Get array of all nth order monomials.
! cycle over nonzero monomials    mad_tpsa_cycle    & mad_ctpsa_...
! map composition                 mad_tpsa_compose
! map inversion                   mad_tpsa_minv
! map partial inversion           mad_tpsa_pminv
! map eval                        mad_tpsa_eval
! map permutation                 mad_tpsa_translate

! Not interfaced:
!   mad_desc_newk
!   mad_desc_newv
!   mad_desc_newn
!   mad_desc_del

module gtpsa_mod

use gtpsa

implicit none

! Defs

integer, parameter :: dp = selected_real_kind(11)
real(c_num_t), parameter :: zero_cn = 0.0d0
real(c_num_t), parameter :: one_cn = 1.0d0
complex(c_cnum_t) :: zero_ccn = cmplx(0.0d0, 0.0d0)
complex(c_cnum_t) :: one_ccn = cmplx(1.0d0, 0.0d0)

! Classes...
! tpsa_class: A real valued taylor series

type tpsa_class
  type(c_ptr) :: ptr = c_null_ptr
contains
  final :: tpsa_final
  procedure :: coef => gtpsa_coef_t                ! Get monomial coef
  procedure :: n_var => gtpsa_n_var_t              ! Total number of variables including knobs
  procedure :: n_main_var => gtpsa_n_main_var_t    ! Total number of non-knob variables (typically the phase-space variables).
  procedure :: n_knobs => gtpsa_n_knobs_t          ! Total number of knob variables
  procedure :: order_max => gtpsa_order_max_t      ! Maximum possible order of any monomial. Order = sum:n_i, i = 1, n_var
  procedure :: order_knobs => gtpsa_order_knobs_t  ! Maximum possible order of knob variables. knob_order = sum:n_i, i = knob variable index.
  procedure :: order_vec => gtpsa_order_vec_t      ! Maximum possible order for a single variable n_i.
  procedure :: order_clear => gtpsa_order_clear_t  ! Clear (zero) monomials based on the monomial order.
  procedure :: n_mono => gtpsa_n_mono_t            ! Number of monomials.
  procedure :: mono_incr => gtpsa_mono_incr_t      ! Increment a monomial given signature.
  procedure :: mono_set => gtpsa_mono_set_t        ! Set a monomial given signature.
  procedure :: mono_setv => gtpsa_mono_setv_t      ! Set a vector of monomials from some starting index.
  procedure :: mono_get_by_index => gtpsa_mono_get_by_index_t   ! Get N^th monomial by index.
  procedure :: mono_get_by_expn => gtpsa_mono_get_by_expn_t     ! Get monomial with given signature.
  procedure :: swap_var_order => gtpsa_swap_var_order_t         ! Swap the order of the variables.
  procedure :: print => gtpsa_print_t              ! Print tpsa
  procedure :: delete => gtpsa_delete_t            ! Delete tpsa
end type

! map_class: An array of tpsas

type map_class
  type (tpsa_class), allocatable :: t(:)
end type

! ctpsa_class: A complex valued Taylor series.

type ctpsa_class
  type(c_ptr) :: ptr = c_null_ptr
contains
  final :: ctpsa_final
  procedure :: coef => gtpsa_coef_z                ! Get monomial coef
  procedure :: n_var => gtpsa_n_var_z              ! Total number of variables including knobs
  procedure :: n_main_var => gtpsa_n_main_var_z    ! Total number of non-knob variables (typically the phase-space variables).
  procedure :: n_knobs => gtpsa_n_knobs_z          ! Total number of knob variables
  procedure :: order_max => gtpsa_order_max_z      ! Maximum possible order of any monomial. Order = sum:n_i, i = 1, n_var
  procedure :: order_knobs => gtpsa_order_knobs_z  ! Maximum possible order of knob variables. knob_order = sum:n_i, i = knob variable index.
  procedure :: order_vec => gtpsa_order_vec_z      ! Maximum possible order for a single variable n_i.
  procedure :: order_clear => gtpsa_order_clear_z  ! Clear (zero) monomials based on the monomial order.
  procedure :: n_mono => gtpsa_n_mono_z            ! Number of monomials.
  procedure :: mono_incr => gtpsa_mono_incr_z      ! Increment a monomial.
  procedure :: mono_set => gtpsa_mono_set_z        ! Set a monomial.
  procedure :: mono_setv => gtpsa_mono_setv_z      ! Set a vector of monomials.
  procedure :: mono_get_by_index => gtpsa_mono_get_by_index_z   ! Get N^th monomial.
  procedure :: mono_get_by_expn => gtpsa_mono_get_by_expn_z     ! Get monomial with given signature.
  procedure :: swap_var_order => gtpsa_swap_var_order_z         ! Swap the order of the variables.
  procedure :: print => gtpsa_print_z              ! Print ctpsa.
  procedure :: delete => gtpsa_delete_z            ! Delete tpsa
  procedure :: real => gtpsa_real_z                ! Real part of ctpsa.
  procedure :: imag => gtpsa_imag_z                ! Imag part of ctpsa.
end type
 
! cmap_class: An array of ctpsas.

type cmap_class
  type (tpsa_class), allocatable :: c(:)
end type

! Routines that are not overloaded

! gtpsa_newd

! Misc

interface gtpsa_new
  module procedure gtpsa_new_t
  module procedure gtpsa_new_z
end interface

interface cmplx
  module procedure gtpsa_cmplx_t
  module procedure gtpsa_cmplx_tt
end interface

! Special functions

interface sin
  module procedure gtpsa_sin_t     ! sin(tpsa)
  module procedure gtpsa_sin_z     ! sin(cpsa)
end interface

interface cos
  module procedure gtpsa_cos_t     ! cos(tpsa)
  module procedure gtpsa_cos_z     ! cos(cpsa)
end interface

interface tan
  module procedure gtpsa_tan_t     ! tan(tpsa)
  module procedure gtpsa_tan_z     ! tan(cpsa)
end interface

interface cot
  module procedure gtpsa_cot_t     ! cot(tpsa)
  module procedure gtpsa_cot_z     ! cot(cpsa)
end interface

interface sinc
  module procedure gtpsa_sinc_t     ! sinc(tpsa)
  module procedure gtpsa_sinc_z     ! sinc(cpsa)
end interface

interface sinh
  module procedure gtpsa_sinh_t     ! sinh(tpsa)
  module procedure gtpsa_sinh_z     ! sinh(cpsa)
end interface

interface cosh
  module procedure gtpsa_cosh_t     ! cosh(tpsa)
  module procedure gtpsa_cosh_z     ! cosh(cpsa)
end interface

interface tanh
  module procedure gtpsa_tanh_t     ! tanh(tpsa)
  module procedure gtpsa_tanh_z     ! tanh(cpsa)
end interface

interface coth
  module procedure gtpsa_coth_t     ! coth(tpsa)
  module procedure gtpsa_coth_z     ! coth(cpsa)
end interface

interface sinhc
  module procedure gtpsa_sinhc_t     ! sinhc(tpsa)
  module procedure gtpsa_sinhc_z     ! sinhc(cpsa)
end interface

interface asin
  module procedure gtpsa_asin_t     ! asin(tpsa)
  module procedure gtpsa_asin_z     ! asin(cpsa)
end interface

interface acos
  module procedure gtpsa_acos_t     ! acos(tpsa)
  module procedure gtpsa_acos_z     ! acos(cpsa)
end interface

interface atan
  module procedure gtpsa_atan_t     ! atan(tpsa)
  module procedure gtpsa_atan_z     ! atan(cpsa)
end interface

interface acot
  module procedure gtpsa_acot_t     ! acot(tpsa)
  module procedure gtpsa_acot_z     ! acot(cpsa)
end interface

interface asinh
  module procedure gtpsa_asinh_t     ! asinh(tpsa)
  module procedure gtpsa_asinh_z     ! asinh(cpsa)
end interface

interface acosh
  module procedure gtpsa_acosh_t     ! acosh(tpsa)
  module procedure gtpsa_acosh_z     ! acosh(cpsa)
end interface

interface atanh
  module procedure gtpsa_atanh_t     ! atanh(tpsa)
  module procedure gtpsa_atanh_z     ! atanh(cpsa)
end interface

interface acoth
  module procedure gtpsa_acoth_t     ! acoth(tpsa)
  module procedure gtpsa_acoth_z     ! acoth(cpsa)
end interface

interface erf
  module procedure gtpsa_erf_t     ! erf(tpsa)
  module procedure gtpsa_erf_z     ! erf(cpsa)
end interface

interface erfc
  module procedure gtpsa_erfc_t     ! erfc(tpsa)
  module procedure gtpsa_erfc_z     ! erfc(cpsa)
end interface

interface sqrt
  module procedure gtpsa_sqrt_t     ! sqrt(tpsa)
  module procedure gtpsa_sqrt_z     ! sqrt(cpsa)
end interface

interface exp
  module procedure gtpsa_exp_t     ! exp(tpsa)
  module procedure gtpsa_exp_z     ! exp(cpsa)
end interface

interface log
  module procedure gtpsa_log_t     ! log(tpsa)
  module procedure gtpsa_log_z     ! log(cpsa)
end interface

interface log_x_div_y  ! log(x/y)
  module procedure gtpsa_log_x_div_y_tt     ! log(tpsa/tpsa)
  module procedure gtpsa_log_x_div_y_zz     ! log(ctpsa/ctpsa)
end interface

interface deriv  ! (tpsa, vec)
  module procedure gtpsa_deriv_t     ! deriv(tpsa, dvec)
  module procedure gtpsa_deriv_z     ! deriv(cpsa, dvec)
end interface

! Map methods

interface map_new
  module procedure map_new_m
  module procedure map_new_c
end interface

interface inv         ! inverse
  module procedure map_inv_m
  module procedure map_inv_c
end interface  

!interface pinv        ! partial inverse
!  module procedure map_pinv_m
!  module procedure map_pinv_c
!end interface  

!interface eval
!  module procedure map_eval_m
!  module procedure map_eval_c
!end interface  

! binary operators
! Rule: Any computation with tpsa and ctpsa will result in a ctpsa.
! Rule: Any computation with tpsa and complex will result in a ctpsa.

interface poisson
  module procedure gtpsa_poisson_tt
  module procedure gtpsa_poisson_zz
end interface

!interface operator(.comp.)     ! Composition
!  module procedure gtpsa_comp_tt
!  module procedure gtpsa_comp_zz
!end interface

interface assignment(=)
  module procedure gtpsa_equal_tt  ! tpsa = tpsa
  module procedure gtpsa_equal_tz  ! tpsa = ctpsa
  module procedure gtpsa_equal_zt  ! ctpsa = tpsa
  module procedure gtpsa_equal_zz  ! ctpsa = ctpsa

  module procedure map_equal_mm  ! map = map
  module procedure map_equal_mc  ! map = cmap
  module procedure map_equal_cm  ! cmap = map
  module procedure map_equal_cc  ! cmap = cmap
end interface

interface operator(+)
  module procedure gtpsa_add_tt    ! tpsa + tpsa
  module procedure gtpsa_add_tz    ! tpsa + ctpsa
  module procedure gtpsa_add_zt    ! ctpsa + tpsa
  module procedure gtpsa_add_zz    ! ctpsa + ctpsa

  module procedure gtpsa_add_tc    ! tpsa + complex
  module procedure gtpsa_add_ct    ! complex + tpsa
  module procedure gtpsa_add_tr    ! tpsa + real
  module procedure gtpsa_add_rt    ! real + tpsa

  module procedure gtpsa_add_zc    ! ctpsa + complex
  module procedure gtpsa_add_cz    ! complex + ctpsa
  module procedure gtpsa_add_zr    ! ctpsa + real
  module procedure gtpsa_add_rz    ! real + ctpsa
end interface

interface operator(-)
  module procedure gtpsa_subtract_tt    ! tpsa - tpsa
  module procedure gtpsa_subtract_tz    ! tpsa - ctpsa
  module procedure gtpsa_subtract_zt    ! ctpsa - tpsa
  module procedure gtpsa_subtract_zz    ! ctpsa - ctpsa

  module procedure gtpsa_subtract_tc    ! tpsa - complex
  module procedure gtpsa_subtract_ct    ! complex - tpsa
  module procedure gtpsa_subtract_tr    ! tpsa - real
  module procedure gtpsa_subtract_rt    ! real - tpsa

  module procedure gtpsa_subtract_zc    ! ctpsa - complex
  module procedure gtpsa_subtract_cz    ! complex - ctpsa
  module procedure gtpsa_subtract_zr    ! ctpsa - real
  module procedure gtpsa_subtract_rz    ! real - ctpsa
end interface

interface operator(*)
  module procedure gtpsa_multiply_tt    ! tpsa * tpsa
  module procedure gtpsa_multiply_tz    ! tpsa * ctpsa
  module procedure gtpsa_multiply_zt    ! ctpsa * tpsa
  module procedure gtpsa_multiply_zz    ! ctpsa * ctpsa

  module procedure gtpsa_multiply_tc    ! tpsa * complex
  module procedure gtpsa_multiply_ct    ! complex * tpsa
  module procedure gtpsa_multiply_tr    ! tpsa * real
  module procedure gtpsa_multiply_rt    ! real * tpsa

  module procedure gtpsa_multiply_zc    ! ctpsa * complex
  module procedure gtpsa_multiply_cz    ! complex * ctpsa
  module procedure gtpsa_multiply_zr    ! ctpsa * real
  module procedure gtpsa_multiply_rz    ! real * ctpsa
end interface

interface operator(/)
  module procedure gtpsa_divide_tt    ! tpsa / tpsa
  module procedure gtpsa_divide_tz    ! tpsa / ctpsa
  module procedure gtpsa_divide_zt    ! ctpsa / tpsa
  module procedure gtpsa_divide_zz    ! ctpsa / ctpsa

  module procedure gtpsa_divide_tc    ! tpsa / complex
  module procedure gtpsa_divide_ct    ! complex / tpsa
  module procedure gtpsa_divide_tr    ! tpsa / real
  module procedure gtpsa_divide_rt    ! real / tpsa

  module procedure gtpsa_divide_zc    ! ctpsa / cmplx
  module procedure gtpsa_divide_cz    ! cmplx / ctpsa
  module procedure gtpsa_divide_zr    ! ctpsa / real
  module procedure gtpsa_divide_rz    ! real / ctpsa
end interface

interface operator(**)
  module procedure gtpsa_power_tt    ! tpsa^tpsa
  module procedure gtpsa_power_ti    ! tpsa^integer
  module procedure gtpsa_power_tr    ! ctpsa^real
  module procedure gtpsa_power_zz    ! ctpsa^ctpsa
  module procedure gtpsa_power_zi    ! ctpsa^integer
  module procedure gtpsa_power_zr    ! ctpsa^real
  module procedure gtpsa_power_zc    ! ctpsa^complex
end interface

contains

!------------------------------------------------------------------------------------------
! Finalizers

subroutine tpsa_final(t)
  type(tpsa_class) :: t
  call mad_tpsa_del(t%ptr)
  t%ptr = c_null  ! Just to make sure.
end subroutine tpsa_final

subroutine ctpsa_final(t)
  type(ctpsa_class) :: t
  call mad_ctpsa_del(t%ptr)
  t%ptr = c_null  ! Just to make sure.
end subroutine ctpsa_final

!------------------------------------------------------------------------------------------
! gtpsa_newd

function gtpsa_newd (desc, mo) result (t_new)
  type(c_ptr), value, intent(in) :: desc
  integer (c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
  type(tpsa_class) t_new
  t_new%ptr = mad_tpsa_newd(desc, mo)
end function gtpsa_newd

!------------------------------------------------------------------------------------------
! interface procedure: gtpsa_new

function gtpsa_new_t (t1, mo) result (t_out)
  type(tpsa_class), value, intent(in) :: t1
  integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mo)
end function gtpsa_new_t

function gtpsa_new_z (z1, mo) result (z_out)
  type(ctpsa_class), value, intent(in) :: z1
  integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mo)
end function gtpsa_new_z

!------------------------------------------------------------------------------------------
! interface procedure: map_new

function map_new_m (m1, m_ord) result (m_out)
  type(map_class), value, intent(in) :: m1
  integer(c_ord_t), value, intent(in) :: m_ord
  type(map_class) m_out
  integer n
  if (allocated(m_out%t)) then
    if (size(m_out%t) /= size(m1%t)) deallocate(m_out%t)
  endif
  allocate(m_out%t(size(m1%t)))
  do n = 1, size(m1%t)
    m_out%t(n)%ptr = mad_tpsa_new(m1%t(n)%ptr, m_ord)
  enddo
end function map_new_m

function map_new_c (c1, m_ord) result (c_out)
  type(cmap_class), value, intent(in) :: c1
  integer(c_ord_t), value, intent(in) :: m_ord
  type(cmap_class) c_out
  integer n
  if (allocated(c_out%c)) then
    if (size(c_out%c) /= size(c1%c)) deallocate(c_out%c)
  endif
  allocate(c_out%c(size(c1%c)))
  do n = 1, size(c1%c)
    c_out%c(n)%ptr = mad_ctpsa_new(c1%c(n)%ptr, m_ord)
  enddo
end function map_new_c

!------------------------------------------------------------------------------------------
! interface inv



!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: coef

function gtpsa_coef_t (t1, evec) result (r_out)
  class(tpsa_class), intent(in) :: t1
  integer, intent(in) :: evec(:)
  real(dp) r_out
  r_out = mad_tpsa_getm(t1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t))
end function gtpsa_coef_t

function gtpsa_coef_z (z1, evec) result (c_out)
  class(ctpsa_class), intent(in) :: z1
  integer, intent(in) :: evec(:)
  complex(dp) c_out
  c_out = mad_ctpsa_getm(z1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t))
end function gtpsa_coef_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: n_var

function gtpsa_n_var_t (t1) result (nv)
  class(tpsa_class), intent(in) :: t1
  integer nv
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr) 
  nv = mad_desc_nvmok(desc)
end function gtpsa_n_var_t

function gtpsa_n_var_z (z1) result (nv)
  class(ctpsa_class), intent(in) :: z1
  integer nv
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr) 
  nv = mad_desc_nvmok(desc)
end function gtpsa_n_var_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: n_main_var

function gtpsa_n_main_var_t (t1) result (nmv)
  class(tpsa_class), intent(in) :: t1
  integer nmv
  integer(c_int) nk
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr) 
  nmv = mad_desc_nvmok(desc, nk_ = nk)
  nmv = nmv - nk
end function gtpsa_n_main_var_t

function gtpsa_n_main_var_z (z1) result (nmv)
  class(ctpsa_class), intent(in) :: z1
  integer nmv
  integer(c_int) nk
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr) 
  nmv = mad_desc_nvmok(desc, nk_ = nk)
  nmv = nmv - nk
end function gtpsa_n_main_var_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: n_knobs

function gtpsa_n_knobs_t (t1) result (nk)
  class(tpsa_class), intent(in) :: t1
  integer nk, nv
  integer(c_int) nkn
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr)
  nv = mad_desc_nvmok(desc, nk_ = nkn)
  nk = nkn
end function gtpsa_n_knobs_t

function gtpsa_n_knobs_z (z1) result (nk)
  class(ctpsa_class), intent(in) :: z1
  integer nk, nv
  integer(c_int) nkn
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr)
  nv = mad_desc_nvmok(desc, nk_ = nkn)
  nk = nkn
end function gtpsa_n_knobs_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: order_max

function gtpsa_order_max_t (t1) result (ord)
  class(tpsa_class), intent(in) :: t1
  integer ord, nv
  integer(c_ord_t) m_ord
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr)
  nv = mad_desc_nvmok(desc, mo_ = m_ord)
  ord = m_ord
end function gtpsa_order_max_t

function gtpsa_order_max_z (z1) result (ord)
  class(ctpsa_class), intent(in) :: z1
  integer ord, nv
  integer(c_ord_t) m_ord
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr)
  nv = mad_desc_nvmok(desc, mo_ = m_ord)
  ord = m_ord
end function gtpsa_order_max_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: order_knobs

function gtpsa_order_knobs_t (t1) result (ordk)
  class(tpsa_class), intent(in) :: t1
  integer ordk, nv
  integer(c_ord_t) ko
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr)
  nv = mad_desc_nvmok(desc, ko_ = ko)
  ordk = ko
end function gtpsa_order_knobs_t

function gtpsa_order_knobs_z (z1) result (ordk)
  class(ctpsa_class), intent(in) :: z1
  integer ordk, nv
  integer(c_ord_t) ko
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr)
  nv = mad_desc_nvmok(desc, ko_ = ko)
  ordk = ko
end function gtpsa_order_knobs_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: order_vec
! It is OK if size(ovec(:)) /= num variables.
! If larger, array will be padded with zeros.
! If smaller, order info for variables with indexes above the array size will be discarded.

subroutine gtpsa_order_vec_t (t1, ovec)
  class(tpsa_class), intent(in) :: t1
  integer ovec(:)
  integer(c_int) nv
  integer(c_ord_t) :: vo(size(ovec))
  type(c_ptr) desc
  desc = mad_tpsa_desc(t1%ptr)
  nv = min(mad_desc_nvmok(desc), size(ovec))
  nv = mad_desc_getvo(desc, nv, vo(1:nv))
  ovec = 0
  ovec(1:nv) = vo(1:nv)
end subroutine gtpsa_order_vec_t

subroutine gtpsa_order_vec_z (z1, ovec)
  class(ctpsa_class), intent(in) :: z1
  integer ovec(:)
  integer(c_int) nv
  integer(c_ord_t) :: vo(size(ovec))
  type(c_ptr) desc
  desc = mad_ctpsa_desc(z1%ptr)
  nv = min(mad_desc_nvmok(desc), size(ovec))
  nv = mad_desc_getvo(desc, nv, vo(1:nv))
  ovec = 0
  ovec(1:nv) = vo(1:nv)
end subroutine gtpsa_order_vec_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: order_clear
! If n1 is negative then clear from |n1| to zero.

subroutine gtpsa_order_clear_t (t1, n1)
  class(tpsa_class) t1
  integer n1
  call mad_tpsa_cutord(t1%ptr, t1%ptr, int(n1, c_int))
end subroutine gtpsa_order_clear_t

subroutine gtpsa_order_clear_z (z1, n1)
  class(ctpsa_class) z1
  integer n1
  call mad_ctpsa_cutord(z1%ptr, z1%ptr, int(n1, c_int))
end subroutine gtpsa_order_clear_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: n_mono
! Number of possible monomials.

function gtpsa_n_mono_t (t1) result (n)
  class(tpsa_class) :: t1
  integer n
  n = mad_tpsa_len(t1%ptr)
end function gtpsa_n_mono_t

function gtpsa_n_mono_z (z1) result (n)
  class(ctpsa_class) :: z1
  integer n
  n = mad_ctpsa_len(z1%ptr)
end function gtpsa_n_mono_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: mono_incr

subroutine gtpsa_mono_incr_t (t1, evec, r)
  class(tpsa_class) :: t1
  integer, intent(in) :: evec(:)
  real(dp), intent(in) :: r
  call mad_tpsa_setm(t1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t), one_cn, real(r, c_num_t))
end subroutine gtpsa_mono_incr_t

subroutine gtpsa_mono_incr_z (z1, evec, c)
  class(ctpsa_class) :: z1
  integer, intent(in) :: evec(:)
  complex(dp), intent(in) :: c
  call mad_ctpsa_setm(z1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t), one_ccn, cmplx(c, kind = c_cnum_t))
end subroutine gtpsa_mono_incr_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: mono_set

subroutine gtpsa_mono_set_t (t1, evec, r)
  class(tpsa_class) :: t1
  integer, intent(in) :: evec(:)
  real(dp), intent(in) :: r
  call mad_tpsa_setm(t1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t), zero_cn, real(r, c_num_t))
end subroutine gtpsa_mono_set_t

subroutine gtpsa_mono_set_z (z1, evec, c)
  class(ctpsa_class) :: z1
  integer, intent(in) :: evec(:)
  complex(dp), intent(in) :: c
  call mad_ctpsa_setm(z1%ptr, int(size(evec), c_ssz_t), int(evec, c_ord_t), zero_ccn, cmplx(c, kind = c_cnum_t))
end subroutine gtpsa_mono_set_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: mono_setv

subroutine gtpsa_mono_setv_t (t1, i, n, v)
  class(tpsa_class) :: t1
  integer(c_idx_t), intent(in) :: i   ! slot index (must be valid)
  integer(c_ssz_t), intent(in) :: n   ! vector length
  real(c_num_t), intent(in) :: v(*)          ! vector to copy
  call mad_tpsa_setv(t1%ptr, i, n, v)
end subroutine gtpsa_mono_setv_t

subroutine gtpsa_mono_setv_z (z1, i, n, v)
  class(ctpsa_class) :: z1
  integer(c_idx_t), intent(in) :: i   ! slot index (must be valid)
  integer(c_ssz_t), intent(in) :: n   ! vector length
  complex(c_cnum_t), intent(in) :: v(*)          ! vector to copy
  call mad_ctpsa_setv(z1%ptr, i, n, v)
end subroutine gtpsa_mono_setv_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: mono_by_index

function gtpsa_mono_get_by_index_t (t1, i) result (r)
  class(tpsa_class) :: t1
  integer(c_idx_t), intent(in) :: i   ! slot index (must be valid)
  real(dp) r
  r = mad_tpsa_geti(t1%ptr, i)
end function gtpsa_mono_get_by_index_t

function gtpsa_mono_get_by_index_z (z1, i) result (c)
  class(ctpsa_class) :: z1
  integer(c_idx_t), intent(in) :: i   ! slot index (must be valid)
  complex(dp) c
  c = mad_ctpsa_geti(z1%ptr, i)
end function gtpsa_mono_get_by_index_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: mono_by_expn
! If expn argument is not present, get coef of scalar monomial (the monomial with all exponents zero).

function gtpsa_mono_get_by_expn_t (t1, expn) result (r)
  class(tpsa_class) :: t1
  integer, optional :: expn(:)
  real(dp) r
  if (present(expn)) then
    r = mad_tpsa_getm(t1%ptr, int(size(expn), c_ssz_t), int(expn, c_ord_t))
  else
    r = mad_tpsa_getm(t1%ptr, int(0, c_ssz_t), [int(0, c_ord_t)])
  endif
end function gtpsa_mono_get_by_expn_t

function gtpsa_mono_get_by_expn_z (z1, expn) result (c)
  class(ctpsa_class) :: z1
  integer, optional :: expn(:)
  complex(dp) c
  if (present(expn)) then
    c = mad_ctpsa_getm(z1%ptr, int(size(expn), c_ssz_t), int(expn, c_ord_t))
  else
    c = mad_ctpsa_getm(z1%ptr, int(0, c_ssz_t), [int(0, c_ord_t)])
  endif
end function gtpsa_mono_get_by_expn_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: swap_var_order
! new_order argument is a vector of the new order of the variables.
! Example: new_order = [1, 2, 4, 3] -> swap variables V4 and V3.
! If any value of new_order is 0, zero any monomial coefs where corresponding exponent for given var is non-zero.
! Example: new_order = [1, 2, 0, 0] -> zero all monomials where V3 or V4 exponent is non-zero.

subroutine gtpsa_swap_var_order_t (t1, new_order)
  class(tpsa_class) t1
  integer new_order(:)
  call mad_tpsa_convert(t1%ptr, t1%ptr, int(size(new_order), c_ssz_t), int(new_order, c_idx_t), 0)
end subroutine gtpsa_swap_var_order_t

subroutine gtpsa_swap_var_order_z (z1, new_order)
  class(ctpsa_class) z1
  integer new_order(:)
  call mad_ctpsa_convert(z1%ptr, z1%ptr, int(size(new_order), c_ssz_t), int(new_order, c_idx_t), 0)
end subroutine gtpsa_swap_var_order_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: print

subroutine gtpsa_print_t (t1, name, eps_, nohdr_, stream_)
  class(tpsa_class), intent(in) :: t1    ! src
  character(*), intent(in) :: name          ! name 
  real(c_num_t), intent(in) :: eps_  ! display precision, e.g. 1d-12
  integer, intent(in) :: nohdr_      ! discard header if not zero
  type(c_ptr) :: stream_             ! dst=c_null_ptr => stdio
  call mad_tpsa_print(t1%ptr, trim(name)//c_eos, eps_, int(nohdr_, c_int), stream_)
end subroutine gtpsa_print_t

subroutine gtpsa_print_z (z1, name, eps_, nohdr_, stream_)
  class(ctpsa_class), intent(in) :: z1    ! src
  character(*), intent(in) :: name          ! name (i.e. null terminated string)
  real(c_num_t), intent(in) :: eps_  ! display precision, e.g. 1d-12
  integer, intent(in) :: nohdr_      ! discard header if not zero
  type(c_ptr) :: stream_             ! dst=c_null_ptr => stdio
  call mad_ctpsa_print(z1%ptr, trim(name)//c_eos, eps_, int(nohdr_, c_int), stream_)
end subroutine gtpsa_print_z

!------------------------------------------------------------------------------------------
! (c)tpsa_class procedure: delete

subroutine gtpsa_delete_t (t1)
  class(tpsa_class), intent(in) :: t1
  call mad_tpsa_del(t1%ptr)
end subroutine gtpsa_delete_t

subroutine gtpsa_delete_z (z1)
  class(ctpsa_class), intent(in) :: z1
  call mad_ctpsa_del(z1%ptr)
end subroutine gtpsa_delete_z

!------------------------------------------------------------------------------------------
! interface: real

function gtpsa_real_z (z1) result (t_out)
  class (ctpsa_class), intent(in) :: z1
  type (tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_real (z1%ptr, t_out%ptr)
end function gtpsa_real_z

!------------------------------------------------------------------------------------------
! interface: imag

function gtpsa_imag_z (z1) result (t_out)
  class (ctpsa_class), intent(in) :: z1
  type (tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_imag (z1%ptr, t_out%ptr)
end function gtpsa_imag_z

!------------------------------------------------------------------------------------------
! interface: cmplx

function gtpsa_cmplx_t (t1) result (z_out)
  type (tpsa_class), intent(in) :: t1
  type (ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr)
end function gtpsa_cmplx_t

function gtpsa_cmplx_tt (t1, t2) result (z_out)
  type (tpsa_class), intent(in) :: t1, t2 
  type (ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, t2%ptr, z_out%ptr)
end function gtpsa_cmplx_tt

!------------------------------------------------------------------------------------------
! cmplx_conj

function cmplx_conj (z1) result (z_out)
  type (ctpsa_class), intent(in) :: z1
  type (ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (z1%ptr, c_null, z_out%ptr)
end function cmplx_conj

!------------------------------------------------------------------------------------------
! interface: inverse

function map_inv_m (m1) result (m_out)
  type (map_class), intent(in) :: m1
  type (map_class) m_out
!  m_out = map_new_m(m1, mad_tpsa_same)
!  call mad_tpsa_inv (m1%t, 1.0_dp, m_out%t)
end function map_inv_m

function map_inv_c (c1) result (c_out)
  type (cmap_class), intent(in) :: c1
  type (cmap_class) c_out
!  c_out = map_new_c(c1, mad_tpsa_same)
!  call mad_tpsa_inv (c1%z, 1.0_dp, c_out%z)
end function map_inv_c

!------------------------------------------------------------------------------------------
! interface: poisson

function gtpsa_poisson_tt (t1, t2) result (t_out)
  type (tpsa_class), intent(in) :: t1, t2 
  type (tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_poisson (t1%ptr, t2%ptr, t_out%ptr, 0_c_int)
end function gtpsa_poisson_tt

function gtpsa_poisson_zz (z1, z2) result (z_out)
  type (ctpsa_class), intent(in) :: z1, z2 
  type (ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_poisson (z1%ptr, z2%ptr, z_out%ptr, 0_c_int)
end function gtpsa_poisson_zz

!------------------------------------------------------------------------------------------
! interface: sin

function gtpsa_sin_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sin(t1%ptr, t_out%ptr)
end function gtpsa_sin_t

function gtpsa_sin_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sin(z1%ptr, z_out%ptr)
end function gtpsa_sin_z

!------------------------------------------------------------------------------------------
! interface: cos

function gtpsa_cos_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_cos(t1%ptr, t_out%ptr)
end function gtpsa_cos_t

function gtpsa_cos_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_cos(z1%ptr, z_out%ptr)
end function gtpsa_cos_z

!------------------------------------------------------------------------------------------
! interface: tan

function gtpsa_tan_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_tan(t1%ptr, t_out%ptr)
end function gtpsa_tan_t

function gtpsa_tan_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_tan(z1%ptr, z_out%ptr)
end function gtpsa_tan_z

!------------------------------------------------------------------------------------------
! interface: cot

function gtpsa_cot_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_cot(t1%ptr, t_out%ptr)
end function gtpsa_cot_t

function gtpsa_cot_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_cot(z1%ptr, z_out%ptr)
end function gtpsa_cot_z

!------------------------------------------------------------------------------------------
! interface: sinc

function gtpsa_sinc_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sinc(t1%ptr, t_out%ptr)
end function gtpsa_sinc_t

function gtpsa_sinc_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sinc(z1%ptr, z_out%ptr)
end function gtpsa_sinc_z

!------------------------------------------------------------------------------------------
! interface: sinh

function gtpsa_sinh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sinh(t1%ptr, t_out%ptr)
end function gtpsa_sinh_t

function gtpsa_sinh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sinh(z1%ptr, z_out%ptr)
end function gtpsa_sinh_z

!------------------------------------------------------------------------------------------
! interface: cosh

function gtpsa_cosh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_cosh(t1%ptr, t_out%ptr)
end function gtpsa_cosh_t

function gtpsa_cosh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_cosh(z1%ptr, z_out%ptr)
end function gtpsa_cosh_z

!------------------------------------------------------------------------------------------
! interface: tanh

function gtpsa_tanh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_tanh(t1%ptr, t_out%ptr)
end function gtpsa_tanh_t

function gtpsa_tanh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_tanh(z1%ptr, z_out%ptr)
end function gtpsa_tanh_z

!------------------------------------------------------------------------------------------
! interface: coth

function gtpsa_coth_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_coth(t1%ptr, t_out%ptr)
end function gtpsa_coth_t

function gtpsa_coth_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_coth(z1%ptr, z_out%ptr)
end function gtpsa_coth_z

!------------------------------------------------------------------------------------------
! interface: sinhc

function gtpsa_sinhc_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sinhc(t1%ptr, t_out%ptr)
end function gtpsa_sinhc_t

function gtpsa_sinhc_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sinhc(z1%ptr, z_out%ptr)
end function gtpsa_sinhc_z

!------------------------------------------------------------------------------------------
! interface: asin

function gtpsa_asin_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_asin(t1%ptr, t_out%ptr)
end function gtpsa_asin_t

function gtpsa_asin_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_asin(z1%ptr, z_out%ptr)
end function gtpsa_asin_z

!------------------------------------------------------------------------------------------
! interface: acos

function gtpsa_acos_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_acos(t1%ptr, t_out%ptr)
end function gtpsa_acos_t

function gtpsa_acos_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_acos(z1%ptr, z_out%ptr)
end function gtpsa_acos_z

!------------------------------------------------------------------------------------------
! interface: atan

function gtpsa_atan_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_atan(t1%ptr, t_out%ptr)
end function gtpsa_atan_t

function gtpsa_atan_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_atan(z1%ptr, z_out%ptr)
end function gtpsa_atan_z

!------------------------------------------------------------------------------------------
! interface: acot

function gtpsa_acot_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_acot(t1%ptr, t_out%ptr)
end function gtpsa_acot_t

function gtpsa_acot_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_acot(z1%ptr, z_out%ptr)
end function gtpsa_acot_z

!------------------------------------------------------------------------------------------
! interface: asinh

function gtpsa_asinh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_asinh(t1%ptr, t_out%ptr)
end function gtpsa_asinh_t

function gtpsa_asinh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_asinh(z1%ptr, z_out%ptr)
end function gtpsa_asinh_z

!------------------------------------------------------------------------------------------
! interface: acosh

function gtpsa_acosh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_acosh(t1%ptr, t_out%ptr)
end function gtpsa_acosh_t

function gtpsa_acosh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_acosh(z1%ptr, z_out%ptr)
end function gtpsa_acosh_z

!------------------------------------------------------------------------------------------
! interface: atanh

function gtpsa_atanh_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_atanh(t1%ptr, t_out%ptr)
end function gtpsa_atanh_t

function gtpsa_atanh_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_atanh(z1%ptr, z_out%ptr)
end function gtpsa_atanh_z

!------------------------------------------------------------------------------------------
! interface: acoth

function gtpsa_acoth_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_acoth(t1%ptr, t_out%ptr)
end function gtpsa_acoth_t

function gtpsa_acoth_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_acoth(z1%ptr, z_out%ptr)
end function gtpsa_acoth_z

!------------------------------------------------------------------------------------------
! interface: erf

function gtpsa_erf_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_erf(t1%ptr, t_out%ptr)
end function gtpsa_erf_t

function gtpsa_erf_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_erf(z1%ptr, z_out%ptr)
end function gtpsa_erf_z

!------------------------------------------------------------------------------------------
! interface: erfc

function gtpsa_erfc_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_erfc(t1%ptr, t_out%ptr)
end function gtpsa_erfc_t

function gtpsa_erfc_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_erfc(z1%ptr, z_out%ptr)
end function gtpsa_erfc_z

!------------------------------------------------------------------------------------------
! interface: sqrt

function gtpsa_sqrt_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sqrt(t1%ptr, t_out%ptr)
end function gtpsa_sqrt_t

function gtpsa_sqrt_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sqrt(z1%ptr, z_out%ptr)
end function gtpsa_sqrt_z

!------------------------------------------------------------------------------------------
! interface: exp

function gtpsa_exp_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_exp(t1%ptr, t_out%ptr)
end function gtpsa_exp_t

function gtpsa_exp_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_exp(z1%ptr, z_out%ptr)
end function gtpsa_exp_z

!------------------------------------------------------------------------------------------
! interface: log

function gtpsa_log_t (t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_log(t1%ptr, t_out%ptr)
end function gtpsa_log_t

function gtpsa_log_z (z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_log(z1%ptr, z_out%ptr)
end function gtpsa_log_z

!------------------------------------------------------------------------------------------
! interface: log_x_div_y

function gtpsa_log_x_div_y_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_logxdy(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_log_x_div_y_tt

function gtpsa_log_x_div_y_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_logxdy(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_log_x_div_y_zz

!------------------------------------------------------------------------------------------
! interface: deriv

function gtpsa_deriv_t (t1, dvec) result (t_out)
  type(tpsa_class), intent(in) :: t1
  integer, intent(in) :: dvec(:)
  type(tpsa_class) t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_derivm(t1%ptr, t_out%ptr, size(dvec), int(dvec, c_ord_t))
end function gtpsa_deriv_t

function gtpsa_deriv_z (z1, dvec) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  integer, intent(in) :: dvec(:)
  type(ctpsa_class) z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_derivm(z1%ptr, z_out%ptr, size(dvec), int(dvec, c_ord_t))
end function gtpsa_deriv_z

!------------------------------------------------------------------------------------------
! interface: assignment(=)

subroutine gtpsa_equal_tt (t_out, t_in)
  type(tpsa_class), intent(in) :: t_in
  type(tpsa_class), intent(inout) :: t_out
  t_out%ptr = mad_tpsa_new(t_in%ptr, mad_tpsa_same)
  call mad_tpsa_copy(t_in%ptr, t_out%ptr)
end subroutine gtpsa_equal_tt

subroutine gtpsa_equal_tz (t_out, z_in)
  type(ctpsa_class), intent(in) :: z_in
  type(tpsa_class), intent(inout) :: t_out
  t_out%ptr = mad_tpsa_new(z_in%ptr, mad_tpsa_same)
  call mad_ctpsa_real(z_in%ptr, t_out%ptr)
end subroutine gtpsa_equal_tz

subroutine gtpsa_equal_zt (z_out, t_in)
  type(tpsa_class), intent(in) :: t_in
  type(ctpsa_class), intent(inout) :: z_out
  z_out%ptr = mad_ctpsa_new(t_in%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx(t_in%ptr, c_null, z_out%ptr)
end subroutine gtpsa_equal_zt

subroutine gtpsa_equal_zz (z_out, z_in)
  type(ctpsa_class), intent(in) :: z_in
  type(ctpsa_class), intent(inout) :: z_out
  z_out%ptr = mad_ctpsa_new(z_in%ptr, mad_tpsa_same)
  call mad_ctpsa_copy(z_in%ptr, z_out%ptr)
end subroutine gtpsa_equal_zz

!------------------------------------------------------------------------------------------
! interface: assignment(=)

subroutine map_equal_mm (m_out, m_in)
  type(map_class), intent(in) :: m_in
  type(map_class), intent(inout) :: m_out
  m_out%t = m_in%t
end subroutine map_equal_mm

subroutine map_equal_mc (m_out, c_in)
  type(cmap_class), intent(in) :: c_in
  type(map_class), intent(inout) :: m_out
  integer n
  do n = 1, size(c_in%c)
    m_out%t(n) = c_in%c(n)
  enddo
end subroutine map_equal_mc

subroutine map_equal_cm (c_out, m_in)
  type(map_class), intent(in) :: m_in
  type(cmap_class), intent(inout) :: c_out
  integer n
  do n = 1, size(m_in%t)
    c_out%c(n) = m_in%t(n)
  enddo
end subroutine map_equal_cm

subroutine map_equal_cc (c_out, c_in)
  type(cmap_class), intent(in) :: c_in
  type(cmap_class), intent(inout) :: c_out
  c_out%c = c_in%c
end subroutine map_equal_cc

!------------------------------------------------------------------------------------------
! interface: operator(+)

function gtpsa_add_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_add(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_add_tt

function gtpsa_add_tz (t1, z2) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class), intent(in) :: z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_addt(z2%ptr, t1%ptr, z_out%ptr)
end function gtpsa_add_tz

function gtpsa_add_zt (z1, t2) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(tpsa_class), intent(in) :: t2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_addt(z1%ptr, t2%ptr, z_out%ptr)
end function gtpsa_add_zt

function gtpsa_add_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_add(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_add_zz

function gtpsa_add_tr (t1, r) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_copy (t1%ptr, t_out%ptr) 
  call mad_tpsa_set0(t_out%ptr, one_cn, real(r, c_num_t))
end function gtpsa_add_tr

function gtpsa_add_rt (r, t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_copy (t1%ptr, t_out%ptr) 
  call mad_tpsa_set0(t_out%ptr, one_cn, real(r, c_num_t))
end function gtpsa_add_rt

function gtpsa_add_tc (t1, c) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_add_tc

function gtpsa_add_ct (c, t1) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_add_ct

function gtpsa_add_zr (z1, r) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(r, kind=c_cnum_t))
end function gtpsa_add_zr

function gtpsa_add_rz (r, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(r, kind=c_cnum_t))
end function gtpsa_add_rz

function gtpsa_add_zc (z1, c) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_add_zc

function gtpsa_add_cz (c, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_add_cz

!------------------------------------------------------------------------------------------
! interface: operator(-)

function gtpsa_subtract_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_sub(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_subtract_tt

function gtpsa_subtract_tz (t1, z2) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class), intent(in) :: z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_tsub(t1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_subtract_tz

function gtpsa_subtract_zt (z1, t2) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(tpsa_class), intent(in) :: t2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_subt(z1%ptr, t2%ptr, z_out%ptr)
end function gtpsa_subtract_zt

function gtpsa_subtract_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_sub(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_subtract_zz

function gtpsa_subtract_tr (t1, r) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_copy (t1%ptr, t_out%ptr) 
  call mad_tpsa_set0(t_out%ptr, one_cn, real(-r, c_num_t))
end function gtpsa_subtract_tr

function gtpsa_subtract_rt (r, t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_copy (t1%ptr, t_out%ptr) 
  call mad_tpsa_set0(t_out%ptr, -one_cn, real(r, c_num_t))
end function gtpsa_subtract_rt

function gtpsa_subtract_tc (t1, c) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(-c, kind=c_cnum_t))
end function gtpsa_subtract_tc

function gtpsa_subtract_ct (c, t1) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, -one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_subtract_ct

function gtpsa_subtract_zr (z1, r) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(-r, kind=c_cnum_t))
end function gtpsa_subtract_zr

function gtpsa_subtract_rz (r, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, -one_ccn, cmplx(r, kind=c_cnum_t))
end function gtpsa_subtract_rz

function gtpsa_subtract_zc (z1, c) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, one_ccn, cmplx(-c, kind=c_cnum_t))
end function gtpsa_subtract_zc

function gtpsa_subtract_cz (c, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, -one_ccn, cmplx(c, kind=c_cnum_t))
end function gtpsa_subtract_cz

!------------------------------------------------------------------------------------------
! interface: operator(*)

function gtpsa_multiply_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_mul(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_multiply_tt

function gtpsa_multiply_tz (t1, z2) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class), intent(in) :: z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_mult(z2%ptr, t1%ptr, z_out%ptr)
end function gtpsa_multiply_tz

function gtpsa_multiply_zt (z1, t2) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(tpsa_class), intent(in) :: t2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_mult(z1%ptr, t2%ptr, z_out%ptr)
end function gtpsa_multiply_zt

function gtpsa_multiply_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_mul(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_multiply_zz

function gtpsa_multiply_tr (t1, r) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_scl(t1%ptr, real(r, c_num_t), t_out%ptr)
end function gtpsa_multiply_tr

function gtpsa_multiply_rt (r, t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_scl(t1%ptr, real(r, c_num_t), t_out%ptr)
end function gtpsa_multiply_rt

function gtpsa_multiply_tc (t1, c) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_scl(z_out%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_tc

function gtpsa_multiply_ct (c, t1) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_scl(z_out%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_ct

function gtpsa_multiply_zr (z1, r) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_scl(z1%ptr, cmplx(r, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_zr

function gtpsa_multiply_rz (r, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_scl(z1%ptr, cmplx(r, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_rz

function gtpsa_multiply_zc (z1, c) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_scl(z1%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_zc

function gtpsa_multiply_cz (c, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_scl(z1%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_multiply_cz

!------------------------------------------------------------------------------------------
! interface: operator(/)

function gtpsa_divide_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_div(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_divide_tt

function gtpsa_divide_tz (t1, z2) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class), intent(in) :: z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_tdiv(t1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_divide_tz

function gtpsa_divide_zt (z1, t2) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(tpsa_class), intent(in) :: t2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_divt(z1%ptr, t2%ptr, z_out%ptr)
end function gtpsa_divide_zt

function gtpsa_divide_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_div(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_divide_zz

function gtpsa_divide_tr (t1, r) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_copy (t1%ptr, t_out%ptr) 
  call mad_tpsa_set0(t_out%ptr, real(one_cn/r, c_num_t), zero_cn)
end function gtpsa_divide_tr

function gtpsa_divide_rt (r, t1) result (t_out)
  type(tpsa_class), intent(in) :: t1
  type(tpsa_class) :: t_out, t_temp
  real(dp), intent(in) :: r
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_inv(t1%ptr, real(r, c_num_t), t_out%ptr)
end function gtpsa_divide_rt

function gtpsa_divide_tc (t1, c) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_inv(z_out%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_divide_tc

function gtpsa_divide_ct (c, t1) result (z_out)
  type(tpsa_class), intent(in) :: t1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(t1%ptr, mad_tpsa_same)
  call mad_ctpsa_cplx (t1%ptr, c_null, z_out%ptr) 
  call mad_ctpsa_inv(z_out%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_divide_ct

function gtpsa_divide_zr (z1, r) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, cmplx(one_ccn/r, kind=c_cnum_t), zero_ccn)
end function gtpsa_divide_zr

function gtpsa_divide_rz (r, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  real(dp), intent(in) :: r
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_inv(z1%ptr, cmplx(r, kind=c_cnum_t), z_out%ptr)
end function gtpsa_divide_rz

function gtpsa_divide_zc (z1, c) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_copy (z1%ptr, z_out%ptr) 
  call mad_ctpsa_set0(z_out%ptr, cmplx(one_ccn/c, kind=c_cnum_t), zero_ccn)
end function gtpsa_divide_zc

function gtpsa_divide_cz (c, z1) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  type(ctpsa_class) :: z_out
  complex(dp), intent(in) :: c
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_inv(z1%ptr, cmplx(c, kind=c_cnum_t), z_out%ptr)
end function gtpsa_divide_cz

!------------------------------------------------------------------------------------------
! interface: operator(**)

function gtpsa_power_tt (t1, t2) result (t_out)
  type(tpsa_class), intent(in) :: t1, t2
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_pow(t1%ptr, t2%ptr, t_out%ptr)
end function gtpsa_power_tt

function gtpsa_power_ti (t1, n) result (t_out)
  type(tpsa_class), intent(in) :: t1
  integer, intent(in) :: n
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_powi(t1%ptr, int(n, c_int), t_out%ptr)
end function gtpsa_power_ti

function gtpsa_power_tr (t1, r) result (t_out)
  type(tpsa_class), intent(in) :: t1
  real(dp), intent(in) :: r
  type(tpsa_class) :: t_out
  t_out%ptr = mad_tpsa_new(t1%ptr, mad_tpsa_same)
  call mad_tpsa_pown(t1%ptr, real(r, c_num_t), t_out%ptr)
end function gtpsa_power_tr

function gtpsa_power_zz (z1, z2) result (z_out)
  type(ctpsa_class), intent(in) :: z1, z2
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_pow(z1%ptr, z2%ptr, z_out%ptr)
end function gtpsa_power_zz

function gtpsa_power_zi (z1, n) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  integer, intent(in) :: n
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_powi(z1%ptr, int(n, c_int), z_out%ptr)
end function gtpsa_power_zi

function gtpsa_power_zr (z1, r) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  real(dp), intent(in) :: r
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_pown(z1%ptr, cmplx(r, kind = c_cnum_t), z_out%ptr)
end function gtpsa_power_zr

function gtpsa_power_zc (z1, c) result (z_out)
  type(ctpsa_class), intent(in) :: z1
  complex(dp), intent(in) :: c
  type(ctpsa_class) :: z_out
  z_out%ptr = mad_ctpsa_new(z1%ptr, mad_tpsa_same)
  call mad_ctpsa_pown(z1%ptr, cmplx(c, kind = c_cnum_t), z_out%ptr)
end function gtpsa_power_zc

end module
