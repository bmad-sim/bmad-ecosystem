!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_key_info_to_str (ix_key, ix_min_key, ix_max_key, key_str, header_str)

use tao_interface, dummy => tao_key_info_to_str

implicit none

type (tao_var_struct), pointer :: var

integer ix_key, ix_min_key, ix_max_key
integer i, n, m, p, ix_var
integer :: j_var1, j_att

real(rp) :: y_here, norm, v, x1, x2, y1, y2

character(*) key_str
character(*) header_str
character(60) fmt, fmt2
character(15) model_str, val0_str, delta_str
character(4) exp_str
character(24) :: r_name = 'tao_key_info_to_str'

! Compute widths of var1 and attrib fields.

j_var1 = 4
j_att = 5

do i = ix_min_key, ix_max_key
  if (i > ubound(s%key, 1)) cycle
  ix_var = s%key(i)
  if (ix_var < 1) cycle
  j_var1 = max(j_var1, len_trim(tao_var1_name(s%var(ix_var))))
  j_att  = max(j_att,  len_trim(tao_var_attrib_name(s%var(ix_var))))
enddo

! Write header 

write (fmt, '(a, i5, a, i2, a)') '(a, ', j_var1-2, 'x, a, ', j_att, 'x, a)'
write (header_str, fmt) 'Name', 'Attrib', 'Value         Value0          Delta'

! Write key info

key_str = ''
if (ix_key > ubound(s%key, 1)) return
  
ix_var = s%key(ix_key)
if (ix_var < 1) return

var => s%var(ix_var)
v = max(abs(var%model_value), abs(var%key_val0), abs(var%key_delta))
if (v == 0) then
  n = 0
  m = 2
else
  m = 1.001 * log10(v)
  n = 3 * floor(m/3.0)
  p = 3 - (m - n)
endif

exp_str = ''

if (m >= -1 .and. m <= 1) then
  fmt2 = '(f11.4, a)'
  n = 0
elseif (m == 2) then
  fmt2 = '(f11.2, a)'
  n = 0
else
  write (fmt2, '(a, i1, a)') '(f7.', p, ', a)'
  write (exp_str, '(a, i3.2)') 'E', n
  if (exp_str(2:2) == ' ') exp_str(2:2) = '+'
endif

model_str = ''; val0_str = ''; delta_str = ''

write (model_str, fmt2) var%model_value / 10.0**n, exp_str
write (val0_str,  fmt2) var%key_val0 / 10.0**n, exp_str
write (delta_str, fmt2) var%key_delta / 10.0**n, exp_str

write (fmt, '(3(a, i2.2))') '(a', j_var1, ', 2x, a', j_att, ', 3a15, 3x, l1)'
write (key_str, fmt) tao_var1_name(var), tao_var_attrib_name(var), model_str, &
                        val0_str, delta_str, var%useit_opt

end subroutine tao_key_info_to_str
