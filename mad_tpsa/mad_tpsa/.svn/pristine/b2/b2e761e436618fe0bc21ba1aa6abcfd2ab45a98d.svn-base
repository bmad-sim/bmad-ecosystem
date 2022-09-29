program tpsa_test

use gtpsa_mod
implicit none

type (tpsa_class) t1, t2, t3
type (ctpsa_class) z2
type (c_ptr) :: d
real(c_num_t)  :: vec(100)
integer(c_idx_t) :: ii(100)
integer n, ovec(10)
!

d = mad_desc_newn(4, 3_1)
t1 = gtpsa_newd(d, mad_tpsa_default)

n = 21
vec(1:n) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0, 20]
call t1%mono_setv(0, n, vec(1:n))
call t1%print("ini t1", 1d-30, 0, c_null);

print *, t1%n_mono()
stop




print *, t1%mono_get_by_expn([1,1,0,0])
print *, t1%mono_get_by_expn()
stop

ii(1:4) = [1,2,4,3]
call mad_tpsa_convert (t1%ptr, t2%ptr, 4, ii(1:4), 0)
call t2%print("t2", 1d-30, 0, c_null);
stop

t2 = gtpsa_newd(d, mad_tpsa_default)

z2 = gtpsa_new(t1, mad_tpsa_same)

print *, 'Num var:', t1%n_var()
print *, 'Num knobs:', t1%n_knobs()
print *, 'Num main var:', t1%n_main_var()
print *, 'Order max:', t1%order_max()
print *, 'Order knobs:', t1%order_knobs()
call t1%order_vec(ovec)
print '(a, 10i4)', 'Order vec:', ovec
stop


!vec = [2.0d0, 0.0d0, 0.0d0]
!call gtpsa_setv (t2, 0, 1+2, vec);

call t1%print("ini t1", 0d0, 0, c_null);
!call gtpsa_print(t2, "ini t1", 0d0, c_null);

print *, 'Here: ', t1%n_var()
print *, 'Here: ', t1%coef([0,1])

t3 = sin(t1)**2 + cos(t1)**2
call t3%print("z2", 0d0, 0, c_null);

print *, 'Here: ', t3%coef([1,2])

t3 = cosh(t1)**2 - sinh(t1)**2
call t3%print("z2", 0d0, 0, c_null);



end program
