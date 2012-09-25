subroutine outline_concat (outline1, outline2, outline3)

  use synrad_struct
  use synrad_interface, except => outline_concat

  implicit none

  type (outline_struct) outline22, outline1, outline2, outline3

  integer n1, n2, n3

!

  if (outline1%zero_is_center .and. outline2%zero_is_center) then
    print *, 'ERROR: CONCATINATING TWO OUTLINES BOTH OF WHICH HAVE ZERO_IS_CENTER'
    print *, '      ', outline1%name, outline2%name
    if (global_com%exit_on_error) call err_exit
  endif

  outline22 = outline2
  outline3 = outline1
  outline3%zero_is_center = outline1%zero_is_center .or. &
                                            outline22%zero_is_center

! negative_x_wall

  n1 = outline1%n_in
  n2 = outline22%n_in
  n3 = outline3%n_in + n2
  outline3%n_in = n3

  outline3%in(n1+1:n3) = outline22%in(1:n2)

  if (outline22%zero_is_center) then
    outline3%in(1:n1)%s = outline3%in(1:n1)%s - &
                                  (outline3%in(n1)%s - outline3%in(n1+1)%s)
  else
    outline3%in(n1+1:n3)%s = outline3%in(n1+1:n3)%s + &
                                  (outline3%in(n1)%s - outline3%in(n1+1)%s)
  endif

  if (outline3%in(n1)%x == outline3%in(n1+1)%x) then
    outline3%in(n1+1:n3-1) = outline3%in(n1+2:n3)
    outline3%n_in = n3 - 1
  endif

! positive_x_wall

  n1 = outline1%n_out
  n2 = outline22%n_out
  n3 = outline3%n_out + n2
  outline3%n_out = n3

  outline3%out(n1+1:n3) = outline22%out(1:n2)

  if (outline22%zero_is_center) then
    outline3%out(1:n1)%s = outline3%out(1:n1)%s - &
                                  (outline3%out(n1)%s - outline3%out(n1+1)%s)
  else
    outline3%out(n1+1:n3)%s = outline3%out(n1+1:n3)%s + &
                                  (outline3%out(n1)%s - outline3%out(n1+1)%s)
  endif

  if (outline3%out(n1)%x == outline3%out(n1+1)%x) then
    outline3%out(n1+1:n3-1) = outline3%out(n1+2:n3)
    outline3%n_out = n3 - 1
  endif

end subroutine
