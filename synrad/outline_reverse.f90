subroutine outline_reverse (outline1, outline2)

  use synrad_struct
  use synrad_interface, except => outline_reverse

  implicit none

  type (outline_struct) outline1, outline2

  integer n

  outline2 = outline1
  outline2%s_center = -outline2%s_center

!

  n = outline2%n_in
  outline2%ix_in_slide = n - outline2%ix_in_slide + 2
  outline2%in(1:n)%x =  outline2%in(n:1:-1)%x
  outline2%in(1:n)%s = -outline2%in(n:1:-1)%s
  outline2%in(2:n)%name      = outline2%in(n:2:-1)%name
  outline2%in(2:n)%blueprint = outline2%in(n:2:-1)%blueprint
  outline2%in(2:n)%phantom   = outline2%in(n:2:-1)%phantom

  n = outline2%n_out
  outline2%ix_out_slide = n - outline2%ix_out_slide + 2
  outline2%out(1:n)%x =  outline2%out(n:1:-1)%x
  outline2%out(1:n)%s = -outline2%out(n:1:-1)%s
  outline2%out(2:n)%name      = outline2%out(n:2:-1)%name
  outline2%out(2:n)%blueprint = outline2%out(n:2:-1)%blueprint
  outline2%out(2:n)%phantom   = outline2%out(n:2:-1)%phantom

end subroutine
