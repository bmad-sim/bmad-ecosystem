module pisa_mod

use bmad, only: rp

implicit none

type pool_struct
  integer name
  real(rp), allocatable :: x(:)  !variables
end type

type pop_struct
  integer name
  real(rp), allocatable :: x(:)  ! variables
  real(rp), allocatable :: x_phys(:)  !variables
  real(rp), allocatable :: o(:)  !objectives
  real(rp), allocatable :: c(:)  !constraints
end type

type breeder_params_struct
  real(rp) :: cross_p = 0.8
  real(rp) :: eta = 1.0
  real(rp) :: mutate_p = 0.1
  real(rp) :: mutate_delta(100)
end type

contains

subroutine find_empty_pop_slot(pop,ix)
  implicit none

  type(pop_struct) :: pop(:)
  integer ix

  do ix=1,size(pop)
    if(pop(ix)%name .eq. -1) exit
  enddo
end subroutine

subroutine delete_the_dead(pop,arc,narc)
  implicit none

  type(pop_struct) :: pop(:)
  integer arc(:)
  integer narc

  integer i,j
  logical dead

  do i=1,size(pop)
    dead = .true.
    do j=1,narc
      if( pop(i)%name == arc(j) ) then
        dead = .false.
        exit
      endif
    enddo
    if( dead ) then
      pop(i)%name = -1
    endif
  enddo
end subroutine

subroutine pisa_cfg_parser(prefix, alpha, mu, lambda, dim, con)
  implicit none

  character(*) prefix
  integer alpha, mu, lambda, dim, con
  integer iostat
  character(30) cfg_str
  integer cfg_int

  alpha=-1
  mu=-1
  lambda=-1
  con=-1

  open(21,file=trim(prefix)//'cfg',action='read',status='old',iostat=iostat)
  if(iostat .ne. 0) then
    write(*,'(a,a,a)') "PISA configuration file ", trim(prefix)//'cfg', " not found.  Terminating."
    error stop
  endif
  do while(.true.)
    read(21,*,iostat=iostat) cfg_str, cfg_int
    if(iostat .ne. 0) exit

    if( (trim(cfg_str) .eq. 'alpha') .or. (trim(cfg_str) .eq. 'initial_population_size') ) then
      alpha = cfg_int
    elseif( (trim(cfg_str) .eq. 'mu') .or. (trim(cfg_str) .eq. 'parent_set_size') ) then
      mu = cfg_int
    elseif( (trim(cfg_str) .eq. 'lambda') .or. (trim(cfg_str) .eq. 'offspring_set_size') ) then
      lambda = cfg_int
    elseif( (trim(cfg_str) .eq. 'dim') .or. (trim(cfg_str) .eq. 'objectives') ) then
      dim = cfg_int
    elseif( (trim(cfg_str) .eq. 'constraints') ) then
      con = cfg_int
    else
      write(*,*) "Unidentifiable content in cfg file.  Terminating program."
      error stop
    endif
  enddo
  close(21)

  if(alpha == -1) then
    write(*,*) "alpha not set properly in "//trim(prefix)//"cfg"
    error stop
  endif
  if(mu == -1) then
    write(*,*) "mu not set properly in "//trim(prefix)//"cfg"
    error stop
  endif
  if(lambda == -1) then
    write(*,*) "lambda not set properly in "//trim(prefix)//"cfg"
    error stop
  endif
  if(dim == -1) then
    write(*,*) "dim not set properly in "//trim(prefix)//"cfg"
    error stop
  endif
  if(con == -1) then
    write(*,*) "con not set properly in "//trim(prefix)//"cfg"
    error stop
  endif
end subroutine

subroutine increment_ptr(ptr,ub)
  implicit none

  integer ptr, ub

  ptr = ptr + 1
  if(ptr .gt. ub) ptr = 1
end subroutine

function name_to_ix(pop,name) result(ix)
  implicit none

  type(pop_struct) pop(:)
  integer name, ix

  do ix=1,size(pop)
    if( name == pop(ix)%name ) exit
  enddo
end function

subroutine read_pisa_indexes(prefix,what,n,items)
  implicit none

  character(20) prefix
  character(3) what
  integer n
  integer items(:)
  integer i
  character(10) end_str

  open(21,file=trim(prefix)//what,action='read')
  read(21,'(i8)') n
  do i=1,n
    read(21,'(i8)') items(i)
  enddo
  read(21,*) end_str
  close(21)
  open(21,file=trim(prefix)//what,action='write')
  write(21,'(i1)') 0
  close(21)
end subroutine

subroutine write_pop_pisa(pop,filename)
  implicit none

  type(pop_struct) pop(:)
  integer i, N
  character(*) filename

  N = size(pop)

  open(21,file=filename,action='write')
  write(21,'(i8)') N*(size(pop(1)%o(:))+size(pop(1)%c(:))+1)
  do i=1,N
    write(21,'(i8,20es18.8)') pop(i)%name, pop(i)%o(:), pop(i)%c(:)
  enddo
  write(21,'(a)') 'END'
  close(21)
end subroutine write_pop_pisa

subroutine write_objective_report(pop,gen_num,filename)
  implicit none

  type(pop_struct) pop(:)
  integer gen_num
  character(*) filename

  integer i
  real feasible

  open(22,file=filename,access='append')
  write(22,'(a,i8)') "# Generation ", gen_num
  do i=1,size(pop)
    if(pop(i)%name .gt. 0) then
      feasible = merge(1.0,0.0,all( pop(i)%c(:) .ge. 0.0d0 ))
      write(22,'(i8,20es18.8)') pop(i)%name, pop(i)%o(:), feasible
    endif
  enddo
  write(22,*)
  write(22,*)
  close(22)
end subroutine write_objective_report

subroutine write_constraint_report(pop,gen_num,filename)
  implicit none

  type(pop_struct) pop(:)
  integer gen_num
  character(*) filename

  integer i

  open(22,file=filename,access='append')
  write(22,'(a,i8)') "# Generation ", gen_num
  do i=1,size(pop)
    if(pop(i)%name .gt. 0) then
      write(22,'(i8,20es18.8)') pop(i)%name, pop(i)%c(:)
    endif
  enddo
  write(22,*)
  write(22,*)
  close(22)
end subroutine write_constraint_report

subroutine block_on_pisa_status(polli,prefix)
  implicit none

  integer polli
  character(*) prefix

  integer sta

  sta = 1
  do while( sta .ne. 2 )
    call milli_sleep(polli)
    sta = poll_state(prefix)
  enddo
end subroutine

function poll_state(prefix)
  implicit none

  character(*) prefix
  integer iostat
  integer poll_state

  poll_state = -1
  open(21,file=trim(prefix)//'sta',action='read',iostat=iostat)
  if(iostat .eq. 0) then
    read(21,'(i1)',iostat=iostat) poll_state
  endif
  close(21)
end function

subroutine write_state(prefix,sta)
  implicit none

  character(*) prefix
  integer sta

  open(21,file=trim(prefix)//'sta',action='write')
  write(21,'(i1)') sta
  close(21)
end subroutine

subroutine kangal_breeder(pop, sel, pool, pool_ptr_b, breeder_params)
  ! http://www.iitk.ac.in/kangal/resources.shtml
  implicit none

  type(pop_struct) pop(:)
  integer sel(:)
  type(pool_struct) pool(:)
  integer pool_ptr_b
  type(breeder_params_struct) breeder_params

  integer mu
  real(rp) expo
  real(rp) u, Bq
  integer i, j
  integer child_a_name, child_b_name
  integer parent_a_ix, parent_b_ix, child_a_ix, child_b_ix
  integer ix, Nvars

  Nvars = size(pop(1)%x(:))
  mu = size(sel)

  do j=1,mu,2
    parent_a_ix = name_to_ix(pop, sel(j))
    parent_b_ix = name_to_ix(pop, sel(j+1))

    child_a_name = pool(pool_ptr_b)%name + 1
    call increment_ptr(pool_ptr_b,size(pool))
    child_a_ix = pool_ptr_b
    pool(child_a_ix)%name = child_a_name

    child_b_name = pool(pool_ptr_b)%name + 1
    call increment_ptr(pool_ptr_b,size(pool))
    child_b_ix = pool_ptr_b
    pool(child_b_ix)%name = child_b_name

    expo = 1.0d0/(breeder_params%eta+1.0d0)

    !cross parents
    do i=1,Nvars
      call random_number(u)
      if( u .le. breeder_params%cross_p ) then
        ! gene crossover
        call random_number(u)
        if (u .le. 0.5) then
          Bq = (2.0d0*u) ** expo
        else
          Bq = (0.5d0/(1.0d0-u)) ** expo
        endif
        pool(child_a_ix)%x(i)=0.5*( (1+Bq)*pop(parent_a_ix)%x(i) + (1-Bq)*pop(parent_b_ix)%x(i) )
        pool(child_b_ix)%x(i)=0.5*( (1-Bq)*pop(parent_a_ix)%x(i) + (1+Bq)*pop(parent_b_ix)%x(i) )
      else
        ! copy the gene
        pool(child_a_ix)%x(i)=pop(parent_a_ix)%x(i)
        pool(child_b_ix)%x(i)=pop(parent_b_ix)%x(i)
      endif
    enddo

    !mutate genes
    ! "Improving NSGA-II with an Adaptive Mutation Operator", Carvalho, Araujo.
    do i=1,Nvars
      call random_number(u)
      if( u .le. breeder_params%mutate_p ) then
        call random_number(u)
        if( u .lt. 0.5 ) then
          Bq = (2.0d0*u) ** expo - 1.0d0
        else
          Bq = 1.0d0 - (2.0d0*(1.0d0-u)) ** expo
        endif
        pool(child_a_ix)%x(i) = pool(child_a_ix)%x(i) + Bq*breeder_params%mutate_delta(i)
      endif
      call random_number(u)
      if( u .le. breeder_params%mutate_p ) then
        call random_number(u)
        if( u .lt. 0.5 ) then
          Bq = (2.0d0*u) ** expo - 1.0d0
        else
          Bq = 1.0d0 - (2.0d0*(1.0d0-u)) ** expo
        endif
        pool(child_b_ix)%x(i) = pool(child_b_ix)%x(i) + Bq*breeder_params%mutate_delta(i)
      endif
    enddo
  enddo
end subroutine

subroutine niave_breeder(pop,parent_a,parent_b,child,mutation_rate)
  implicit none

  type(pop_struct) pop(:)
  integer parent_a, parent_b, child
  real mutation_rate

  real r
  integer i
  integer child_ix, ix

  !name child
  do child_ix=1,size(pop)
    if( pop(child_ix)%name == -1 ) then
      pop(child_ix)%name = child
      exit
    endif
  enddo

  !cross parents
  do i=1,size(pop(1)%x(:))
    call random_number(r)
    if( r .lt. 0.5 ) then
      pop(child_ix)%x(i)=pop(name_to_ix(pop,parent_a))%x(i)
    else
      pop(child_ix)%x(i)=pop(name_to_ix(pop,parent_b))%x(i)
    endif
  enddo

  !add mutation
  call random_number(r)
  ix = ceiling(r*size(pop(1)%x(:)))  !pick a random gene to mutate
  call random_number(r)
  pop(child_ix)%x(ix) = pop(child_ix)%x(ix) + mutation_rate*(2.0*r-1.0)  !mutate that gene
end subroutine

end module
