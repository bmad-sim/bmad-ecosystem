program coarray_example

  use bmad
  use ptc_layout_mod

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: co(:)
  type (coord_struct), allocatable :: orb(:)
  type (coord_struct) vecs(100)
  type(taylor_struct) A(1:6), A_inverse(1:6), dhdj(1:6)
  type(taylor_struct) A1(1:6), A1_inverse(1:6)
  type(taylor_struct) one_turn_map(6)

  ! declare variables as coarrays using [:]
  type(coord_struct), allocatable :: lab_vecs(:,:)[:]
  type(coord_struct), allocatable ::   A_vecs(:,:)[:]
  type(coord_struct), allocatable ::  A1_vecs(:,:)[:]

  type(coord_struct), allocatable :: lab_vecs_buf(:)
  type(coord_struct), allocatable ::   A_vecs_buf(:)
  type(coord_struct), allocatable ::  A1_vecs_buf(:)

  integer i, j
  integer vec_num
  integer n_vecs
  integer n_turns, dims
  integer tracking_method

  real(rp) x_init, y_init
  real(rp) A_coords(6)
  real(rp) A1_coords(6)
  real(rp) pz

  ! coarray stuff: used to build a schedule of which worker processes which vector
  integer my_worker_num
  integer num_workers
  integer counter
  integer, allocatable :: work_schedule(:)

  character*60 in_file
  character*60 vec_file
  character*100 lat_file
    
  namelist / parameters / lat_file, &
                          tracking_method, &
                          vec_file, &
                          n_turns, &          !Number of turns particle must survive
                          dims               !either 4 or 6

  ! coarray stuff
  num_workers = num_images()    !query number of process images that were spawned
  my_worker_num = this_image()  !query what is my unique worker number

  call getarg(1,in_file)
  open(10, file = in_file)
  read(10, nml = parameters)
  close(10)

  open(11, file = vec_file)
  read(11, '(I3)') n_vecs
  IF(my_worker_num == 1) WRITE(*,*) "Reading ", n_vecs, " initial vectors."
  do i=1,n_vecs
    read(11, *) vecs(i)%vec
  enddo
  close(11)

  !coarray stuff
  !Have worker 1 parse lattice first, to avoid multiple processes writing digested file at once.
  IF(my_worker_num .eq. 1) THEN
    call bmad_parser(lat_file, lat)
  ENDIF
  sync all  !make all processes wait for worker 1 to finish parsing lattice.
  IF(my_worker_num .ne. 1) THEN
    call bmad_parser(lat_file, lat)
  ENDIF


  IF( dims == 4 ) THEN
    CALL set_on_off(rfcavity$, lat, off$)
  ELSEIF( dims == 6 ) THEN
    CALL set_on_off(rfcavity$, lat, on$)
  ELSE
    WRITE(*,*) "Error setting dims.  Terminate."
    STOP
  ENDIF
  CALL closed_orbit_calc(lat,co,dims)
  CALL lat_make_mat6(lat, -1, co)
  CALL twiss_at_start(lat)
  CALL twiss_propagate_all(lat)

  CALL lat_to_ptc_layout(lat)
  CALL ptc_layouts_resplit( 0.04d0, 0.0d0, .false., 0.0d0, 0.0d0)

  pz = 0.0
  CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, pz = pz)
  CALL normal_form_taylors(one_turn_map, .false., dhdj, A, A_inverse)
  CALL truncate_taylor_to_order(A_inverse,1,A1_inverse)

  if( tracking_method .gt. 0 ) THEN
    do i = 1, lat%n_ele_max
      lat%ele(i)%tracking_method = tracking_method
    enddo
  endif

  ALLOCATE(orb(0:lat%n_ele_track))

  ! coarray stuff:  notice the [*] required for allocatable coarrays.
  ALLOCATE(lab_vecs(n_vecs, n_turns)[*])
  ALLOCATE(  A_vecs(n_vecs, n_turns)[*])
  ALLOCATE( A1_vecs(n_vecs, n_turns)[*])

  ALLOCATE(lab_vecs_buf(n_turns))
  ALLOCATE(  A_vecs_buf(n_turns))
  ALLOCATE( A1_vecs_buf(n_turns))

  ! coarray stuff:  build the work schedule
  IF(my_worker_num == 1) write(*,*) "Number of Images: ", num_workers
  allocate(work_schedule(num_workers))
  work_schedule = 0
  counter = 1
  do i=1,n_vecs
    work_schedule(counter) = work_schedule(counter) + 1 
    counter = counter + 1
    IF(counter .gt. num_workers) counter=1
  enddo

  do i=1,work_schedule(my_worker_num)
    vec_num = my_worker_num + (i-1)*num_workers
    orb(0) = vecs(vec_num)
    do j=1,n_turns
      CALL track_all(lat, orb)
      A_coords = track_taylor(orb(lat%n_ele_track)%vec, A_inverse)
      A1_coords = track_taylor(orb(lat%n_ele_track)%vec, A1_inverse)

      lab_vecs_buf(j)%vec = orb(lat%n_ele_track)%vec
        A_vecs_buf(j)%vec = A_coords
       A1_vecs_buf(j)%vec = A1_coords

      orb(0) = orb(lat%n_ele_track)
    enddo

    ! coarray stuff:  copy the trajectory over to worker number 1
    lab_vecs(vec_num,:)[1] = lab_vecs_buf(:)
    A_vecs(vec_num,:)[1] = A_vecs_buf(:)
    A1_vecs(vec_num,:)[1] = A1_vecs_buf(:)
  enddo

  ! coarray stuff:  the syntax for sync directives is a bit unconventional.
  ! 'sync all' requires that all processes reach this point before any continue on.
  sync all

  IF( my_worker_num .eq. 1 ) THEN  !worker 1 write the output files
    open(20,file='track.out')
    open(21,file='track_A.out')
    open(22,file='track_A1.out')

    write(20,'(A8,6A14)') "# turn", "x", "px", "y", "py", "z", "pz"
    write(21,'(A8,6A14)') "# turn", "x", "px", "y", "py", "z", "pz"
    write(22,'(A8,6A14)') "# turn", "x", "px", "y", "py", "z", "pz"

    do i=1,n_vecs
      write(20,'(A,6ES14.4)') "# vec = ", vecs(i)%vec
      write(21,'(A,6ES14.4)') "# vec = ", vecs(i)%vec
      write(22,'(A,6ES14.4)') "# vec = ", vecs(i)%vec
      do j=1,n_turns
        write(20,'(I8,6ES14.4)') j, lab_vecs(i,j)%vec
        write(21,'(I8,6ES14.4)') j,   A_vecs(i,j)%vec
        write(22,'(I8,6ES14.4)') j,  A1_vecs(i,j)%vec
      enddo
      write(20,*)
      write(20,*)
      write(21,*)
      write(21,*)
      write(22,*)
      write(22,*)
    enddo

    close(20)
    close(21)
    close(22)
  endif

end program
                              
