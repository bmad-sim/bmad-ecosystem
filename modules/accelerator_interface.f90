module accelerator_interface

  use taylor_mod

  interface
    subroutine ele_to_fibre (bmad_ele, fibre_ele, param, &
                                                  integ_order, steps)
      use accelerator_struct
      implicit none
      type (ele_struct) bmad_ele
      type (fibre) fibre_ele
      type (param_struct) :: param
      integer, optional :: integ_order, steps
    end subroutine
  end interface

  interface
    subroutine energy_to_kinetic (energy, particle, &
                                       gamma, kinetic, beta, p0c, brho)
      use accelerator_struct
      implicit none
      real*8, intent(in) :: energy
      real*8, intent(out), optional :: gamma
      real*8, intent(out), optional :: kinetic
      real*8, intent(out), optional :: beta
      real*8, intent(out), optional :: p0c
      real*8, intent(out), optional :: brho
      integer, intent(in), optional :: particle
    end subroutine
  end interface

  interface
    function kind_name (this_kind)
      use accelerator_struct
      implicit none
      integer this_kind
      character*20 kind_name
    end function
  end interface
 
  interface
    function map_coef (y, i, j, k, l, style)
      use accelerator_struct
      implicit none
      type (real_8) y(:)
      integer i
      integer, optional :: j
      integer, optional :: k
      integer, optional :: l
      integer, optional :: style
      real*8 map_coef
    end function
  end interface

  interface
    subroutine map_track1 (y, ele, param)
      use accelerator_struct
      implicit none
      type (real_8) y(:)
      type (ele_struct) ele
      type (param_struct) param
    end subroutine
  end interface

  interface
    subroutine ring_to_layout (ring, ptc_layout)
      use accelerator_struct
      implicit none
      type (ring_struct), intent(in) :: ring
      type (layout), intent(inout) :: ptc_layout
    end subroutine
  end interface

  interface
    subroutine sort_universal_terms (ut_in, ut_sorted)
      use s_tracking
      implicit none
      type (universal_taylor), intent(in) :: ut_in
      type (universal_taylor) :: ut_sorted
    end subroutine
  end interface

  interface
    subroutine type_layout (lay)
      use accelerator_struct
      implicit none
      type (layout) lay
    end subroutine
  end interface
 
  interface
    subroutine type_map (y)
      use accelerator_struct
      implicit none
      type (real_8) :: y(:)
    end subroutine
  end interface

  interface
    subroutine type_map1 (y, type0, n_dim, style)
      use accelerator_struct
      implicit none
      type (real_8), intent(in) :: y(:)
      integer n_dim
      logical type0
      integer, optional :: style
    end subroutine
  end interface

  interface
    subroutine type_fibre (fib)
      use accelerator_struct
      implicit none
      type (fibre), intent(in) :: fib
    end subroutine
  end interface

  interface
    subroutine type_real_8_taylors (y, switch_z)
      use accelerator_struct
      implicit none
      type (real_8) y(:)
      logical, optional :: switch_z
    end subroutine
  end interface

end module
