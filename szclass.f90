module szclass
  use precision1, only : dp

  implicit none

  type spline_data
      integer :: n
      integer :: status
      real(dp), dimension(:), allocatable :: x
      real(dp), dimension(:), allocatable :: f
      real(dp), dimension(:), allocatable :: ddf
  end type spline_data

contains

    subroutine init_spline_data(npts,spl)
        integer, intent(in) :: npts
        type(spline_data), intent(out) :: spl
        allocate(spl%x(npts))
        allocate(spl%f(npts))
        allocate(spl%ddf(npts),STAT=spl%status)
    end subroutine init_spline_data

    subroutine del_spline_data(spl)
        type(spline_data), intent(inout) :: spl
        deallocate(spl%ddf)
        deallocate(spl%f)
        deallocate(spl%x)
    end subroutine del_spline_data

end module szclass
