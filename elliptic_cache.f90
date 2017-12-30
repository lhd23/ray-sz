module elliptic_cache
  use precision1, only : dp
  implicit none

  type ellip_obj
      real(dp) :: upper_lim
      integer, dimension(3) :: indices
      complex(dp) :: val
  end type ellip_obj
  
  type(ellip_obj), dimension(:), allocatable, save :: cache
!$omp threadprivate(cache)

  integer, save :: calls=0,count=0

  interface assignment (=) !overload assignment
      module procedure equal_ellip_obj
      module procedure equal_ellip_obj_arr
  end interface
  
contains

    subroutine search_cache(ulim,p_indices,val,exists)
! Search cache and return integral value if found
        implicit none
        real(dp), intent(in) :: ulim
        integer, intent(in) :: p_indices(3)
        complex(dp), intent(out) :: val
        logical, intent(out) :: exists
        integer iell
        val=0._dp !default value
        exists=.false.
        calls=calls+1
        if (allocated(cache)) then
            do iell=1,size(cache)
                if (all(cache(iell)%indices == p_indices) .and. &
                    abs(ulim-cache(iell)%upper_lim) < 1e-15) then
                    exists=.true.
                    val=cache(iell)%val
                    count=count+1
                    return
                end if
            end do
        end if
    end subroutine search_cache

    subroutine update_cache(ulim,p_indices,val)
        implicit none
        real(dp), intent(in) :: ulim
        integer, intent(in) :: p_indices(3)
        complex(dp), intent(in) :: val
        type(ellip_obj), allocatable :: temp(:)
        integer ni
        if (.not. allocated(cache)) then
            allocate(cache(1))
            ni=0
        else !extend the size of cache
            ni=size(cache)
            allocate(temp(ni))
            temp=cache
            deallocate(cache)
            allocate(cache(ni+1))
            cache(1:ni)=temp
            deallocate(temp)
        end if
        cache(ni+1)%upper_lim=ulim
        cache(ni+1)%indices=p_indices
        cache(ni+1)%val=val
    end subroutine update_cache

    subroutine equal_ellip_obj(new,old)
        implicit none
        type(ellip_obj), intent(out) :: new
        type(ellip_obj), intent(in) :: old
        new%upper_lim=old%upper_lim
        new%indices=old%indices
        new%val=old%val
    end subroutine equal_ellip_obj

    subroutine equal_ellip_obj_arr(new,old)
        implicit none
        type(ellip_obj), intent(out), dimension(:) :: new
        type(ellip_obj), intent(in), dimension(:) :: old
        integer i
        do i=1,size(old)
            new(i)%upper_lim=old(i)%upper_lim
            new(i)%indices=old(i)%indices
            new(i)%val=old(i)%val
        end do
    end subroutine equal_ellip_obj_arr
    
end module elliptic_cache
