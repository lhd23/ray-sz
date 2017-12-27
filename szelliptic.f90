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

module szelliptic
  use precision1, only : dp
  use constants
  use utils, only : fcubic_roots,zbrent, &
                    rj,rd,rf
  implicit none

  type szcmplx
      real(dp) :: p !cubic coefficient
      real(dp) :: q
      complex(dp), dimension(3) :: roots
      complex(dp), dimension(3) :: roots_inv
      complex(dp), dimension(3) :: f
      real(dp) :: prd_roots !x1.x2.x3
      real(dp) :: xi !sqrt(-x1.x2.x3)      
      real(dp) :: p3 !1/sqrt(-x1.x2.x3)
  end type szcmplx

  interface szellip_1111
      module procedure szellip_1111
      module procedure szellip_1111_roots
  end interface szellip_1111

  private ! szellip_113_,szellip_115_
  public :: szcmplx,init_szcmplx,szellip_1111,rf44, &
          szellip_1311,szellip_1131,szellip_1113, &
          szellip_1511,szellip_1151,szellip_1115, &
          szellip_3333,szellip_1333, &
          szellip_5555,szellip_3555,szellip_1555

contains

    subroutine init_szcmplx(p,q,szc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes common complex quantities needed to evaluate   !
! certain formulas involving Carlson R-functions.         !
!                                                         !
! Checks for unphysical case: Rdot^2 is non-negative the  !
! quantities k and M in the equation                      !
!                                                         !
!   Rdot^2 = 2M/R - k + Lambda * R^2/3                    !
!                                                         !
! must be such that it is also non-negative.              !
! Alternatively, the square root of the RHS must be real. !
! This is certain when the discriminant -4p^3-27q^2 of    !
! the depressed cubic                                     !
!                                                         !
!   R^3 + p * R + q                                       !
!                                                         !
!     p = -k / (Omegal*H0^2)                              !
!     q = 2M / (Omegal*H0^2)                              !
!                                                         !
! is negative and the real root is negative. The real     !
! root R1 is given by                                     !
!                                                         !
!   R1 = S + T                                            !
!                                                         !
!     S = (-q/2 + sqrt(D))^1/3                            !
!     T = (-q/2 - sqrt(D))^1/3                            !
!                                                         !
! where D = -Discriminant/108. Note R1 > 0 when           !
!                                                         !
!   k^3 < 27 * Omegal * H0^2 * M^2                        !
!                                                         !
! Given that M >= 0 and for a void almost always          !
! k < 0 this mostly hold                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: p,q
        type(szcmplx), intent(out) :: szc
        complex(dp), dimension(3) :: r !roots
        real(dp) :: prd
! Clear cache for every new r (or p and q)
        if (allocated(cache)) deallocate(cache) 
        call fcubic_roots(p,q,r)
! Should be one real root only and less than 0
        if (real(r(1)) > 0._dp) then
            write(*,'(A,X,A,X,A,X,A)') &
              'warning: real root > 0.', &
              'Root should be negative to ensure', &
              'square roots in integrands', &
              'are real valued'
        end if
        szc%p=p
        szc%q=q
        szc%roots=r
        szc%roots_inv=1._dp/szc%roots
        szc%f(1)=1._dp/((r(1)-r(2))*(r(1)-r(3)))
        szc%f(2)=1._dp/((r(2)-r(1))*(r(2)-r(3)))
        szc%f(3)=1._dp/((r(3)-r(1))*(r(3)-r(2)))
        prd=real(product(szc%roots))
        if (prd < 0._dp) then
            szc%prd_roots=prd
            szc%xi=sqrt(abs(prd))
            szc%p3=1._dp/szc%xi
        else
            STOP &
             'init_szcmplx: sqrt(-x1.x2.x3) is imag'
        end if
    end subroutine init_szcmplx

    real(dp) function szellip_1111_roots(ulim,roots)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the integral from 0 to ulim of                          !
!                                                                 !
!    x^1/2 . (x-x1)^{-1/2} . (x-x2)^{-1/2} . (x-x3)^{-1/2} dx     !
!                                                                 !
! normalised by 1 / sqrt(Omegal * H0^2)                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: ulim
        complex(dp), dimension(3), &
                intent(in) :: roots
        complex(dp), dimension(3) :: x
        complex(dp) :: ulim_inv
        real(dp) :: w,f        
! Carlson arguments are of form:  1/ulim - 1/root
        ulim_inv=dcmplx(1._dp/ulim)
        x=ulim_inv-1._dp/roots
        w=real(product(roots))
        f=TWOTHIRD/sqrt(abs(w))
        szellip_1111_roots=f*real(rj(x(1),x(2),x(3),ulim_inv))
    end function szellip_1111_roots

    real(dp) function szellip_1111(ulim,szc)
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        complex(dp), dimension(3) :: x
        complex(dp) ulim_inv
        real(dp) f
        ulim_inv=dcmplx(1._dp/ulim)
        x=ulim_inv-szc%roots_inv
        f=TWOTHIRD*szc%p3
        szellip_1111=f*real(rj(x(1),x(2),x(3),ulim_inv))
    end function szellip_1111

    real(dp) function rf44(a2,b2,a4,ulim)
        implicit none
        real(dp), intent(in) :: a2,b2,a4,ulim
        real(dp), parameter :: two=2._dp
        real(dp) :: c2,ulim_inv,xi,a4_inv,M2,Lp2,Lm2
        c2=a2**2+b2**2
        ulim_inv=1._dp/ulim
        a4_inv=1._dp/a4
        xi=a4*c2
! Carlson 1991 (4.4)
        M2=two*xi*((a2/c2-ulim_inv) &
           +sqrt((a2/c2-ulim_inv)**2+(b2/c2)**2))
        Lp2=M2+two*xi*(a4_inv-a2/c2 &
           +sqrt((a2/c2-a4_inv)**2+(b2/c2)**2))
        Lm2=M2+two*xi*(a4_inv-a2/c2 &
           -sqrt((a2/c2-a4_inv)**2+(b2/c2)**2))
        rf44=two*rf(M2,Lm2,Lp2)
    end function rf44

    complex(dp) function szellip_1311(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-3,-1,-1/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_113_(szc,ulim,(/2,3,1/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1311

    complex(dp) function szellip_1131(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-1,-3,-1/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_113_(szc,ulim,(/3,1,2/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1131

    complex(dp) function szellip_1113(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-1,-1,-3/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_113_(szc,ulim,(/1,2,3/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1113

    complex(dp) function szellip_113_(szc,ulim,perm)
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, intent(in), dimension(3) :: perm
        complex(dp) :: ulim_inv,x2,x3,x4,xi,a2,a3,a4
        x2=szc%roots(perm(1))
        x3=szc%roots(perm(2))
        x4=szc%roots(perm(3))
        ulim_inv=dcmplx(1._dp/ulim)
        xi=szc%xi
        a2=ulim_inv-1._dp/x2
        a3=ulim_inv-1._dp/x3
        a4=ulim_inv-1._dp/x4
        szellip_113_=-TWOTHIRD/(x4*xi)*rd(a2,a3,a4)
    end function szellip_113_

    complex(dp) function szellip_1511(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-5,-1,-1/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_115_(szc,ulim,(/2,3,1/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1511

    complex(dp) function szellip_1151(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-1,-5,-1/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_115_(szc,ulim,(/3,1,2/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1151

    complex(dp) function szellip_1115(szc,ulim) result(y)
        use elliptic_cache
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        integer, dimension(3) :: p_ind
        complex(dp) :: val
        logical exists
        p_ind=(/-1,-1,-5/)
        call search_cache(ulim,p_ind,val,exists)
        if (exists) then
            y=val
            return
        else
            y=szellip_115_(szc,ulim,(/1,2,3/))
            call update_cache(ulim,p_ind,y)
            return
        end if
    end function szellip_1115

    complex(dp) function szellip_115_(szc,ulim,perm)
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim        
        integer, intent(in), dimension(3) :: perm
        complex(dp) :: ulim_inv,x2,x3,x4,a2,a3,a4
        complex(dp) :: s,t1,t2,t3
        real(dp) :: xi
        x2=szc%roots(perm(1))
        x3=szc%roots(perm(2))
        x4=szc%roots(perm(3))
        ulim_inv=dcmplx(1._dp/ulim)
        a2=ulim_inv-1._dp/x2
        a3=ulim_inv-1._dp/x3
        a4=ulim_inv-1._dp/x4
        xi=szc%xi
        s=TWOTHIRD*xi*szc%f(perm(3))/x4**2
        t1=THIRD*(1._dp/x2+1._dp/x3+1._dp/x4 &
                -3._dp*x4/(x2*x3))*rd(a2,a3,a4)
        t2=rf(a2,a3,a4)
        t3=-sqrt(a2)*sqrt(a3)/sqrt(a4)**3
        szellip_115_=s*(t1+t2+t3)
    end function szellip_115_
    
    real(dp) function szellip_3333(ulim,szc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the REAL integral from 0 to ulim of                      !
!                                                                  !
!    x^3/2 . (x-x1)^{-3/2} . (x-x2)^{-3/2} . (x-x3)^{-3/2} dx      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        complex(dp) :: s1
        complex(dp), dimension(3) :: ell,coeff
        ell(1)=szellip_1311(szc,ulim)
        ell(2)=szellip_1131(szc,ulim)
        ell(3)=szellip_1113(szc,ulim)
        coeff=szc%roots*szc%f
        s1=sum(coeff*ell)
        if (abs(aimag(s1)/s1) < 1.e-12_dp) then
            szellip_3333=real(s1)
        else
            STOP 'szellip_3333: result should be REAL'
        end if
    end function szellip_3333

    real(dp) function szellip_1333(ulim,szc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the REAL integral from 0 to ulim of                      !
!                                                                  !
!    x^1/2 . (x-x1)^{-3/2} . (x-x2)^{-3/2} . (x-x3)^{-3/2} dx      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: ulim
        type(szcmplx), intent(in) :: szc
        complex(dp) :: s1
        complex(dp), dimension(3) :: ell,coeff
        ell(1)=szellip_1311(szc,ulim)
        ell(2)=szellip_1131(szc,ulim)
        ell(3)=szellip_1113(szc,ulim)
        coeff=szc%f
        s1=sum(coeff*ell)
        if (abs(aimag(s1)/s1) < 1.e-12_dp) then
            szellip_1333=real(s1)
        else
            STOP 'szellip_1333: result should be REAL'
        end if
    end function szellip_1333
    
    real(dp) function szellip_5555(szc,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes  the elliptic integral [p]=[3,-5,-5,-5]   !
! by decomposing into 6 smaller integrals:           !
!                                                    !
! [p1,p2,p3,p4] = c1 . [1,-3,-1,-1]                  !
!               + c2 . [1,-1,-3,-1]                  !
!               + c3 . [1,-1,-1,-3]                  !
!               + c4 . [1,-5,-1,-1]                  !
!               + c5 . [1,-1,-5,-1]                  !
!               + c6 . [1,-1,-1,-5]                  !
! where                                              !
!                                                    !
!   [p1,p2,p3,p4] = integral from 0 to ulim of       !
!          x^p1/2 . (x-x2)^p2/2 . (x-x3)^p3/2        !
!                 . (x-x3)^p3/2 . (x-x4)^p4/2 dx     !
! and                                                !
!                                                    !
!    c1 = 2(x2.x3.x4 - x2^3)                         !
!                / ((x2 - x3) (x2 - x4))^3           !
!                                                    !
!    c2 = 2(x2.x3.x4 - x3^3)                         !
!                / ((x3 - x2) (x3 - x4))^3           !
!                                                    !
!    c3 = 2(x2.x3.x4 - x4^3)                         !
!                / ((x4 - x2) (x4 - x3))^3           !
!                                                    !
!    c4 = x2 / ((x2 - x3) (x2 - x4))^2               !
!                                                    !
!    c5 = x3 / ((x3 - x2) (x3 - x4))^2               !
!                                                    !
!    c6 = x4 / ((x4 - x2) (x4 - x3))^2               !
!                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
        complex(dp), dimension(6) :: c,y
        complex(dp) :: x2,x3,x4,prx,z
        real(dp), parameter :: ztol=1.e-12_dp
        x2 = szc%roots(1)
        x3 = szc%roots(2)
        x4 = szc%roots(3)
        prx = dcmplx(szc%prd_roots)
        c(1) = 2._dp * szc%f(1)**3 * (prx - x2**3)
        c(2) = 2._dp * szc%f(2)**3 * (prx - x3**3)
        c(3) = 2._dp * szc%f(3)**3 * (prx - x4**3)
        c(4) = (x2 * szc%f(1))**2
        c(5) = (x3 * szc%f(2))**2
        c(6) = (x4 * szc%f(3))**2
! Use simplified formula (faster)
        y(1) = c(1) * szellip_1311(szc,ulim)
        y(2) = c(2) * szellip_1131(szc,ulim)
        y(3) = c(3) * szellip_1113(szc,ulim)
        y(4) = c(4) * szellip_1511(szc,ulim)
        y(5) = c(5) * szellip_1151(szc,ulim)
        y(6) = c(6) * szellip_1115(szc,ulim)
        z = sum(y)
        if (aimag(z)/abs(z) < ztol) then
            szellip_5555 = real(z)
        else
            STOP 'szellip_5555: integral not real'
        end if
    end function szellip_5555

    real(dp) function szellip_3555(szc,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! **Differs from szellip_5555 in coefficients only** !
! Computes  the elliptic integral [p]=[3,-5,-5,-5]   !
! by decomposing into 6 smaller integrals:           !
!                                                    !
! [p1,p2,p3,p4] = c1 . [1,-3,-1,-1]                  !
!               + c2 . [1,-1,-3,-1]                  !
!               + c3 . [1,-1,-1,-3]                  !
!               + c4 . [1,-5,-1,-1]                  !
!               + c5 . [1,-1,-5,-1]                  !
!               + c6 . [1,-1,-1,-5]                  !
! where                                              !
!                                                    !
!   [p1,p2,p3,p4] = integral from 0 to ulim of       !
!          x^p1/2 . (x-x2)^p2/2 . (x-x3)^p3/2        !
!                 . (x-x3)^p3/2 . (x-x4)^p4/2 dx     !
! and                                                !
!                                                    !
!    c1 = (x2.x3 + x2.x4 + x3.x4 - 3.x2^3)           !
!                / ((x2 - x3) (x2 - x4))^3           !
!                                                    !
!    c2 = (x2.x3 + x2.x4 + x3.x4 - 3.x3^3)           !
!                / ((x3 - x2) (x3 - x4))^3           !
!                                                    !
!    c3 = (x2.x3 + x2.x4 + x3.x4 - 3.x4^3)           !
!                / ((x4 - x2) (x4 - x3))^3           !
!                                                    !
!    c4 = x2 / ((x2 - x3) (x2 - x4))^2               !
!                                                    !
!    c5 = x3 / ((x3 - x2) (x3 - x4))^2               !
!                                                    !
!    c6 = x4 / ((x4 - x2) (x4 - x3))^2               !
!                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
        complex(dp), dimension(6) :: c,y
        complex(dp) :: x2,x3,x4,prx,z        
        real(dp), parameter :: ztol=1.e-12_dp
        x2 = szc%roots(1)
        x3 = szc%roots(2)
        x4 = szc%roots(3)
        prx = x2 * x3 + x2 * x4 + x3 * x4
        c(1) = szc%f(1)**3 * (prx - 3._dp * x2**2)
        c(2) = szc%f(2)**3 * (prx - 3._dp * x3**2)
        c(3) = szc%f(3)**3 * (prx - 3._dp * x4**2)
        c(4) = x2 * szc%f(1)**2
        c(5) = x3 * szc%f(2)**2
        c(6) = x4 * szc%f(3)**2
! Use simplified formula (faster)
        y(1) = c(1) * szellip_1311(szc,ulim)
        y(2) = c(2) * szellip_1131(szc,ulim)
        y(3) = c(3) * szellip_1113(szc,ulim)
        y(4) = c(4) * szellip_1511(szc,ulim)
        y(5) = c(5) * szellip_1151(szc,ulim)
        y(6) = c(6) * szellip_1115(szc,ulim)
        z = sum(y)
        if (aimag(z)/abs(z) < ztol) then
            szellip_3555 = real(z)
        else
            STOP 'szellip_3555: integral not real'
        end if
    end function szellip_3555

    real(dp) function szellip_1555(szc,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! **Differs from szellip_5555 in coefficients only** !
! Computes  the elliptic integral [p]=[3,-5,-5,-5]   !
! by decomposing into 6 smaller integrals:           !
!                                                    !
! [p1,p2,p3,p4] = c1 . [1,-3,-1,-1]                  !
!               + c2 . [1,-1,-3,-1]                  !
!               + c3 . [1,-1,-1,-3]                  !
!               + c4 . [1,-5,-1,-1]                  !
!               + c5 . [1,-1,-5,-1]                  !
!               + c6 . [1,-1,-1,-5]                  !
! where                                              !
!                                                    !
!   [p1,p2,p3,p4] = integral from 0 to ulim of       !
!          x^p1/2 . (x-x2)^p2/2 . (x-x3)^p3/2        !
!                 . (x-x3)^p3/2 . (x-x4)^p4/2 dx     !
! and                                                !
!                                                    !
!    c1 = 2(x2 + x3 + x4 - 3.x2)                     !
!                / ((x2 - x3) (x2 - x4))^3           !
!                                                    !
!    c2 = 2(x2 + x3 + x4 - 3.x3)                     !
!                / ((x3 - x2) (x3 - x4))^3           !
!                                                    !
!    c3 = 2(x2 + x3 + x4 - 3.x4)                     !
!                / ((x4 - x2) (x4 - x3))^3           !
!                                                    !
!    c4 = x2 / ((x2 - x3) (x2 - x4))^2               !
!                                                    !
!    c5 = x3 / ((x3 - x2) (x3 - x4))^2               !
!                                                    !
!    c6 = x4 / ((x4 - x2) (x4 - x3))^2               !
!                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
        complex(dp), dimension(6) :: c,y
        complex(dp) :: x2,x3,x4,sumx,z
        real(dp), parameter :: ztol=1.e-12_dp
        x2 = szc%roots(1)
        x3 = szc%roots(2)
        x4 = szc%roots(3)
        sumx = x2 + x3 + x4
        c(1) = szc%f(1)**3 * (sumx - 3._dp * x2) * 2._dp
        c(2) = szc%f(2)**3 * (sumx - 3._dp * x3) * 2._dp
        c(3) = szc%f(3)**3 * (sumx - 3._dp * x4) * 2._dp
        c(4) = szc%f(1)**2
        c(5) = szc%f(2)**2
        c(6) = szc%f(3)**2
! Use simplified formula (faster)
        y(1) = c(1) * szellip_1311(szc,ulim)
        y(2) = c(2) * szellip_1131(szc,ulim)
        y(3) = c(3) * szellip_1113(szc,ulim)
        y(4) = c(4) * szellip_1511(szc,ulim)
        y(5) = c(5) * szellip_1151(szc,ulim)
        y(6) = c(6) * szellip_1115(szc,ulim)
        z = sum(y)
        if (aimag(z)/abs(z) < ztol) then
            szellip_1555 = real(z)
        else
            STOP 'szellip_1555: integral not real'
        end if
        
    end function szellip_1555    

end module szelliptic
