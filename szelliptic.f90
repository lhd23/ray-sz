module szelliptic
  use precision1, only : dp
  use constants
  use utils, only : fcubic_roots,zbrent, &
                    rj,rd,rf,rc
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
  public :: szcmplx,init_szcmplx,szellip_1111, &
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

    real(dp) function szellip_1111_roots(roots,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the integral from 0 to ulim of                          !
!                                                                 !
!    x^1/2 . (x-x1)^{-1/2} . (x-x2)^{-1/2} . (x-x3)^{-1/2} dx     !
!                                                                 !
! normalised by 1 / sqrt(Omegal * H0^2)                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        complex(dp), dimension(3), &
                intent(in) :: roots
        real(dp), intent(in) :: ulim
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

    real(dp) function szellip_1111(szc,ulim)
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
        complex(dp), dimension(3) :: x
        complex(dp) :: ulim_inv
        real(dp) :: f
        ulim_inv=dcmplx(1._dp/ulim)
        x=ulim_inv-szc%roots_inv
        f=TWOTHIRD*szc%p3
        szellip_1111=f*real(rj(x(1),x(2),x(3),ulim_inv))
    end function szellip_1111

    real(dp) function I3prime(a1,a4,b1,b4,f,g,h,x,y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the integral in terms of Carlson R-functions of             !
!                                                                     !
!    I3'=[1,-1,-1,-1]                                                 !
!       =\int_y^x (f+g t+h t^2) (a_1+b_1 t)^{-1/2} (a_4+b_4 t)^{-1/2} !
!                                                                     !
! where f+g t+ h t^2 > 0 and x-y > 0    [Carlson (1991) eq (2.17)]    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: a1,a4,b1,b4
        real(dp), intent(in) :: f,g,h
        real(dp), intent(in) :: x,y
        real(dp), parameter :: TWO=2._dp,THREE=3._dp,FOUR=4._dp
        real(dp), parameter :: SIX=6._dp,NINE=9._dp
        real(dp) :: X1,X4,Y1,Y4,d14,beta1
        real(dp) :: c112,c442,c142,c11,c44,xi,eta
        real(dp) :: M2,Lp2,Lm2,U,U2,W12
        real(dp) :: Q12,P12,rho,rj_,rf_,rc_
! Carlson 1991 (2.1)
        X1=sqrt(a1+b1*x)
        X4=sqrt(a4+b4*x)
        Y1=sqrt(a1+b1*y)
        Y4=sqrt(a4+b4*y)
        d14=a1*b4-a4*b1
! Carlson 1991 (2.2)
        beta1=g*b1-TWO*h*a1
! Carlson 1991 (2.3)
        c112=TWO*f*b1*b1-g*(a1*b1+a1*b1)+TWO*h*a1*a1
        c442=TWO*f*b4*b4-g*(a4*b4+a4*b4)+TWO*h*a4*a4
        c142=TWO*f*b1*b4-g*(a1*b4+a4*b1)+TWO*h*a1*a4
        c11=sqrt(c112)
        c44=sqrt(c442)
! (2.4)
        xi=sqrt(f+x*(g+h*x))
        eta=sqrt(f+y*(g+h*y))
! (2.6) xi,eta,x,y,f,g,h all > 0
        M2=(TWO*xi*eta+TWO*f+g*(x+y)+TWO*h*x*y) &
            *((X1*Y4+Y1*X4)/(x-y))**2
! (2.7)
        Lp2=M2+c142+c11*c44
        Lm2=M2+c142-c11*c44
! (2.8)
        U=(X1*X4*eta+Y1*Y4*xi)/(x-y)
        U2=U**2
! (2.10)
        W12=U2-c112*b4/(TWO*b1)
        Q12=W12/(X1*Y1)**2
        P12=Q12+h*b4/b1
! (2.11)
        rho=d14*(beta1-sqrt(TWO*h)*c11)/b1
        rj_=FOUR*rho*rj(M2,Lm2,Lp2,M2+rho)
        rf_=-SIX*rf(M2,Lm2,Lp2)
        rc_=THREE*rc(U2,W12)
        I3prime=sqrt(TWO*c112/(NINE*h))*(rj_+rf_+rc_) &
                +TWO*rc(P12,Q12)
    end function I3prime

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
    
    real(dp) function szellip_3333(szc,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the REAL integral from 0 to ulim of                      !
!                                                                  !
!    x^3/2 . (x-x1)^{-3/2} . (x-x2)^{-3/2} . (x-x3)^{-3/2} dx      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
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

    real(dp) function szellip_1333(szc,ulim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the REAL integral from 0 to ulim of                      !
!                                                                  !
!    x^1/2 . (x-x1)^{-3/2} . (x-x2)^{-3/2} . (x-x3)^{-3/2} dx      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        type(szcmplx), intent(in) :: szc
        real(dp), intent(in) :: ulim
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
        x2=szc%roots(1)
        x3=szc%roots(2)
        x4=szc%roots(3)
        prx=dcmplx(szc%prd_roots)
        c(1)=2._dp*szc%f(1)**3*(prx-x2**3)
        c(2)=2._dp*szc%f(2)**3*(prx-x3**3)
        c(3)=2._dp*szc%f(3)**3*(prx-x4**3)
        c(4)=(x2*szc%f(1))**2
        c(5)=(x3*szc%f(2))**2
        c(6)=(x4*szc%f(3))**2
! Use simplified formula (faster)
        y(1)=c(1)*szellip_1311(szc,ulim)
        y(2)=c(2)*szellip_1131(szc,ulim)
        y(3)=c(3)*szellip_1113(szc,ulim)
        y(4)=c(4)*szellip_1511(szc,ulim)
        y(5)=c(5)*szellip_1151(szc,ulim)
        y(6)=c(6)*szellip_1115(szc,ulim)
        z=sum(y)
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
        x2=szc%roots(1)
        x3=szc%roots(2)
        x4=szc%roots(3)
        prx=x2*x3+x2*x4+x3*x4
        c(1)=szc%f(1)**3*(prx-3._dp*x2**2)
        c(2)=szc%f(2)**3*(prx-3._dp*x3**2)
        c(3)=szc%f(3)**3*(prx-3._dp*x4**2)
        c(4)=x2*szc%f(1)**2
        c(5)=x3*szc%f(2)**2
        c(6)=x4*szc%f(3)**2
! Use simplified formula (faster)
        y(1)=c(1)*szellip_1311(szc,ulim)
        y(2)=c(2)*szellip_1131(szc,ulim)
        y(3)=c(3)*szellip_1113(szc,ulim)
        y(4)=c(4)*szellip_1511(szc,ulim)
        y(5)=c(5)*szellip_1151(szc,ulim)
        y(6)=c(6)*szellip_1115(szc,ulim)
        z=sum(y)
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
        x2=szc%roots(1)
        x3=szc%roots(2)
        x4=szc%roots(3)
        sumx=x2+x3+x4
        c(1)=szc%f(1)**3*(sumx-3._dp*x2)*2._dp
        c(2)=szc%f(2)**3*(sumx-3._dp*x3)*2._dp
        c(3)=szc%f(3)**3*(sumx-3._dp*x4)*2._dp
        c(4)=szc%f(1)**2
        c(5)=szc%f(2)**2
        c(6)=szc%f(3)**2
! Use simplified formula (faster)
        y(1)=c(1)*szellip_1311(szc,ulim)
        y(2)=c(2)*szellip_1131(szc,ulim)
        y(3)=c(3)*szellip_1113(szc,ulim)
        y(4)=c(4)*szellip_1511(szc,ulim)
        y(5)=c(5)*szellip_1151(szc,ulim)
        y(6)=c(6)*szellip_1115(szc,ulim)
        z=sum(y)
        if (aimag(z)/abs(z) < ztol) then
            szellip_1555 = real(z)
        else
            STOP 'szellip_1555: integral not real'
        end if
        
    end function szellip_1555    

end module szelliptic
