module szlocal
  use precision1, only : dp
  use constants
  use cosmo_params, only : H0,ct2,FLRW,r0,dlr,get_age
  use szclass
  use szelliptic
  use szfuncs, only : k_approx,Omr0,Mfunc,Mprime, &
                      Mpprime
  use utils, only : fcubic_roots,splint,spline, &
                    zbrent
  implicit none

  real(dp), parameter :: OH2=FLRW%Ol*H0**2
  real(dp), parameter :: C0=1._dp/OH2
    
  type szshell
     real(dp) :: r
     real(dp) :: M
     real(dp) :: k
     type(szcmplx) :: szc
  end type szshell
  
  type(spline_data), pointer :: kspl => null()
!$omp threadprivate(kspl)

contains

    subroutine init_model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate arrays from last call to init_model then !
! recompute age and splines with new void parameters  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer, parameter :: npts=700
        call get_age !Mpc
        r0=r0/FLRW%h
        dlr=0.1_dp*r0
        if (associated(kspl)) then
            call del_spline_data(kspl)
            deallocate(kspl)
        end if
        allocate(kspl)
        call etbcubic(npts,kspl) !get & store splines
    end subroutine init_model
        
    subroutine init_shell(r,shell)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Must be called for each new r before below routines can   !
! be called. Shells labelled only by coord r (time not req) !
! from which p,q,roots can be determined                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: r
        type(szshell), intent(out) :: shell
        type(szcmplx) :: szc
        real(dp) :: M,k,p,q
        call get_k(r,k)
        M=Mfunc(r)
        p=-k*C0
        q=2._dp*M*C0
        call init_szcmplx(p,q,szc)
        shell=szshell(r,M,k,szc)
    end subroutine init_shell

    subroutine get_k(r,k)
        implicit none
        real(dp), intent(in) :: r
        real(dp), intent(out) :: k
        real(dp), parameter :: acc=1.e-9_dp
        real(dp), parameter :: eps=0.2_dp
        real(dp) :: ul,p,q,pmid,plo,phi
        real(dp), save :: plast
!$omp threadprivate(plast)
        ul=r
        q=2._dp*Mfunc(r)*C0
        pmid=-k_approx(r)*C0
        if (abs(pmid) < 1.e-16) pmid=plast
        plo=pmid-eps*abs(pmid)
        phi=pmid+eps*abs(pmid)
        p=zbrent(p_eqn,plo,phi,acc)
        plast=p
        k=-p*OH2
    contains
        real(dp) function p_eqn(pin)
! f(p) = c*t0 - [1,-1,-1,-1]/sqrt(Omegal*H0^2)
            implicit none
            real(dp), intent(in) :: pin
            complex(dp), dimension(3) :: roots
            call fcubic_roots(pin,q,roots)
            p_eqn=ct2-sqrt(C0)*szellip_1111(ul,roots)
        end function p_eqn
    end subroutine get_k
    
    subroutine get_kprime(shell,kprime)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(out) :: kprime
        real(dp) :: C1,e2,e3,Mp,Rd,f1,f2
        real(dp), parameter :: t=0._dp
        C1=sqrt(C0)*C0
        e2=C1*szellip_3333(shell%szc,shell%r)
        e3=C1*szellip_1333(shell%szc,shell%r)
        Mp=Mprime(shell%r)
        call get_Rdot(shell,t,Rd)
        f1=2._dp/e2
        f2=e3*Mp-1._dp/Rd
        kprime=f1*f2
    end subroutine get_kprime
    
    subroutine etbcubic(npts,spl)
        implicit none
        integer, intent(in) :: npts
        type(spline_data), intent(out) :: spl
        real(dp), dimension(:), allocatable :: rarr
        real(dp), dimension(:), allocatable :: karr
        real(dp), dimension(:), allocatable :: k2arr
        real(dp), parameter :: tol=1.0e-9
        real(dp), dimension(5) :: ro
        real(dp) :: h,kp1,kpn,r1,rn,ri
        integer :: i,m,n
        kp1=0._dp
        kpn=0._dp
        !origin needs to be handled separately:
        ro=(/1.0e-5_dp, 1.0e-4_dp, &
             1.0e-3_dp, 1.0e-2_dp, 0.1_dp/)
        allocate(rarr(npts))
        allocate(karr(npts))
        m=npts-size(ro)
        r1=1._dp
        rn=5._dp*r0
        h=(rn-r1)/(m-1)
        rarr(1:size(ro))=ro
        rarr(size(ro)+1:npts)=(/(r1+h*(i-1),i=1,m)/)
        n=0
        do i=1,npts
            ri=rarr(i)
            call get_k(ri,karr(i))
            if (abs(Omr0(ri)-FLRW%Om) < tol) EXIT
            n=i
        end do
        allocate(k2arr(n))
        call spline(rarr(1:n),karr(1:n),n,kp1,kpn,k2arr)
        call init_spline_data(n,spl)
        spl%n=n
        spl%x=rarr(1:n)
        spl%f=karr(1:n)
        spl%ddf=k2arr(1:n)
        deallocate(k2arr)
        deallocate(karr)
        deallocate(rarr)
    end subroutine etbcubic

    subroutine get_R(shell,t,R)
        implicit none
        type(szshell), intent(in) :: shell      
        real(dp), intent(in) :: t
        real(dp), intent(out) :: R
        real(dp), parameter :: acc=1.e-9_dp
        real(dp), parameter :: h=0.2_dp
        real(dp) :: SQRTC0,t1,Rlo,Rhi
        if (abs(t) < 1.e-14_dp) then
            R=shell%r
            return
        end if
        SQRTC0=sqrt(C0)
! t0=0 so shift time coordinate so t0=age_of_univ
        t1=ct2+t
        Rlo=shell%r*(1._dp-h)
        Rhi=shell%r
        R=zbrent(R_eqn,Rlo,Rhi,acc)
    contains
        real(dp) function R_eqn(Rin)
            implicit none
            real(dp), intent(in) :: Rin
            real(dp) ul
            ul=Rin
            R_eqn=t1-SQRTC0*szellip_1111(ul,shell%szc)
        end function R_eqn
    end subroutine get_R

    subroutine get_Rdot(shell,t,Rdot,R_)
        implicit none
        type(szshell), intent(in) :: shell 
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Rdot
        real(dp), intent(in), optional :: R_
        real(dp) :: R,p,q,s2,Dt
        if (present(R_)) then; R=R_; else; call get_R(shell,t,R); end if
        p=shell%szc%p
        q=shell%szc%q
        s2=R**2+p+q/R
        if (s2 < 0._dp) then
            Dt=-(4._dp*p**3+27._dp*q**2)
            write(*,'(A,X,A,E11.2,X,A,X,A)') &
                  'get_Rdot: Rdot^2 is negative.', &
                  'Discriminant is', Dt, &
                  'Cubic has more than one real root', &
                  'when discriminant is positive'
            STOP
        end if
        Rdot=sqrt(OH2)*sqrt(s2)
    end subroutine get_Rdot

    subroutine get_Rddot(shell,t,Rddot,R_)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Rddot
        real(dp), intent(in), optional :: R_
        real(dp) :: R,M
        if (present(R_)) then; R=R_; else; call get_R(shell,t,R); end if
        M=shell%M
        Rddot=-M/R**2+OH2*R
    end subroutine get_Rddot
    
    subroutine get_Rprime(shell,t,Rprime,kp_,R_,Rd_)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Rprime
        real(dp), intent(in), optional :: kp_
        real(dp), intent(in), optional :: R_
        real(dp), intent(in), optional :: Rd_
        real(dp) :: C1,R,Rd,kp,Mp,e2,e3
        if (abs(t) < 1.e-14_dp) then
            Rprime=1._dp
            return
        end if
        if (present(kp_)) then; kp=kp_; else; call get_kprime(shell,kp); end if
        if (present(R_)) then; R=R_; else; call get_R(shell,t,R); end if
        if (present(Rd_)) then; Rd=Rd_; else; call get_Rdot(shell,t,Rd,R); end if
        C1=sqrt(C0)*C0
        Mp=Mprime(shell%r)
        e2=C1*szellip_3333(shell%szc,ulim=R)
        e3=C1*szellip_1333(shell%szc,ulim=R)
        Rprime=Rd*(Mp*e3-0.5_dp*kp*e2)
    end subroutine get_Rprime

    subroutine get_Rdotprime(shell,t,Rdotprime,kp_,R_,Rd_,Rp_)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Rdotprime
        real(dp), intent(in), optional :: kp_
        real(dp), intent(in), optional :: R_
        real(dp), intent(in), optional :: Rd_
        real(dp), intent(in), optional :: Rp_      
        real(dp) :: R,Rd,Rp,kp,M,Mp
        real(dp) :: f1,f2,g1,g2
        if (present(kp_)) then
            kp=kp_; else; call get_kprime(shell,kp)
        end if
        if (present(R_)) then
            R=R_; else; call get_R(shell,t,R)
        end if
        if (present(Rd_)) then
            Rd=Rd_; else; call get_Rdot(shell,t,Rd,R)
        end if
        if (present(Rp_)) then
            Rp=Rp_; else; call get_Rprime(shell,t,Rp,kp,R,Rd)
        end if
        M=shell%M
        Mp=Mprime(shell%r)
        f1=1._dp/Rd
        f2=Mp/R-0.5_dp*kp
        g1=Rp*f1
        g2=OH2*R-M/R**2
        Rdotprime=f1*f2+g1*g2
    end subroutine get_Rdotprime

    subroutine get_kpprime(shell,kpprime,kp_)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(out) :: kpprime
        real(dp), intent(in), optional :: kp_
        real(dp) :: C1,C2,t,Rd,Rdp
        real(dp) :: kp,Mp,Mpp,e2,e3,e4
        real(dp) :: e5,e6,d2,d3,f1,f2
        if (present(kp_)) then
            kp=kp_
        else
            call get_kprime(shell,kp)
        end if
        C1=sqrt(C0)*C0
        C2=C1*C0
        t=0._dp
        call get_Rdot(shell,t,Rd)
        call get_Rdotprime(shell,t,Rdp)
        e2=C1*szellip_3333(shell%szc,ulim=shell%r)
        e3=C1*szellip_1333(shell%szc,ulim=shell%r)
        e4=C2*szellip_5555(shell%szc,ulim=shell%r)
        e5=C2*szellip_3555(shell%szc,ulim=shell%r)
        e6=C2*szellip_1555(shell%szc,ulim=shell%r)
        Mp=Mprime(shell%r)
        Mpp=Mpprime(shell%r)
        d2=1._dp/Rd**3-3._dp*Mp*e5+1.5_dp*kp*e4
        d3=1._dp/(shell%r*Rd**3)-3._dp*Mp*e6+1.5_dp*kp*e5
        f1=-(e3*Mp-1._dp/Rd)*d2/e2
        f2=d3*Mp+e3*Mpp+Rdp/Rd**2
        kpprime=2._dp*(f1+f2)/e2
    end subroutine get_kpprime

    subroutine get_Rpprime(shell,t,Rpprime,kp_,R_,Rd_,Rp_,Rdp_,kpp_)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Rpprime
        real(dp), intent(in), optional :: kp_
        real(dp), intent(in), optional :: R_
        real(dp), intent(in), optional :: Rd_
        real(dp), intent(in), optional :: Rp_
        real(dp), intent(in), optional :: Rdp_
        real(dp), intent(in), optional :: kpp_
        real(dp) :: C1,C2,Mp,Mpp
        real(dp) :: e2,e3,e4,e5,e6,d2,d3,f1,f2
        real(dp) :: kp,R,Rd,Rp,Rdp,kpp,Rp_Rd3
        real(dp), parameter :: HALF=0.5_dp,THREE=3._dp
        if (abs(t) < 1.e-14_dp) then
            Rpprime=0._dp
            return
        end if
        if (present(kp_)) then
            kp=kp_; else; call get_kprime(shell,kp)
        end if
        if (present(R_)) then
            R=R_; else; call get_R(shell,t,R)
        end if
        if (present(Rd_)) then
            Rd=Rd_; else; call get_Rdot(shell,t,Rd,R)
        end if
        if (present(Rp_)) then
            Rp=Rp_; else; call get_Rprime(shell,t,Rp,kp,R,Rd)
        end if
        if (present(Rdp_)) then
            Rdp=Rdp_; else; call get_Rdotprime(shell,t,Rdp,kp,R,Rd,Rp)
        end if
        if (present(kpp_)) then
            kpp=kpp_; else; call get_kpprime(shell,kpp,kp)
        end if
        C1=sqrt(C0)*C0
        C2=C1*C0
        Mp=Mprime(shell%r)
        Mpp=Mpprime(shell%r)
        Rp_Rd3=Rp/Rd**3
        e2=C1*szellip_3333(shell%szc,ulim=R)
        e3=C1*szellip_1333(shell%szc,ulim=R)
        e4=C2*szellip_5555(shell%szc,ulim=R)
        e5=C2*szellip_3555(shell%szc,ulim=R)
        e6=C2*szellip_1555(shell%szc,ulim=R)
        d2=Rp_Rd3-THREE*(Mp*e5-HALF*kp*e4)
        d3=Rp_Rd3/R-THREE*(Mp*e6-HALF*kp*e5)
        f1=-HALF*e2*(kp*Rdp+kpp*Rd)
        f2=e3*(Mp*Rdp+Mpp*Rd)
        Rpprime=f1+f2+Rd*(Mp*d3-HALF*kp*d2)
    end subroutine get_Rpprime
    
    subroutine Hubble(shell,t,H)
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: t
        real(dp), intent(out) :: H
        real(dp) :: R,Rd,Rp,Rdp
        call get_R(shell,t,R)
        call get_Rdot(shell,t,Rd)
        call get_Rprime(shell,t,Rp)
        call get_Rdotprime(shell,t,Rdp)
        H=THIRD*(2._dp*Rd/R+Rdp/Rp)
    end subroutine Hubble
    
end module szlocal
