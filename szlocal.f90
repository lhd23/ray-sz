! A module for computing the radial and temporal functions
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

  type szderivs
      real(dp), pointer :: R   => null()
      real(dp), pointer :: Rp  => null()
      real(dp), pointer :: Rd  => null()
      real(dp), pointer :: Rpp => null()
      real(dp), pointer :: Rdd => null()
      real(dp), pointer :: Rdp => null()
      real(dp), pointer :: k   => null()
      real(dp), pointer :: kp  => null()
      real(dp), pointer :: kpp => null()
      real(dp), pointer :: M   => null()
  end type szderivs

  type szlocal_class
      real(dp) :: t
      real(dp) :: r
      type(szderivs) :: f
      type(szcmplx) :: szc
  end type szlocal_class
  
  type(spline_data), pointer :: kspl => null()
!$omp threadprivate(kspl)

contains

    subroutine init_model(interp) !init global model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate arrays from last call to init_model then !
! recompute age and splines with new void parameters  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        logical, optional, intent(in) :: interp
        logical do_interp
        integer, parameter :: npts=700
        if (present(interp)) then
            do_interp=interp
        else
            do_interp=.false.
        end if
        call get_age !Mpc
        r0=r0/FLRW%h
        dlr=0.1_dp*r0
        if (do_interp) then
            if (associated(kspl)) then
                call del_spline_data(kspl)
                deallocate(kspl)
            else !get & store splines  
                allocate(kspl)
                call etbcubic(npts,kspl)
            end if
        end if
    end subroutine init_model
        
    subroutine init_shell(r,shell) !legacy
        implicit none
        real(dp), intent(in) :: r
        type(szshell), intent(out) :: shell
        type(szcmplx) :: szc
        real(dp) :: M,k,p,q 
        k=get_k(r)
        M=Mfunc(r)
        p=-k*C0
        q=2._dp*M*C0
        call init_szcmplx(p,q,szc)
        shell=szshell(r,M,k,szc)
    end subroutine init_shell

    function init_szlocal_class(t,r) result(ctr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For each new t and r this routine must be called          !
! before below routines can be called                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(in) :: r
        type(szlocal_class) :: ctr
        type(szcmplx) :: szc
        real(dp) :: p,q,k
        if (r < 0.) write(*,*) 'warning: r is negative'
        ctr%t=t
        ctr%r=r
        allocate(ctr%f%k)
        allocate(ctr%f%M)
        ! call splint(kspl%x,kspl%f,kspl%ddf,kspl%n,r,k)
        ! ctr%f%k=k
        ctr%f%k=get_k(r)
        ctr%f%M=Mfunc(r)
        p=-ctr%f%k*C0
        q=2._dp*ctr%f%M*C0
        call init_szcmplx(p,q,ctr%szc)
    end function init_szlocal_class

    function get_k(r) result(k)
        implicit none
        real(dp), intent(in) :: r
        real(dp), parameter :: acc=1.e-9_dp
        real(dp), parameter :: eps=0.2_dp
        real(dp) :: k,p,q,pmid,plo,phi
        real(dp), save :: plast
!$omp threadprivate(plast)
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
            p_eqn=ct2-sqrt(C0)*szellip_1111(roots,ulim=r)
        end function p_eqn
    end function get_k
    
    subroutine get_kprime(szloc)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp) :: C1,e2,e3,Mp,Rd,f1,f2
        if (associated(szloc%f%kp)) then
            return
        else
            allocate(szloc%f%kp)
            C1=sqrt(C0)*C0
            e2=C1*szellip_3333(szloc%szc,ulim=szloc%r)
            e3=C1*szellip_1333(szloc%szc,ulim=szloc%r)
            call get_Rdot(szloc,Rd)
            Mp=Mprime(szloc%r)
            f1=2._dp/e2
            f2=e3*Mp-1._dp/Rd
            szloc%f%kp=f1*f2
            return
        end if
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
            karr(i)=get_k(ri)
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

    subroutine get_R(szloc)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp), parameter :: acc=1.e-9_dp
        real(dp), parameter :: h=0.2_dp
        real(dp) :: SQRTC0,t1,Rlo,Rhi
        if (associated(szloc%f%R)) then
            return
        else
            allocate(szloc%f%R)
            if (abs(szloc%t) < 1.e-14_dp) then
                szloc%f%R=szloc%r
                return
            end if
            SQRTC0=sqrt(C0)
! t0=0 so shift time coordinate so t0=age_of_univ
            t1=ct2+szloc%t
            Rlo=szloc%r*(1._dp-h)
            Rhi=szloc%r
            szloc%f%R=zbrent(R_eqn,Rlo,Rhi,acc)
            return
        end if
    contains
        real(dp) function R_eqn(Rin)
            implicit none
            real(dp), intent(in) :: Rin
            R_eqn=t1-SQRTC0*szellip_1111(szloc%szc,ulim=Rin)
        end function R_eqn
    end subroutine get_R

    subroutine get_Rdot(szloc,Rd0)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp), intent(out), optional :: Rd0
        real(dp) :: R,p,q,s2,Dt,Rd
! R will either be r if t=t0 or R=R(t,r)
        if (present(Rd0)) then !compute Rdot(t=t0,r) 
            R=szloc%r
        else !compute Rdot(t,r) where t=szloc%t
            if (associated(szloc%f%Rd)) then
                return !no need to continue
            else
                allocate(szloc%f%Rd)
                if (.not. associated(szloc%f%R)) &
                        call get_R(szloc)
                R=szloc%f%R
            end if
        end if
        p=szloc%szc%p
        q=szloc%szc%q
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
        Rd=sqrt(OH2)*sqrt(s2)
        if (present(Rd0)) then
            Rd0=Rd
        else
            szloc%f%Rd=Rd
        end if
    end subroutine get_Rdot

    subroutine get_Rddot(szloc)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp) :: R,M
        if (associated(szloc%f%Rdd)) then
            return
        else
            allocate(szloc%f%Rdd)
            if (.not. associated(szloc%f%R)) &
                    call get_R(szloc)
            R=szloc%f%R
            M=szloc%f%M
            szloc%f%Rdd=-M/R**2+OH2*R
            return
        end if
    end subroutine get_Rddot
    
    subroutine get_Rprime(szloc,Rp0)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp), intent(out), optional :: Rp0
        real(dp) :: C1,R,Rd,kp,Mp,e2,e3
        if (present(Rp0)) then
            Rp0=1._dp
            return
        else
            if (associated(szloc%f%Rp)) then
                return
            else
                allocate(szloc%f%Rp)
                if (abs(szloc%t) < 1.e-14_dp) then
                    szloc%f%Rp=1._dp
                    return
                end if
                if (.not. associated(szloc%f%kp)) &
                        call get_kprime(szloc)
                if (.not. associated(szloc%f%R)) &
                        call get_R(szloc)
                if (.not. associated(szloc%f%Rd)) &
                        call get_Rdot(szloc)
! note it is important get_R is called before get_Rdot 
! for Rdot can be quickly evaluated once R is gotten
                kp=szloc%f%kp
                R=szloc%f%R
                Rd=szloc%f%Rd
                C1=sqrt(C0)*C0
                Mp=Mprime(szloc%r)
                e2=C1*szellip_3333(szloc%szc,ulim=R)
                e3=C1*szellip_1333(szloc%szc,ulim=R)
                szloc%f%Rp=Rd*(Mp*e3-0.5_dp*kp*e2)
                return
            end if
        end if
    end subroutine get_Rprime

    subroutine get_Rdotprime(szloc,Rdp0)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp), intent(out), optional :: Rdp0
        real(dp) :: R,Rd,Rp,kp,M,Mp,Rdp
        real(dp) :: f1,f2,g1,g2
        if (present(Rdp0)) then
            if (.not. associated(szloc%f%kp)) &
                    call get_kprime(szloc)
            kp=szloc%f%kp
            R=szloc%r
            call get_Rdot(szloc,Rd)
            call get_Rprime(szloc,Rp)
        else
            if (associated(szloc%f%Rdp)) then
                return
            else
                allocate(szloc%f%Rdp)
                if (.not. associated(szloc%f%kp)) &
                        call get_kprime(szloc)
                if (.not. associated(szloc%f%R)) &
                        call get_R(szloc)
                if (.not. associated(szloc%f%Rd)) &
                        call get_Rdot(szloc)
                if (.not. associated(szloc%f%Rp)) &
                        call get_Rprime(szloc)
                kp=szloc%f%kp
                R=szloc%f%R
                Rd=szloc%f%Rd
                Rp=szloc%f%Rp
            end if
        end if
        M=szloc%f%M
        Mp=Mprime(szloc%r)
        f1=1._dp/Rd
        f2=Mp/R-0.5_dp*kp
        g1=Rp*f1
        g2=OH2*R-M/R**2
        Rdp=f1*f2+g1*g2
        if (present(Rdp0)) then
            Rdp0=Rdp
        else
            szloc%f%Rdp=Rdp
        end if
    end subroutine get_Rdotprime

    subroutine get_kpprime(szloc)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp) :: C1,C2,t,Rd,Rdp
        real(dp) :: kp,Mp,Mpp,e2,e3,e4
        real(dp) :: e5,e6,d2,d3,f1,f2
        if (associated(szloc%f%kpp)) then
            return
        else
            allocate(szloc%f%kpp)
            if (.not. associated(szloc%f%kp)) &
                    call get_kprime(szloc)
            kp=szloc%f%kp
            C1=sqrt(C0)*C0
            C2=C1*C0
            call get_Rdot(szloc,Rd)
            call get_Rdotprime(szloc,Rdp)
            e2=C1*szellip_3333(szloc%szc,ulim=szloc%r)
            e3=C1*szellip_1333(szloc%szc,ulim=szloc%r)
            e4=C2*szellip_5555(szloc%szc,ulim=szloc%r)
            e5=C2*szellip_3555(szloc%szc,ulim=szloc%r)
            e6=C2*szellip_1555(szloc%szc,ulim=szloc%r)
            Mp=Mprime(szloc%r)
            Mpp=Mpprime(szloc%r)
            d2=1._dp/Rd**3-3._dp*Mp*e5+1.5_dp*kp*e4
            d3=1._dp/(szloc%r*Rd**3)-3._dp*Mp*e6+1.5_dp*kp*e5
            f1=-(e3*Mp-1._dp/Rd)*d2/e2
            f2=d3*Mp+e3*Mpp+Rdp/Rd**2
            szloc%f%kpp=2._dp*(f1+f2)/e2
            return
        end if
    end subroutine get_kpprime

    subroutine get_Rpprime(szloc)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp) :: C1,C2,Mp,Mpp
        real(dp) :: e2,e3,e4,e5,e6,d2,d3,f1,f2
        real(dp) :: kp,R,Rd,Rp,Rdp,kpp,Rp_Rd3
        real(dp), parameter :: HALF=0.5_dp,THREE=3._dp
        if (associated(szloc%f%Rpp)) then
            return
        else
            allocate(szloc%f%Rpp)
            if (abs(szloc%t) < 1.e-14_dp) then
                szloc%f%Rpp=0._dp
                return
            end if
            if (.not. associated(szloc%f%kp)) &
                    call get_kprime(szloc)
            if (.not. associated(szloc%f%R)) &
                    call get_R(szloc)
            if (.not. associated(szloc%f%Rd)) &
                    call get_Rdot(szloc)
            if (.not. associated(szloc%f%Rp)) &
                    call get_Rprime(szloc)
            if (.not. associated(szloc%f%Rdp)) &
                    call get_Rdotprime(szloc)
            if (.not. associated(szloc%f%kpp)) &
                    call get_kpprime(szloc)
! call above functions in this order is most efficient
            kp=szloc%f%kp
            R=szloc%f%R
            Rd=szloc%f%Rd
            Rp=szloc%f%Rp
            Rdp=szloc%f%Rdp
            kpp=szloc%f%kpp
            C1=sqrt(C0)*C0
            C2=C1*C0
            Mp=Mprime(szloc%r)
            Mpp=Mpprime(szloc%r)
            Rp_Rd3=Rp/Rd**3
            e2=C1*szellip_3333(szloc%szc,ulim=R)
            e3=C1*szellip_1333(szloc%szc,ulim=R)
            e4=C2*szellip_5555(szloc%szc,ulim=R)
            e5=C2*szellip_3555(szloc%szc,ulim=R)
            e6=C2*szellip_1555(szloc%szc,ulim=R)
            d2=Rp_Rd3-THREE*(Mp*e5-HALF*kp*e4)
            d3=Rp_Rd3/R-THREE*(Mp*e6-HALF*kp*e5)
            f1=-HALF*e2*(kp*Rdp+kpp*Rd)
            f2=e3*(Mp*Rdp+Mpp*Rd)
            szloc%f%Rpp=f1+f2+Rd*(Mp*d3-HALF*kp*d2)
            return
        end if
    end subroutine get_Rpprime
    
    subroutine Hubble(szloc,H)
        implicit none
        type(szlocal_class), intent(inout) :: szloc
        real(dp), intent(out) :: H
        real(dp) :: R,Rd,Rp,Rdp
        call get_R(szloc)
        call get_Rdot(szloc)
        call get_Rprime(szloc)
        call get_Rdotprime(szloc)
        R=szloc%f%R
        Rd=szloc%f%Rd
        Rp=szloc%f%Rp
        Rdp=szloc%f%Rdp
        H=THIRD*(2._dp*Rd/R+Rdp/Rp)
    end subroutine Hubble

end module szlocal
