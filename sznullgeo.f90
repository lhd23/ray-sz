module sznullgeo
  use precision1, only : dp
  use constants
  use szgeom
  use pix_tools, only : ang2vec
  implicit none

  private

  public :: rayshoot_dr
  
contains

    subroutine rayshoot_dr(initial_pos,RA,DEC,yout,iexit,val_exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The exit options are                                                        !
!                                                                             !
!      iexit = 1    compute luminosity distance out to a fixed redshift       !
!      iexit = 2    compute redshift out to a fixed luminosity distance       !
!      otherwise    compute redshift out to a fixed comoving time in the past !
!                                                                             !
! the fixed value is specified by val_exit                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in)  :: initial_pos(3)
        real(dp), intent(in)  :: RA,DEC
        real(dp), intent(out) :: yout
        integer,  intent(in)  :: iexit
        real(dp), intent(in)  :: val_exit
        real(dp), dimension(4) :: pvi,nvi
        real(dp), dimension(:), allocatable :: yi,yf
        real(dp) :: kti,ktf,si,ds,t_init,zp1
        t_init=0._dp
        si=0._dp
        ds=1.e-2_dp
        call construct_ray(t_init,initial_pos,RA,DEC,pvi,nvi)
        if (iexit == 1 .or. iexit == 2) then
            if (val_exit < 0.) STOP 'exit value must be positive'
            !null geodesics and DA
            allocate(yi(1:10))
            allocate(yf(1:10))
            yi(1:4)=nvi
            yi(5:8)=pvi
            yi(9)=1._dp  !dDA/ds
            yi(10)=0._dp !DA
        else
            if (val_exit > 0.) STOP 'exit value must be negative'
            !null geodesics only 
            allocate(yi(1:8))
            allocate(yf(1:8))
            yi(1:4)=nvi
            yi(5:8)=pvi
        end if
        call rayshoot(si,yi,ds,yf,iexit,val_exit)
        kti=yi(1) !k^t at observer
        ktf=yf(1) !k^t at source
        zp1=abs(ktf/kti)
        if (iexit == 1) then
            !return dL at given z
            yout=yf(10)*zp1*zp1
        else
            !return z (either at given t or dL)
            yout=zp1-1._dp 
        end if
        deallocate(yf)
        deallocate(yi)
    end subroutine rayshoot_dr
    
    subroutine rayshoot(si,yinit,ds,yf,iexit,val_exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the null geodesic equation as a first order system     !
! until synchronous time t reached (some specified point in    !
! the past)                                                    !
!                                                              !
! note                                                         !
!       y = (k^\mu, x^\mu [, dDA/ds, DA])                      !
!       dy = (dk^\mu/ds, dx^\mu/ds [, d2DA/ds2, dDA/ds])       !
!       i.e.  y(5)=t, y(6)=r, y(7)=\theta, y(8)=\phi           !
!             [y(9)=dDA/ds, y(10)=DA]                          !
!                                                              !
! where s is affine parameter                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use utils, only : dopri_stepper,controller,ydout,zbrent
        implicit none
        real(dp), intent(inout) :: si,ds
        real(dp), dimension(:), intent(in) :: yinit
        real(dp), dimension(size(yinit)), intent(out) :: yf
        integer, intent(in) :: iexit
        real(dp), intent(in) :: val_exit
        real(dp), dimension(size(yinit)) :: yi,yo,dyi,dyo
        real(dp), dimension(7,size(yinit)) :: k
        real(dp), parameter :: sacc=1.e-9_dp
        real(dp) :: dsi,so,sf,zo,zp1,dL
        real(dp) :: atol,rtol,errold,err
        logical reject
        atol=1.e-8_dp
        rtol=1.e-10_dp
        reject=.true.
        errold=1.e-4_dp
  1     yi=yinit
        call ode_null(si,yi,dyi)
        do 
            call dopri_stepper(ode_null,ds,si,so,k,yi,yo,dyi,dyo,err,atol,rtol)
            call controller(err,errold,reject,ds)
! Guard against solution getting stuck
            if (ds < 1.e-14) then
                write(*,*) 'stuck at (t,r,theta,phi) = ', yo(5:8)
                write(*,*) 'restarting do loop with lower tolerance'
                atol=atol*10._dp
                rtol=rtol*10._dp
                goto 1
            end if
            dsi=ds
            if (reject) then
                cycle !retry with new ds 
            else
                errold=err
! Exit strategy given by iexit and exit value = val_exit
                if (iexit == 1) then
                    zo=abs(yo(1)/yi(1))-1._dp
                    if (zo >= val_exit) EXIT !redshift exceeds exit value
                else if (iexit == 2) then
                    zp1=abs(yo(1)/yi(1))
                    dL=yo(10)*zp1*zp1
                    if (dL >= val_exit) EXIT !dL exceeds exit value
                else
                    if (yo(5) <= val_exit) EXIT !time precedes exit value
                end if
                si=so
                yi=yo
                dyi=dyo
            end if            
        end do
! Find s=sf corresponding to exit value
        sf=zbrent(f,si,so,sacc)
! Find soln at sf by interpolating between yi and yo
        yf=ydout(dsi,si,sf,k,yi,yo)
    contains 
        real(dp) function f(s)
            implicit none
            real(dp), intent(in) :: s
            real(dp), dimension(size(yinit)) :: yint
            real(dp) zint
            yint=ydout(dsi,si,s,k,yi,yo)
            if (iexit == 1) then
                zint=abs(yint(1)/yi(1))-1._dp
                f=zint-val_exit
            else if (iexit == 2) then
                f=yint(10)-val_exit
            else
                f=yint(5)-val_exit
            end if
        end function f
    end subroutine rayshoot

    subroutine ode_null(s,y,dy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The null geodesic equations as a first-order ODE solved    !
! simultaneously with the angular diameter distance D_A      !
! This routine takes as input: (1) the indepedent variable s !
! which, although is not explicit, is generically required   !
! by the ODE solver and (2) the dependent variable y at s    !
!                                                            !
!    s    affine parameter (independent variable)            !
!    y    (k^mu, x^mu, dDA/ds, DA); size(y)=10               !
!    dy   (dk^mu/ds, dx^mu/ds, d2DA/ds2, dDA/ds)             !
!                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real(dp), intent(in) :: s
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(size(y)), intent(out) :: dy
        real(dp), dimension(4) :: x,v
        real(dp), dimension(4,4) :: gam_ij
        real(dp), dimension(4,4,4) :: gam
        type(four_position) :: xc
        real(dp) :: gkk,rho
        integer :: i,j,k
        x=y(5:8)
        v=y(1:4)
        if (x(2) < 0.) write(*,*) 'negative r coord: ', x(2)
        xc=init_four_position(x)
        call christoffel(xc,gam)
        do k=1,4
            gam_ij=gam(k,:,:)
            gkk=dot_product(matmul(gam_ij,v),v)
            dy(k)=-gkk
        end do
        dy(5:8)=y(1:4)
! compute if solving for angular diameter distance
        if (size(y) == 10) then
            call density(xc,rho)
            dy(9)=-0.5_dp*rho*y(1)*y(1)
            dy(10)=y(9)
        end if
        call del_four_position(xc)
    end subroutine ode_null
    
    subroutine construct_ray(ti,obs_pos,ra,dec,pv,nv) !fix v3
        implicit none
        real(dp), intent(in) :: ti
        real(dp), dimension(3), intent(in) :: obs_pos
        real(dp), intent(in) :: ra,dec !radians
        real(dp), dimension(4), intent(out) :: pv,nv
        real(dp), dimension(4) :: nvo
        real(dp), dimension(3) :: dv
        real(dp), dimension(4,4) :: g,Ec
        real(dp), parameter :: tol=1.e-8_dp
        real(dp) :: alpha,psi,thetac,k2
        integer :: i,j
        pv(1)=ti
        pv(2:4)=obs_pos
        call metric(pv,g)
        call tetrad(pv,Ec)
!   nv - null vector
!   dv - directional vector in local orthonormal frame
        call ang2vec(dec,ra,dv)
! Re-orient healpix coordinates so direction of north pole 
! coincides with radial direction BV(2,:). MAKE SURE phi=PI/2
        ! thetac=pi-theta
        ! psi=thetac
        ! alpha=phi-0.5*PI
        ! call rotatex(psi,dv)
        ! call rotatez(alpha,dv)
        nvo(1)=-1._dp
        nvo(2)=dv(1)
        nvo(3)=dv(2)
        nvo(4)=dv(3)
        nv=matmul(Ec,nvo)
        ! k2=dot_product(matmul(g,nv),nv)
        ! print *, k2
! perturb if ray passes through coord centre
        if (abs(nv(3)) < tol .and. abs(nv(4)) < tol) then
            write(*,*) 'warning: ray crosses r=0'
            nv(3)=sign(tol,nv(3))
            nv(4)=sign(tol,nv(4))
! need to re-normalise it
            k2=0._dp
            do j=2,4
                do i=2,4
                    k2=k2+g(i,j)*nv(i)*nv(j)
                enddo
            enddo
            nv=nv/sqrt(abs(k2))
            nv(1)=-sqrt(abs(k2)) 
        end if
    end subroutine construct_ray

end module sznullgeo
