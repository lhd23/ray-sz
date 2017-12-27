module sznullgeo
  use precision1, only : dp
  use constants
  use szgeom, only : metric,tetrad,christoffel
  use pix_tools, only : ang2vec
  implicit none

  private

  public :: propagation,ode_null
  
contains

    subroutine propagation(initial_pos,RA,DEC,redshift)
! If raytracing over the sky then phi is arbitrary for axisymmetry
        implicit none
        real(dp), intent(in), dimension(3) :: initial_pos
        real(dp), intent(in) :: RA,DEC
        real(dp), intent(out) :: redshift
        real(dp), dimension(4) :: pvi,nvi,pv,nv
        real(dp), parameter :: t_init=0._dp
        real(dp), parameter :: tf=-400._dp 
        real(dp) :: si,ds,kti,ktf,nv_(4),g(4,4),k2
        si=0._dp
        ds=1.e-2_dp
        call construct_geodesic(t_init,initial_pos,RA,DEC,pvi,nvi)
        call solve_nullgeo_eqns(si,pvi,nvi,ds,tf,pv,nv)
        kti=nvi(1) !k^t
        ktf=nv(1)
        redshift=abs(ktf/kti)-1._dp
        call metric(pv,g)
        nv_=matmul(g,nv)
        k2=dot_product(nv_,nv); print *, 'k^mu k_mu: ', k2
        return
    end subroutine propagation
    
    subroutine solve_nullgeo_eqns(si,pvi,nvi,ds,tf,pv,nv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the null geodesic equation as a first order system     !
! until synchronous time t reached (some specified point in    !
! the past)                                                    !
!                                                              !
! note:                                                        !
!       y = (k^\mu, x^\mu)                                     !
!       dy = \frac{dy}{ds}                                     !
!          = (\frac{dk^\mu}{ds}, \frac{dx^\mu}{ds},            !
!       i.e.  y(5)=t, y(6)=r, y(7)=\theta, y(8)=\phi           !
!                                                              !
! where s is affine parameter                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use utils, only : dopri_stepper,controller,ydout,zbrent
        implicit none
        real(dp), intent(in) :: tf
        real(dp), intent(inout) :: si,ds
        real(dp), dimension(4), intent(in) :: pvi,nvi
        real(dp), dimension(4), intent(out) :: pv,nv
        real(dp), dimension(8) :: yi,yo,yf,dyi,dyo
        real(dp), dimension(7,size(yi)) :: k
        real(dp), parameter :: atol=1.e-9_dp !1e-14_dp
        real(dp), parameter :: rtol=1.e-11_dp !1e-15_dp
        real(dp), parameter :: sacc=1.e-9_dp
        real(dp) :: dsi,so,sf
        real(dp) :: errold,err
        logical reject
        reject=.true.
        errold=1.e-4_dp
        yi(1:4)=nvi
        yi(5:8)=pvi
        call ode_null(si,yi,dyi)
        do 
            call dopri_stepper(ode_null,ds,si,so,k,yi,yo,dyi,dyo,err,atol,rtol)
            call controller(err,errold,reject,ds)
            dsi=ds
            if (reject) then
                cycle !retry with new ds 
            else
                errold=err
                if (yo(5) <= tf) then !yo(5)=to
                    EXIT
                else
                    si=so
                    yi=yo
                    dyi=dyo
                end if
            end if            
        end do
! Find point s=sf such that t=tf
        sf=zbrent(f,si,so,sacc)
! Interpolate between yi and yo
        yf=ydout(dsi,si,sf,k,yi,yo)
        pv=yf(5:8)
        nv=yf(1:4)
    contains 
        real(dp) function f(s)
            implicit none
            real(dp), intent(in) :: s
            real(dp), dimension(size(yi)) :: yint
            yint=ydout(dsi,si,s,k,yi,yo)
            f=yint(5)-tf
        end function f
    end subroutine solve_nullgeo_eqns

    subroutine ode_null(s,y,dy)
        implicit none
        real(dp), intent(in) :: s
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(size(y)), intent(out) :: dy
        real(dp), dimension(4) :: x,v
        real(dp), dimension(4,4) :: gamj
        real(dp), dimension(4,4,4) :: gam
        real(dp) :: gkk
        integer j
        x=y(5:8)
        v=y(1:4)
        call christoffel(x,gam)
        do j=1,4
            gamj=gam(j,:,:)
            gkk=dot_product(matmul(gamj,v),v)
            dy(j)=-gkk
        end do
        dy(5:8)=y(1:4)
    end subroutine ode_null
    
    subroutine construct_geodesic(ti,obs_pos,ra,dec,pv,nv) !fix v3
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
    end subroutine construct_geodesic

end module sznullgeo
