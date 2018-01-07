module szcmb
  use precision1, only : dp
  use constants
  use utils ,only : vec2ang,ang2vec,rotate
  use pix_tools, only : pix2ang_ring,nside2npix
  use alm_tools, only : map2alm,alm2cl
  implicit none

contains

    subroutine cmbcal(nside,dtt,iwrite)
!
!   Calculate CMB sky at observer
!   and return (T0-<T0>)/<T0> as a healpix map
!   
        use cosmo_params, only : r_obs,theta_obs,model
        use sznullgeo, only : rayshoot_dr
        implicit none
        integer, intent(in) :: nside
        integer, intent(in) :: iwrite
        real(dp), dimension(0:12*nside*nside-1), intent(out) :: dtt
        real(dp), dimension(0:12*nside*nside-1) :: zp1_inv !1/(z+1)
        real(dp) :: tht_ls,phi_ls,z,zp1_inv_ave,obs_pos(3)
        integer :: npix,ipix
        character(len=28) :: filename

        obs_pos=(/r_obs,theta_obs,0.543*PI/)
! sky appears the same for any azimuthal angle (phi arbitrary)

        npix=nside2npix(nside)
!$omp parallel do private(z,tht_ls,phi_ls) &
!$omp shared(obs_pos,zp1_inv) &
!$omp schedule(guided)
        do ipix=0,npix-1
            call pix2ang_ring(nside,ipix,theta=tht_ls,phi=phi_ls)
            call rayshoot_dr(obs_pos,tht_ls,phi_ls,z,iexit=0,val_exit=-400._dp)
            zp1_inv(ipix)=1._dp/(1._dp+z)
        end do
!$omp end parallel do
        
        zp1_inv_ave=sum(zp1_inv)/real(npix,kind=dp)
        dtt=zp1_inv/zp1_inv_ave-1._dp

        if (iwrite == 1) then
            write(filename,'(a,a,a,i4.4,a)') &
                    'map_',model,'_Nside',nside,'.txt'
            open(unit=44, file=filename, status='replace')
            do ipix=1,npix
                write(44,*) zp1_inv(ipix), dtt(ipix)
            end do
            close(44)
        end if
        return
    end subroutine cmbcal

    subroutine cmbstats(nside,map)
!
!   Compute dipole and quadrupole
!
        use cosmo_params, only : T0CMB
        implicit none
        real(dp), dimension(:), intent(in) :: map
        integer, intent(in) :: nside
        integer, parameter :: lmax=2,mmax=2
        real(dp), dimension(1:2*nside,1) :: dw8
        complex(dp), dimension(1,0:lmax,0:mmax) :: alm
        real(dp), dimension(0:lmax,1) :: cl
        real(dp) :: dipole,quadrupole,map_max,map_min
        real(dp) :: z,zbounds(2)
        dw8=1._dp
        z=1._dp
        zbounds=(/-z,z/) !(/-1,1/) means no mask
        call map2alm(nside,lmax,mmax,map,alm,zbounds,dw8)
        call alm2cl(lmax,mmax,alm,cl)

        map_max=maxval(map)
        map_min=minval(map)

        write(*,*) 'C_0 :', cl(0,1)
        write(*,*) 'C_1 :', cl(1,1)
        write(*,*) 'C_2 :', cl(2,1)
        write(*,*) 'DT/T dpl :', 1.5_dp*sqrt(cl(1,1)/PI)
        write(*,*) 'map maxval :', map_max
        write(*,*) 'map minval :', map_min
        dipole=1.5_dp*sqrt(cl(1,1)/PI)*T0CMB
        quadrupole=(3._dp/PI)*cl(2,1)*T0CMB**2

        write(*,*) '---------------------------------------------'
        write(*,*) 'Dipole [mK] =', dipole * 1e3_dp
        write(*,*) 'Quad [uK^2] =', quadrupole * 1e12_dp
        write(*,*) '---------------------------------------------'

    end subroutine cmbstats

    function chi_squared(nside)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine computes the chi^2 of the composite sample for the    !
! input Szekeres model:                                              !
!                                                                    !
!     chi^2 = \sum_i [(di - dL) / \sigma di]^2                       !
!                                                                    !
! where di is the observed luminosity distance, dL is the            !
! model luminosity distance and depends on the model parameters, and !
! \sigma di is the error.                                            !
! The luminosity distance as a power series is given by              !
!                                                                    !
!     dL = c/H_0  [z + 1/2 (1-q0)z^2 + ...]                          !
!                                                                    !
! where q0 is the deceleration parameter see Weinberg 2008 (1.4.9),  !
! Visser 2004 [gr-qc/0309109].                                       !
! Then the error (to leading order) is                               !
!                                                                    !
!    \sigma dL \sim \sigma (c z / H_0) = \sigma_v / H_0              !
!                                                                    !
! where cz is the velocity and \sigma_v = c \sigma z.                !
!                                                                    !
! This routine computes dL for the Szekeres model such that the      !
! computed redshift coincides with the observed redshift.            !
! The variables tht \in [0,pi] and phi \in [0,2pi] are the local sky !
! angular coordinates related to the observer's celestial sphere     !
! NOT spacetime theta,phi.                                           !
!                                                                    !
! Composite_LG.txt has 4534 lines and 6 columns corresponding to     !
!                                                                    !
!   v[km/s], dL[Mpc/h], v_pec[km/s], sigma_v[km/s], l[deg], b[deg]   !
!                                                                    !
! where v=cz                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use cosmo_params, only : r_obs,theta_obs,q0,j0,FLRW, &
                ldipole,bdipole
        use sznullgeo, only : rayshoot_dr
        implicit none
        real(dp) :: chi_squared
        integer,  intent(in) :: nside
        real(dp), dimension(0:12*nside*nside-1) :: dtt
        real(dp), dimension(:), allocatable :: chi2_arr
        real(dp), dimension(3) :: dpl_vec
        real(dp), dimension(3) :: obs_pos
        real(dp) :: vi,di,vpeci,sig_vi,li,bi,sig_di
        real(dp) :: tht_dpl,phi_dpl,tht_dpl0,phi_dpl0
        real(dp) :: tht,phi,dL,z_exit,a1,a2,b1,b2
        integer  :: i,n,io
        real(dp), parameter :: H_0=100._dp*FLRW%h !km/s/Mpc
        real(dp), parameter :: &
                VEL2RED=1._dp/(SPEED_OF_LIGHT*1.0e-3_dp)
        !q0 j0 what for?
! observed dipole angle
        tht_dpl0=HALFPI-bdipole
        phi_dpl0=ldipole

        call cmbcal(nside,dtt,iwrite=0)
        call dipole_vec(nside,dtt,dpl_vec)
        call vec2ang(dpl_vec,tht_dpl,phi_dpl)

        obs_pos=(/r_obs,theta_obs,0.543*PI/)

        open(unit=77,file='Composite_LG.txt',status='OLD',action='READ')
        n=0
        do
            read(77,*,iostat=io) vi,di,vpeci,sig_vi,li,bi
            if (io /= 0) EXIT
            n=n+1
        end do
        rewind(77)
        allocate(chi2_arr(1:n))

        a1=phi_dpl0
        a2=HALFPI-tht_dpl0
        b1=HALFPI-tht_dpl
        b2=phi_dpl

!$omp parallel do private(vi,di,vpeci,sig_vi,li,bi,sig_di,tht,phi,dL,z_exit) &
!$omp shared(obs_pos,chi2_arr,a1,a2,b1,b2)
        do i=1,n
            read(77,*) vi,di,vpeci,sig_vi,li,bi
            di=di/FLRW%h !evaluate h to get in units Mpc
            sig_di=sig_vi/H_0
            z_exit=vi*VEL2RED
            call galactic2szekeres(li,bi,tht,phi,a1,a2,b1,b2)
            call rayshoot_dr(obs_pos,tht,phi,dL,iexit=1,val_exit=z_exit)
            ! print *, z_exit,di,dL
            chi2_arr(i)=((di-dL)/sig_di)**2
        end do
!$omp end parallel do

        chi_squared=sum(chi2_arr)

        close(77)
        deallocate(chi2_arr)
        return
    end function chi_squared

    subroutine dipole_vec(nside,map,dpl_vec)
        implicit none
        integer, intent(in) :: nside
        real(dp), dimension(:), intent(in) :: map
        real(dp), dimension(3), intent(out) :: dpl_vec
        real(dp), dimension(1:2*nside,1) :: dw8
        real(dp), dimension(2) :: zbounds=(/-1._dp,1._dp/)
        integer, parameter :: lmax=2,mmax=2
        complex(dp), dimension(1,0:lmax,0:mmax) :: alm
        real(dp) :: a10,a11r,a11i,norm
        dw8=1._dp
        call map2alm(nside,lmax,mmax,map,alm,zbounds,dw8)
        a10=real(alm(1,1,0))
        a11r=real(alm(1,1,1))
        a11i=aimag(alm(1,1,1))
        dpl_vec(1)=-SQRT2*a11r
        dpl_vec(2)=+SQRT2*a11i
        dpl_vec(3)=+a10
        norm=sqrt(a10*a10+2._dp*(a11r**2+a11i**2))
        dpl_vec=dpl_vec/norm
    end subroutine dipole_vec


    subroutine galactic2szekeres(l,b,theta,phi,a1,a2,b1,b2) !formerly invrotmap
!
!   Maps galactic coordinates (l,b) [deg] to szekeres local
!   sky coordinates (theta,phi) (NOT spacetime theta,phi)
!
        implicit none
        real(dp), intent(in) :: l,b
        real(dp), intent(out) :: theta,phi
        real(dp), intent(in) :: a1,a2,b1,b2
        real(dp), dimension(3) :: vector
        theta=HALFPI-b*DEG2RAD
        phi=l*DEG2RAD
        call ang2vec(theta,phi,vector)
! rotate dipole vector onto ex=(1,0,0)
        call rotate(vector,psi=-a1,axis='z')
        call rotate(vector,psi=-a2,axis='y')
! rotate ex onto the found dipole in szekeres
        call rotate(vector,psi=+b1,axis='y')
        call rotate(vector,psi=+b2,axis='z')
        call vec2ang(vector,theta,phi)
    end subroutine galactic2szekeres

end module szcmb
