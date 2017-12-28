module szcmb
  use precision1, only : dp
  use constants, only : PI
  implicit none

contains

    subroutine cmbcal(nside)
!
!   Calculate CMB sky at observer
!   and return (T0-<T0>)/<T0> as a healpix map
!   
        use cosmo_params, only : r_obs,theta_obs,model
        use szlocal, only : init_model
        use sznullgeo, only : propagation
        use pix_tools, only : pix2ang_ring,nside2npix
        implicit none
        integer, intent(in) :: nside
        real(dp), dimension(:), allocatable :: dtt
        real(dp), dimension(:), allocatable :: hpt 
        real(dp) :: l,b,red,tav,obs_pos(3)
        integer :: npix,ipix,IERR,JERR
        character(len=28) :: filename

        npix=nside2npix(nside)
        allocate(hpt(0:npix-1),stat=IERR)
        allocate(dtt(0:npix-1),stat=JERR)
        if ((IERR+JERR) > 0) &
                STOP 'cmbcal: HEALPIX map allocation failed'

! Initialise model; compute splines; calculate age of universe
        call init_model
        obs_pos=(/r_obs,theta_obs,0.543*PI/)

!$omp parallel do private(red,l,b) &
!$omp shared(obs_pos,hpt) &
!$omp schedule(guided)
        do ipix=0,npix-1
            call pix2ang_ring(nside,ipix,theta=b,phi=l)
            call propagation(obs_pos,l,b,red)
            hpt(ipix) = 1._dp/(1._dp+red)
        end do
!$omp end parallel do
        
        tav=sum(hpt)/real(npix,kind=dp)
        dtt=hpt/tav-1._dp
        write(filename,'(a,a,a,i4.4,a)') &
                'map_',model,'_Nside',nside,'.txt'
        open(unit=44, file=filename, status='replace')
        do ipix=0,npix-1
            write(44,*) hpt(ipix), dtt(ipix)
        end do
        close(44)

        call cmbstats(nside,dtt)

        deallocate(dtt)
        deallocate(hpt)

    end subroutine cmbcal

    subroutine cmbstats(nside,map)
        use cosmo_params, only : T0CMB
        use alm_tools, only : map2alm,alm2cl
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
! Compute dipole and quadrupole
        call map2alm(nside,lmax,mmax,map,alm,zbounds,dw8)
        call alm2cl(lmax,mmax,alm,cl)

        map_max=maxval(map)
        map_min=minval(map)
        print *, 'C_0 =', cl(0,1)
        print *, 'C_1 =', cl(1,1)
        print *, 'C_2 =', cl(2,1)
        print *, 'Tdpl :', 1.5_dp*sqrt(cl(1,1)/PI)
        print *, 'map maxval: ', map_max
        print *, 'map minval: ', map_min
        dipole=1.5_dp*sqrt(cl(1,1)/PI)*T0CMB
        quadrupole = (3._dp/PI)*cl(2,1)*T0CMB**2

        print *, '---------------------------------------------'
        print *, 'Dipole [mK] =', dipole * 1e3_dp, &
                 'Quad [uK^2] =', quadrupole * 1e12_dp
        print *, '---------------------------------------------'

    end subroutine cmbstats

end module szcmb

!-----------------------------------------------------------------

program raytrace
    use precision1, only : dp
    use constants, only : PI
    use cosmo_params, only : model,alpha,amp,r0,r_obs,theta_obs,FLRW
    use szcmb
    use elliptic_cache, only : count,calls
    implicit none
    integer, parameter :: nside=8

!======== Foreground model parameters ===========
!
!   alpha  :   dipole parameter
!   amp    :   void amplitude  
!   r0     :   void transition radius (Mpc)
!   model  :   FLRW, LTB, or Szekeres
!
!== Observer parameters
!   r_obs, thta_obs : observer coordinates
!   Note thta_obs in units of pi i.e.
!
!           0 <= thta_obs <= 1
!
    open(unit=9,file='params')
    read(9,*) model
    read(9,*) alpha
    read(9,*) amp
    read(9,*) r0
    read(9,*) r_obs
    read(9,*) theta_obs
    close(9)
    r_obs=r_obs/FLRW%h
    theta_obs=theta_obs*PI

    call cmbcal(nside)

    print *, 'count/calls: ', count, calls

end program raytrace
