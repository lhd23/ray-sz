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
