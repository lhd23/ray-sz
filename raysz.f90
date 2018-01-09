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
subroutine raytrace(p_alpha,p_amp,p_r0,p_r_obs,p_theta_obs,yout)
    use cosmo_params, only : model,alpha,amp,r0,r_obs,theta_obs,FLRW
    use constants, only : PI
    use szlocal, only : init_model
    use szcmb
    implicit none
!============= INPUT PARAMETERS ===============
    real(kind=8),     intent(in) :: p_alpha
    real(kind=8),     intent(in) :: p_amp
    real(kind=8),     intent(in) :: p_r0
    real(kind=8),     intent(in) :: p_r_obs
    real(kind=8),     intent(in) :: p_theta_obs
!==============================================
    real(kind=8),     intent(out) :: yout
    integer, parameter :: nside=8
    real(dp), dimension(0:12*nside*nside-1) :: dtt

    model     = 'Szekeres'
    alpha     = p_alpha
    amp       = p_amp
    r0        = p_r0
    r_obs     = p_r_obs/FLRW%h
    theta_obs = p_theta_obs*PI

! Initialise model; compute splines; calculate age of universe
    call init_model

    call cmbcal(nside,dtt,iwrite=0)
    yout=0.

end subroutine raytrace

subroutine loglike(p_alpha,p_amp,p_r0,p_r_obs,p_theta_obs,yout)
    use cosmo_params, only : model,alpha,amp,r0,r_obs,theta_obs,FLRW
    use constants, only : PI
    use szlocal, only : init_model
    use szcmb
    implicit none
!============= INPUT PARAMETERS ===============
    real(kind=8),     intent(in) :: p_alpha
    real(kind=8),     intent(in) :: p_amp
    real(kind=8),     intent(in) :: p_r0
    real(kind=8),     intent(in) :: p_r_obs
    real(kind=8),     intent(in) :: p_theta_obs
!==============================================
    real(kind=8),     intent(out) :: yout
    integer,          parameter :: nside=8

    model     = 'Szekeres'
    alpha     = p_alpha
    amp       = p_amp
    r0        = p_r0
    r_obs     = p_r_obs/FLRW%h
    theta_obs = p_theta_obs*PI

! Initialise model; compute splines; calculate age of universe
    call init_model

    yout=-0.5_dp*chi_squared(nside)

end subroutine loglike
