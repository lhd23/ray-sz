subroutine raytrace(x_model,x_alpha,x_amp,x_r0,x_r_obs,x_theta_obs,yout)
    use cosmo_params, only : model,alpha,amp,r0,r_obs,theta_obs,FLRW
    use constants, only : PI
    use szcmb
    implicit none
    character(len=8), intent(in) :: x_model
    real(kind=8), intent(in) :: x_alpha
    real(kind=8), intent(in) :: x_amp
    real(kind=8), intent(in) :: x_r0
    real(kind=8), intent(in) :: x_r_obs
    real(kind=8), intent(in) :: x_theta_obs
    real(kind=8), intent(out) :: yout
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
    model=x_model
    alpha=x_alpha
    amp=x_amp
    r0=x_r0
    r_obs=x_r_obs/FLRW%h
    theta_obs=x_theta_obs*PI

    call cmbcal(nside)
    yout=0.
end subroutine raytrace
