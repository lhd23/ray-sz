module cosmo_params
  use precision1, only : dp
  use constants
  implicit none

  type FLRW_pars
      real(dp) :: h
      real(dp) :: Om
      real(dp) :: Ok
      real(dp) :: Ol
  end type FLRW_pars

! Planck 2013 parameters
  real(dp), parameter :: h_plk  = 0.673_dp
  real(dp), parameter :: Om_plk = 0.315_dp
  real(dp), parameter :: Ok_plk = 0._dp
  real(dp), parameter :: Ol_plk = 0.685_dp

! Background model parameters
  type(FLRW_pars), parameter :: &
          FLRW=FLRW_pars(h_plk,Om_plk,Ok_plk,Ol_plk)

  real(dp), parameter ::                        &
       H0    = 1.e5_dp * FLRW%h * TU / LU / CS, & !1/Mpc
       rho_c = 3._dp * H0**2 / KAPPA_C2,        &
       lb    = 3._dp * FLRW%Ol * H0**2,         &
       gkr   = KAPPA_C2 * rho_c * FLRW%Om,      &
       M0    = 0.5_dp * FLRW%Om * H0**2,        &
       T0CMB = 2.7255_dp,                       &
       q0    = -0.55_dp,                        &
       j0    = 1.0_dp

  real(dp), parameter :: ldipole  = 276.4_dp
  real(dp), parameter :: bdipole  = 29.3_dp

  real(dp), pointer :: ct2 => null() !age of universe in Mpc

! Model specific parameters
  character(len=8), save :: model
  real(dp), save :: alpha,amp,r0,dlr

! Observer position
  real(dp), save :: r_obs,theta_obs

contains

  subroutine get_age
    implicit none
    real(dp) :: q
    if (.not. associated(ct2)) then
        allocate(ct2)
        if (abs(FLRW%Ok) < 1.e-8_dp) then
            q=sqrt(abs(FLRW%Ol))
            ct2=log((1._dp+q)/(1._dp-q))/(3._dp*q*H0)
        else
            STOP 'Hardcoded to work for Ok=0 only!'
        end if
    else
        return
    end if
  end subroutine get_age
  
end module cosmo_params

