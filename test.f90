program test
  use precision1, only : dp
  use constants
  use cosmo_params, only : model,alpha,amp,r0,dlr,r_obs,theta_obs,FLRW
  use szlocal
  use test_module
  implicit none

  open(unit=9,file='params')
  read(9,*) model
  read(9,*) alpha
  read(9,*) amp
  read(9,*) r0
  read(9,*) r_obs
  read(9,*) theta_obs
  close(9)

  r_obs = r_obs/FLRW%h
  theta_obs = theta_obs*PI !thta_obs is fraction of PI

  call init_model

  ! call test_get_k
  ! call test_get_kprime
  ! call test_get_Rprime
  ! call test_Rdotprime
  ! call test_get_kpprime
  ! call test_Hubble
  ! call test_pX555
  call test_propagation
  
end program test
