module constants
  use precision1, only : dp
  implicit none
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: HALFPI = 1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: PIO3 = 1.047197551196597746154214461093167628066_dp
  real(dp), parameter :: THIRD = 0.333333333333333333333333333333333333333_dp
  real(dp), parameter :: TWOTHIRD = 0.66666666666666666666666666666666666666_dp
  real(dp), parameter :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
  real(dp), parameter :: SQRT3 = 1.732050807568877293527446341505872366943_dp
  real(dp), parameter :: DEG2RAD = PI/180._dp
  real(dp), parameter :: RAD2DEG = 180._dp/PI
  
  real(dp), parameter :: MPC_IN_M = 3.085678e22_dp
  real(dp), parameter :: SOLAR_MASS_IN_KG = 1.989e30_dp
  real(dp), parameter :: YEAR_IN_SEC = 31557600._dp
  
  real(dp), parameter :: MU = 1.0e15_dp * SOLAR_MASS_IN_KG !1 Mass unit   = 10^15 M_sun
  real(dp), parameter :: LU = MPC_IN_M                     !1 Length unit = Mpc
  real(dp), parameter :: TU = 1.0e6_dp * YEAR_IN_SEC       !1 Time unit   = 10^6 years

  real(dp), parameter :: SPEED_OF_LIGHT = 299792458._dp
  real(dp), parameter :: CS = SPEED_OF_LIGHT * TU / LU
  real(dp), parameter :: GCONS = 6.6742e-11_dp * MU * TU**2 / LU**3
  real(dp), parameter :: KAPPA = 8._dp * PI * GCONS / CS**4
  real(dp), parameter :: KAPPA_C2 = KAPPA * CS**2
end module constants
