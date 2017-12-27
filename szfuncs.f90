module szfuncs 
  use precision1, only : dp
  implicit none

contains

    real(dp) function delta(r)
        use cosmo_params, only : amp,dlr,r0
        implicit none
        real(dp), intent(in) :: r
        real(dp) :: x
        x=0.5_dp*(r-r0)/dlr
        delta=0.5_dp*amp*(1._dp-tanh(x))
    end function delta

    real(dp) function deltaprime(r)
        use cosmo_params, only : amp,dlr,r0
        implicit none
        real(dp), intent(in) :: r
        real(dp) :: x,f
        x=0.5_dp*(r-r0)/dlr
        f=1._dp/cosh(x)**2
        deltaprime=-0.25_dp*amp*f/dlr
    end function deltaprime

    real(dp) function deltapprime(r)
        use cosmo_params, only : amp,dlr,r0
        implicit none
        real(dp), intent(in) :: r
        real(dp) :: x,f
        x=0.5_dp*(r-r0)/dlr
        f=sinh(x)/cosh(x)**3
        deltapprime=0.25_dp*amp*f/dlr**2
    end function deltapprime
    
    real(dp) function Mfunc(r)
        use cosmo_params, only : M0
        implicit none
        real(dp), intent(in) :: r
        Mfunc=M0*(1._dp+delta(r))*r**3
    end function Mfunc

    real(dp) function Mprime(r)
        use cosmo_params, only : M0
        implicit none
        real(dp), intent(in) :: r
        real(dp) :: delp1,ddel
        delp1=delta(r)+1._dp
        ddel=deltaprime(r)
        Mprime=M0*(3._dp*delp1+ddel*r)*r**2
    end function Mprime

    real(dp) function Mpprime(r)
        use cosmo_params, only : M0
        real(dp), intent(in) :: r
        real(dp) :: delp1,ddel,dddel,a
        delp1=delta(r)+1._dp
        ddel=deltaprime(r)
        dddel=deltapprime(r)
        a=6._dp*(delp1+ddel*r)
        Mpprime=M0*(a+dddel*r**2)*r
    end function Mpprime
    
    real(dp) function Omr0(r)
        use cosmo_params, only : FLRW
        implicit none
        real(dp), intent(in) :: r
        Omr0=FLRW%Om*(1._dp+delta(r))
    end function Omr0

    real(dp) function k_approx(r)
        use cosmo_params, only : H0,FLRW
        implicit none
        real(dp), intent(in) :: r
        k_approx=2.455*FLRW%Om*delta(r)*(H0*r)**2
    end function k_approx

    real(dp) function S_sz(r)
        use cosmo_params, only : alpha
        implicit none
        real(dp), intent(in) :: r
        S_sz=r**alpha
    end function S_sz

    real(dp) function Sp_sz(r)
        use cosmo_params, only : alpha
        implicit none
        real(dp), intent(in) :: r
        Sp_sz=alpha*S_sz(r)/r
    end function Sp_sz

    real(dp) function Spp_sz(r)
        use cosmo_params, only : alpha
        implicit none
        real(dp), intent(in) :: r
        Spp_sz=(alpha-1._dp)*Sp_sz(r)/r
    end function Spp_sz

    real(dp) function P_sz(r)
        implicit none
        real(dp), intent(in) :: r
        P_sz=1._dp
    end function P_sz

    real(dp) function Pp_sz(r)
        implicit none
        real(dp), intent(in) :: r
        Pp_sz=0._dp
    end function Pp_sz

    real(dp) function Ppp_sz(r)
        implicit none
        real(dp), intent(in) :: r
        Ppp_sz=0._dp
    end function Ppp_sz

    real(dp) function Q_sz(r)
        implicit none
        real(dp), intent(in) :: r
        Q_sz=0._dp
    end function Q_sz

    real(dp) function Qp_sz(r)
        implicit none
        real(dp), intent(in) :: r
        Qp_sz=0._dp
    end function Qp_sz

    real(dp) function Qpp_sz(r)
        implicit none
        real(dp), intent(in) :: r
        Qpp_sz=0._dp
    end function Qpp_sz

end module szfuncs
